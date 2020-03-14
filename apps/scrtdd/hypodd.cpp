/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU Affero General Public License as published*
 * by the Free Software Foundation, either version 3 of the License, or    *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This program is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#include "hypodd.h"

#include <seiscomp3/core/datetime.h>
#include <seiscomp3/core/strings.h>
#include <seiscomp3/core/typedarray.h>
#include <seiscomp3/io/recordinput.h>
#include <seiscomp3/io/records/mseedrecord.h>
#include <seiscomp3/processing/operator/transformation.h>
#include <seiscomp3/processing/operator/ncomps.h>
#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/filter.h>
#include <seiscomp3/client/inventory.h>
#include <seiscomp3/datamodel/station.h>
#include <seiscomp3/utils/files.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <set>
#include <regex>
#include <boost/filesystem.hpp>
#include <boost/bind.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <sys/wait.h>
#include <unistd.h>

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp3/logging/log.h>
 
using namespace std;
using namespace Seiscomp;
using namespace Seiscomp::Processing;
using Seiscomp::Core::stringify;
using DataModel::ThreeComponents;

namespace {

pid_t startExternalProcess(const vector<string> &cmdparams,
                           bool waitChild,
                           const string& workingDir="")
{
    pid_t pid;
    string cmdline;
    vector<char *> params(cmdparams.size());
    for ( size_t i = 0; i < cmdparams.size(); ++i )
    {
        params[i] = (char*)cmdparams[i].c_str();
        if ( i > 0 ) cmdline += " ";
        cmdline += cmdparams[i];
    }
    params.push_back(nullptr);

    if ( ! workingDir.empty() )
        SEISCOMP_INFO("Working directory %s", workingDir.c_str());
    SEISCOMP_INFO("Executing command: %s ", cmdline.c_str());

    pid = fork();

    if ( pid < 0 ) // fork error
    {
        SEISCOMP_ERROR("Error (%d) in fork()", pid);
        return pid;
    }

    if ( pid == 0 ) // child
    {
        if ( ! workingDir.empty() )
        {
            if( chdir(workingDir.c_str()) != 0 )
            {
                exit(1);
            }
        }

        execv(params[0], &params[0]);
        exit(1);
    }
    else // parent
    {
        if (waitChild) // wait for the child to complete
        {
            int   status;
            do {
                pid = waitpid(pid, &status, 0);
            } while (pid == -1 && errno == EINTR);

            if (status != 0)
                SEISCOMP_ERROR("Command exited with non zero value (%d)", status);
        }
    }

    return pid;
}



void copyFileAndReplaceLines(const string& srcFilename,
                             const string& destFilename,
                             map<int,string> linesToReplace,
                             const string& comment="*")
{
    ifstream srcFile(srcFilename);
    ofstream destFile(destFilename);
    if ( ! srcFile.is_open() || ! destFile.is_open() )
    {
        string msg = stringify("Cannot copy %s to %s", srcFilename.c_str(), destFilename.c_str());
        throw runtime_error(msg);
    }

    string line;
    int lineNum = 0;
    while( std::getline(srcFile, line) )
    {
        // increase line number when not a comment
        if ( line.rfind(comment, 0) != 0 )
            lineNum++;

        // replace line
        if ( linesToReplace.find(lineNum) != linesToReplace.end())
        {
            line = linesToReplace[lineNum];
            linesToReplace.erase(lineNum);
        }

        // copy line to output
        destFile << line << std::endl;
    }
}


/*
 * Compute distance in km between two points
 */
double computeDistance(double lat1, double lon1, double depth1,
                       double lat2, double lon2, double depth2,
                       double *azimuth = nullptr, double *backAzimuth = nullptr)
{
    double distance, az, baz;
    Math::Geo::delazi(lat1, lon1, lat2, lon2, &distance, &az, &baz);

    if (azimuth) *azimuth = az;
    if (backAzimuth) *backAzimuth = baz;

    double Hdist = Math::Geo::deg2km(distance);
    double Vdist = abs(depth1 - depth2);
    // this is an approximation that works when the distance is small
    // and the Earth curvature can be assumed flat
    return std::sqrt( std::pow(Hdist,2) + std::pow(Vdist,2) );
}


double computeDistance(const HDD::Catalog::Event& ev1, const HDD::Catalog::Event& ev2,
                       double *azimuth = nullptr, double *backAzimuth = nullptr)
{
    return computeDistance(ev1.latitude, ev1.longitude, ev1.depth,
                           ev2.latitude, ev2.longitude, ev2.depth,
                           azimuth, backAzimuth);
}


double computeDistance(const HDD::Catalog::Event& event, const HDD::Catalog::Station& station,
                       double *azimuth = nullptr, double *backAzimuth = nullptr)
{
    return computeDistance(event.latitude, event.longitude, event.depth,
                           station.latitude, station.longitude, -(station.elevation/1000.),
                           azimuth, backAzimuth);
}


/*
 * Ellipsoid standard equation:
 *
 *      (x-xo)^2 / axix_a^2 + (y-yo)^2 / axis_b^2 + (z-zo)^2 / axis_c^2 = 1
 *
 */
struct Ellipsoid
{
    bool isInside(double lat, double lon, double depth) const
    {
        double distance, az, baz;
        Math::Geo::delazi(lat, lon, this->lat, this->lon, &distance, &az, &baz);

        distance = Math::Geo::deg2km(distance); // distance to km
        az += orientation; // correct azimuth by the orientation of a and b axes

        double dist_x = distance * std::cos(az);
        double dist_y = distance * std::sin(az);
        double dist_z = depth - this->depth;

        double one = std::pow( dist_x / axis_a, 2) +
                     std::pow( dist_y / axis_b, 2) +
                     std::pow( dist_z / axis_c, 2);
        return one <= 1;
    }

    double axis_a=0, axis_b=0, axis_c=0; // axis in km
    double lat=0, lon=0, depth=0;        // origin
    double orientation = 0;              // degress: when 0 -> axis_a is East-West and axis_b is North-South
};


/*
 * Helper class to implement Waldhauser's paper method of neighboring events selection
 * based on 5 concentric ellipsoidal layers
 *
 * Quadrants (1-4 above depth, 5-8 below depth):
 *
 *        lat
 *         ^
 *         |
 *    2/6  |   1/5
 *         |
 * -----------------> lon
 *         |
 *    3/7  |   4/8
 *         |
 */
class HddEllipsoid : public Core::BaseObject
{

public:

    HddEllipsoid(double vertical_axis_len, double lat, double lon, double depth)
        : HddEllipsoid(vertical_axis_len/2., vertical_axis_len, lat, lon, depth)
    { }

    HddEllipsoid(double inner_vertical_axis_len, double outer_vertical_axis_len,
                 double lat, double lon, double depth)
    {
        _outerEllipsoid.axis_a = outer_vertical_axis_len / 2.;
        _outerEllipsoid.axis_b = _outerEllipsoid.axis_a;
        _outerEllipsoid.axis_c = outer_vertical_axis_len;

        _innerEllipsoid.axis_a = inner_vertical_axis_len / 2.;
        _innerEllipsoid.axis_b = _innerEllipsoid.axis_a;
        _innerEllipsoid.axis_c = inner_vertical_axis_len;

        _outerEllipsoid.lat   = _innerEllipsoid.lat   = lat;
        _outerEllipsoid.lon   = _innerEllipsoid.lon   = lon;
        _outerEllipsoid.depth = _innerEllipsoid.depth = depth;
    }

    const Ellipsoid& getInnerEllipsoid() const { return _innerEllipsoid; }
    const Ellipsoid& getOuterEllipsoid() const { return _outerEllipsoid; }

    bool isInside(double lat, double lon, double depth, int quadrant/* 1-8 */) const
    {
        // be in the right quadrant and inside the outer layer and outside of the inner one
        return isInQuadrant(_innerEllipsoid, lat, lon, depth, quadrant) &&  isInside(lat, lon, depth);
    }

    bool isInside(double lat, double lon, double depth) const
    {
        // be inside the outer layer and outside of the inner one
        return _outerEllipsoid.isInside(lat, lon, depth) && ! _innerEllipsoid.isInside(lat, lon, depth);
    }

    static bool
    isInQuadrant(const Ellipsoid& ellip, double lat, double lon, double depth, int quadrant /* 1-8 */)
    {
        if (quadrant < 1 || quadrant > 8)
            throw std::invalid_argument( "quadrant should be between 1 and 8");

        if (depth < ellip.depth && set<int>{1,2,3,4}.count(quadrant) != 0) return false;
        if (depth > ellip.depth && set<int>{5,6,7,8}.count(quadrant) != 0) return false;

        if (lon < ellip.lon && set<int>{1,4,5,8}.count(quadrant) != 0) return false;
        if (lon > ellip.lon && set<int>{2,3,6,7}.count(quadrant) != 0) return false;

        if (lat < ellip.lat && set<int>{1,2,5,6}.count(quadrant) != 0) return false;
        if (lat > ellip.lat && set<int>{3,4,7,8}.count(quadrant) != 0) return false;

        return true;
    }

private:

    Ellipsoid _outerEllipsoid;
    Ellipsoid _innerEllipsoid;
};

DEFINE_SMARTPOINTER(HddEllipsoid);



class Randomer {

public:

    Randomer(size_t min, size_t max, unsigned int seed = std::random_device{}())
        : gen_{seed}, dist_{min, max}
    { }

    // if you want predictable numbers
    void setSeed(unsigned int seed)
    {
        gen_.seed(seed);
    }

    size_t next()
    {
        return dist_(gen_);
    }

private:

    // random seed by default
    std::mt19937 gen_;
    std::uniform_int_distribution<size_t> dist_;
};

}


namespace Seiscomp {
namespace HDD {


HypoDD::HypoDD(const CatalogCPtr& catalog, const Config& cfg, const string& workingDir)
       : _workingDir(workingDir), _cfg(cfg)
{
    setCatalog(catalog);

    if ( ! Util::pathExists(_workingDir) )
    {
        if ( ! Util::createPath(_workingDir) )
        {
            string msg = "Unable to create working directory: " + _workingDir;
            throw runtime_error(msg);
        }
    }

    _cacheDir = (boost::filesystem::path(_workingDir)/"wfcache").string();
    if ( ! Util::pathExists(_cacheDir) )
    {
        if ( ! Util::createPath(_cacheDir) )
        {
            string msg = "Unable to create cache directory: " + _cacheDir;
            throw runtime_error(msg);
        }
    }

    _wfDebugDir = (boost::filesystem::path(_workingDir)/"wfdebug").string();
    if ( _cfg.wfFilter.dump )
    {
        if ( ! Util::pathExists(_wfDebugDir) )
        {
            if ( ! Util::createPath(_wfDebugDir) )
            {
                string msg = "Unable to create waveform debug directory: " + _wfDebugDir;
                throw runtime_error(msg);
            }
        }
    }
    
    _wf = new WfMngr(_cfg.step2Clustering.recordStreamURL, _cacheDir, _wfDebugDir);
    _wf->setProcessing(_cfg.wfFilter.filterStr, _cfg.wfFilter.resampleFreq);
    _wf->setSnr(_cfg.snr.minSnr, _cfg.snr.noiseStart, _cfg.snr.noiseEnd, _cfg.snr.signalStart, _cfg.snr.signalEnd);
    _wf->setWfebug(_cfg.wfFilter.dump);
}



HypoDD::~HypoDD()
{
    //
    // delete all in working directory except the cache directory
    //
    if ( _workingDirCleanup )
    {
        for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(_workingDir), {}))
        {
            if ( ! boost::filesystem::equivalent(entry, _cacheDir) && 
                 ! boost::filesystem::equivalent(entry, _wfDebugDir )  )
            {
                SEISCOMP_INFO("Deleting %s", entry.path().string().c_str());
                try { boost::filesystem::remove_all(entry); } catch ( ... ) { }
            }
        }
    } 
}



void HypoDD::setCatalog(const CatalogCPtr& catalog)
{
    _srcCat = catalog;
    _ddbgc = filterPhasesAndSetWeights(_srcCat, Catalog::Phase::Source::CATALOG,
                                       _cfg.validPphases, _cfg.validSphases);
}



// Creates dir name from event. This id has the following format:
// OriginTime_Lat_Lon_CreationDate_Random
// eg 20111210115715_46343_007519_20111210115740_6666
string HypoDD::generateWorkingSubDir(const Catalog::Event& ev) const
{
    static Randomer ran(0, 1000);
    string id = stringify("%s_%05d_%06d_%s_%04zu",
                          ev.time.toString("%Y%m%d%H%M%S").c_str(), // origin time
                          int(ev.latitude*1000), // Latitude
                          int(ev.longitude*1000), // Longitude 
                          Core::Time::GMT().toString("%Y%m%d%H%M%S").c_str(), // creation time
                          ran.next() // random number
                          );
    return id;
}



void HypoDD::preloadData()
{
    _wf->resetCounters();

    unsigned numPhases = 0, numSPhases = 0;

    //
    // Preload waveforms on disk and cache them in memory (pre-processed)
    //
    for (const auto& kv : _ddbgc->getEvents() )
    {
        const Catalog::Event& event = kv.second;
        auto eqlrng = _ddbgc->getPhases().equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
            const Catalog::Phase& phase = it->second;
            Core::TimeWindow tw = xcorrTimeWindowLong(phase);
            const auto xcorrCfg = _cfg.xcorr.at(phase.procInfo.type);

            for (string component : xcorrCfg.components )
            {
                Catalog::Phase tmpPh = phase;
                tmpPh.channelCode = WfMngr::getBandAndInstrumentCodes(tmpPh.channelCode) + component;
                _wf->getWaveform(tw, event, tmpPh, &_wfCache, _useCatalogDiskCache, true); 
            }

            numPhases++;
            if (  phase.procInfo.type == "S" ) numSPhases++;
        }
    }

    unsigned snr_low, wf_no_avail, wf_cached, wf_downloaded;
    _wf->getCounters(snr_low, wf_no_avail, wf_cached, wf_downloaded);

    SEISCOMP_INFO("Finished preloading catalog waveform data: total phases %u (P %.f%%, S %.f%%) "
                  "phases with Signal to Noise ratio too low %u (%.f%%), "
                  "phases data not available %u (%.f%%), "
                  "(waveforms downloaded %u, waveforms loaded from disk cache %u)",
                  numPhases, ((numPhases-numSPhases)* 100. / numPhases), 
                  (numSPhases* 100. / numPhases), snr_low, (snr_low * 100. / numPhases),
                  wf_no_avail, (wf_no_avail * 100. / numPhases), wf_downloaded, wf_cached);
}



/*
 * Fixed weighting scheme based on pick time uncertainties
 * Class 0: 0     - 0.025  sec
 *       1: 0.025 - 0.050  sec
 *       2: 0.050 - 0.100  sec
 *       3: 0.100 - 0.200  sec
 *       4: 0.200 - 0.400  sec
 *       5: 0.400 -        sec
 *  weight = 1 / 2^class
 */
double HypoDD::computePickWeight(double uncertainty /* secs */ ) const
{
    double weight = 0;

    if      ( uncertainty >= 0.000 && uncertainty <= 0.025 ) weight = 1.00;
    else if ( uncertainty >  0.025 && uncertainty <= 0.050 ) weight = 0.80;
    else if ( uncertainty >  0.050 && uncertainty <= 0.100 ) weight = 0.60;
    else if ( uncertainty >  0.100 && uncertainty <= 0.200 ) weight = 0.40;
    else if ( uncertainty >  0.200 && uncertainty <= 0.400 ) weight = 0.20;
    else                                                    weight = 0.10;

    return weight;
}



double HypoDD::computePickWeight(const Catalog::Phase& phase) const
{
    return computePickWeight( (phase.lowerUncertainty + phase.upperUncertainty) / 2. );
}


/*
 * Build a catalog with requested phases only and for the same event/station pair
 * make sure to have only one phase. If multiple phases are found, keep the highest
 * priority one
 */
CatalogPtr HypoDD::filterPhasesAndSetWeights(const CatalogCPtr& catalog, const Catalog::Phase::Source& source,
                                             const std::vector<std::string>& PphaseToKeep,
                                             const std::vector<std::string>& SphaseToKeep) const
{
    SEISCOMP_INFO("Selecting preferred phases from catalog");

    multimap<unsigned,Catalog::Phase> filteredS;
    multimap<unsigned,Catalog::Phase> filteredP;

    // loop through each event
    for (const auto& kv :  catalog->getEvents() )
    {
        const Catalog::Event& event = kv.second;

        // loop through all phases of current event
        const auto& eqlrng = catalog->getPhases().equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
            // keep the phase only if it has a higher priority of an existing one
            // or if this is the only one for a given station
            const Catalog::Phase& phase = it->second;

            // P phase
            auto itpp = find(PphaseToKeep.begin(), PphaseToKeep.end(), phase.type);
            if ( itpp != PphaseToKeep.end() )
            {
                auto priority = std::distance(PphaseToKeep.begin(), itpp);

                // fetch already selected P phases for current event, and
                // check if there is already a P phase for the same station
                bool inserted = false;
                auto eqlrng2 = filteredP.equal_range(event.id);
                for (auto it2 = eqlrng2.first; it2 != eqlrng2.second; ++it2)
                {
                    Catalog::Phase& existingPhase = it2->second;
                    auto existingPriority = std::distance(
                        PphaseToKeep.begin(),
                        find(PphaseToKeep.begin(), PphaseToKeep.end(), existingPhase.type)
                    );
                    if ( existingPhase.type == phase.type           &&
                         existingPhase.stationId == phase.stationId &&
                         existingPriority < priority )
                    {
                        existingPhase = phase;
                        inserted = true;
                        break;
                    }
                }
                if ( ! inserted )
                    filteredP.emplace(phase.eventId, phase);
                continue;
            }

            // S phase
            auto itsp = find(SphaseToKeep.begin(), SphaseToKeep.end(), phase.type);
            if ( itsp != SphaseToKeep.end() )
            {
                auto priority = std::distance(SphaseToKeep.begin(), itsp);

                // fetch already selected S phases for current event, and
                // check if there is already a S phase for the same station
                bool inserted = false;
                auto eqlrng2 = filteredS.equal_range(event.id);
                for (auto it2 = eqlrng2.first; it2 != eqlrng2.second; ++it2)
                {
                    Catalog::Phase& existingPhase = it2->second;
                    auto existingPriority = std::distance(
                        SphaseToKeep.begin(),
                        find(SphaseToKeep.begin(), SphaseToKeep.end(), existingPhase.type)
                    );
                    if ( existingPhase.type == phase.type           &&
                         existingPhase.stationId == phase.stationId &&
                         existingPriority < priority )
                    {
                        existingPhase = phase;
                        inserted = true;
                        break;
                    }
                }
                if ( ! inserted )
                    filteredS.emplace(phase.eventId, phase);
                continue;
            }

            SEISCOMP_DEBUG("Discard phase (%s), the type is not among the selected ones",
                           string(phase).c_str());
        }
    }

    // loop through selected phases and replace actual phase name
    //  with a generic P or S
    multimap<unsigned,Catalog::Phase> filteredPhases;
    for (auto& it : filteredP)
    {
        Catalog::Phase& phase = it.second;
        phase.procInfo.weight = computePickWeight(phase);
        phase.procInfo.type = "P";
        phase.procInfo.source = source;

        filteredPhases.emplace(phase.eventId, phase);
    }
    for (auto& it : filteredS)
    {
        Catalog::Phase& phase = it.second;
        phase.procInfo.weight = computePickWeight(phase);
        phase.procInfo.type = "S";
        phase.procInfo.source = source;

        filteredPhases.emplace(phase.eventId, phase);
    }

    return new Catalog(catalog->getStations(), catalog->getEvents(), filteredPhases);
}


CatalogPtr HypoDD::relocateCatalog(bool force, bool usePh2dt)
{
    SEISCOMP_INFO("Starting HypoDD relocator in multiple events mode");

    CatalogPtr catToReloc = new Catalog(*_ddbgc);

    // Create working directory 
    string catalogWorkingDir = (boost::filesystem::path(_workingDir)/"catalog").string(); 
    if ( !Util::pathExists(catalogWorkingDir) )
    {
        if ( !Util::createPath(catalogWorkingDir) )
        {
            string msg = "Unable to create working directory: " + catalogWorkingDir;
            throw runtime_error(msg);
        }
    }

    // write catalog for debugging purpose
    if ( ! _workingDirCleanup )
    {
        catToReloc->writeToFile(
            (boost::filesystem::path(catalogWorkingDir)/"starting-event.csv").string(),
            (boost::filesystem::path(catalogWorkingDir)/"starting-phase.csv").string(),
            (boost::filesystem::path(catalogWorkingDir)/"starting-station.csv").string() );
    }

    // Create station.dat for hypodd (if not already generated)
    string stationFile = (boost::filesystem::path(catalogWorkingDir)/"station.dat").string();
    if ( force || ! Util::fileExists(stationFile) )
    {
        createStationDatFile(catToReloc, stationFile);
    }

    string eventFile = (boost::filesystem::path(catalogWorkingDir)/"event.dat").string();
    string dtctFile = (boost::filesystem::path(catalogWorkingDir)/"dt.ct").string();
    string dtccFile = (boost::filesystem::path(catalogWorkingDir)/"dt.cc").string(); 

    set<unsigned> selectedEvents;

    if ( ! usePh2dt )
    {
        // Create event.dat for hypodd (if not already generated)
        if ( force || ! Util::fileExists(eventFile) )
        {
            createEventDatFile(catToReloc, eventFile);
        }

        // Find Neighbouring Events in the catalog
        map<unsigned,CatalogPtr> neighbourCats = selectNeighbouringEventsCatalog(
            catToReloc, _cfg.step2Clustering.minWeight,
            _cfg.step2Clustering.minESdist, _cfg.step2Clustering.maxESdist,
            _cfg.step2Clustering.minEStoIEratio, _cfg.step2Clustering.minDTperEvt,
            _cfg.step2Clustering.maxDTperEvt, _cfg.step2Clustering.minNumNeigh,
            _cfg.step2Clustering.maxNumNeigh, _cfg.step2Clustering.numEllipsoids,
            _cfg.step2Clustering.maxEllipsoidSize, true
        );

        // calculate cross correlated differential travel times
        // Create dt.cc (if not already generated)
        if ( force || ! Util::fileExists(dtccFile) )
        {
            // Perform cross correlation, which also detects picks around theoretical
            // arrival times. The catalog will be updated with those theoretical phases 
            const XCorrCache xcorr = buildXCorrCache(neighbourCats, _cfg.artificialPhases.enable);
            createDtCc(neighbourCats, dtccFile, xcorr);
        }

        // calculate absolute travel times from catalog phases
        // Create dt.ct (if not already generated)
        if ( force || ! Util::fileExists(dtctFile) )
        {
            createDtCt(neighbourCats, dtctFile);
        }

        //
        // update selected event list, numNeighbours information and artifical phases
        //
        multimap<unsigned,Catalog::Phase> newPhases;

        for (const auto& kv : neighbourCats)
        {
            unsigned evId = kv.first;
            const CatalogPtr& neighbourCat = kv.second;

            // update selected events
            selectedEvents.insert(evId);

            // update event number of neighbouring information
            const Catalog::Event& ev1 = neighbourCat->getEvents().at(evId);
            Catalog::Event ev2 = catToReloc->getEvents().at(evId);
            ev2.relocInfo.numNeighbours = ev1.relocInfo.numNeighbours;
            catToReloc->updateEvent(ev2);

            // Save new phases (artificial phases too)
            auto eqlrng = neighbourCat->getPhases().equal_range(evId);
            newPhases.insert(eqlrng.first, eqlrng.second);
        }

        // Replace phases with new ones
        catToReloc = new Catalog(catToReloc->getStations(), catToReloc->getEvents(), newPhases);
    }
    else
    {
        // Create phase.dat for ph2dt (if not already generated)
        string phaseFile = (boost::filesystem::path(catalogWorkingDir)/"phase.dat").string();
        if ( force || ! Util::fileExists(phaseFile) )
        {
            createPhaseDatFile(catToReloc, phaseFile);
        }

        // run ph2dt
        // input files: ph2dt.inp station.dat phase.dat
        // output files: station.sel event.sel event.dat dt.ct
        if ( force || !Util::fileExists(dtctFile) )
        {
            runPh2dt(catalogWorkingDir, stationFile, phaseFile);
            string stationSelFile = (boost::filesystem::path(catalogWorkingDir)/"station.sel").string();
            if ( Util::fileExists(stationSelFile) )
                boost::filesystem::copy_file(stationSelFile, stationFile, boost::filesystem::copy_option::overwrite_if_exists);
            string eventSelfile = (boost::filesystem::path(catalogWorkingDir)/"event.sel").string();
            if ( Util::fileExists(eventSelfile) )
                boost::filesystem::copy_file(eventSelfile, eventFile, boost::filesystem::copy_option::overwrite_if_exists);
        }

        // Reads the event pairs matched in dt.ct which are selected by ph2dt and
        // calculate cross correlated differential travel_times for every pair.
        // input dt.ct
        // output dt.cc
        if ( force || ! Util::fileExists(dtccFile) )
        {
            selectedEvents = createDtCcPh2dt(catToReloc, dtctFile, dtccFile);
        }
    }

    // run hypodd
    // input : dt.cc dt.ct event.sel station.sel hypoDD.inp
    // output : hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
    string ddrelocFile = (boost::filesystem::path(catalogWorkingDir)/"hypoDD.reloc").string();
    string ddresidualFile = (boost::filesystem::path(catalogWorkingDir)/"hypoDD.res").string();
    if ( force || ! Util::fileExists(ddrelocFile) || ! Util::fileExists(ddresidualFile) )
    {
        boost::filesystem::remove(ddrelocFile);
        boost::filesystem::remove(ddresidualFile);
        runHypodd(catalogWorkingDir, dtccFile, dtctFile, eventFile, stationFile, _cfg.hypodd.step2CtrlFile);
    }

    // load a catalog from hypodd output file
    // input: hypoDD.reloc
    CatalogPtr relocatedCatalog = loadRelocatedCatalog(catToReloc, ddrelocFile, ddresidualFile);

    // write catalog for debugging purpose
    if ( ! _workingDirCleanup )
    {
        relocatedCatalog->writeToFile(
            (boost::filesystem::path(catalogWorkingDir)/"relocated-event.csv").string(),
            (boost::filesystem::path(catalogWorkingDir)/"relocated-phase.csv").string(),
            (boost::filesystem::path(catalogWorkingDir)/"relocated-station.csv").string());
    }

    if ( _workingDirCleanup ) boost::filesystem::remove_all(catalogWorkingDir);

    // Remove not relocated events or events that were selected only as neighbour
    vector<unsigned> evIdsToRemove;
    for (const auto& kv : relocatedCatalog->getEvents() )
    {
        const Catalog::Event& ev = kv.second;
        if ( ! ev.relocInfo.isRelocated || 
             ( selectedEvents.count(ev.id) == 0 && selectedEvents.size() > 0 ) )
        {
            evIdsToRemove.push_back(ev.id);
        }
    }
    for (unsigned evId : evIdsToRemove ) relocatedCatalog->removeEvent(evId);

    return relocatedCatalog;
}



CatalogPtr HypoDD::relocateSingleEvent(const CatalogCPtr& singleEvent)
{
    // there must be only one event in the catalog, the origin to relocate
    Catalog::Event evToRelocate = singleEvent->getEvents().begin()->second;
    const auto& evToRelocatePhases = singleEvent->getPhases().equal_range(evToRelocate.id);

    SEISCOMP_INFO("Starting HypoDD relocator in single event mode: event %s (%ld phases)",
                  string(evToRelocate).c_str(), std::distance(evToRelocatePhases.first, evToRelocatePhases.second));

    // Create working directory
    string subFolder = generateWorkingSubDir(evToRelocate);
    subFolder = (boost::filesystem::path(_workingDir)/subFolder).string();
    if ( Util::pathExists(subFolder) )
    {
        boost::filesystem::remove_all(subFolder);
    }

    //
    // Step 1: refine location without cross correlation
    //
    SEISCOMP_INFO("Performing step 1: initial location refinement (no cross correlation)");

    string eventWorkingDir = (boost::filesystem::path(subFolder)/"step1").string();

    CatalogPtr evToRelocateCat = filterPhasesAndSetWeights(singleEvent, Catalog::Phase::Source::RT_EVENT,
                                                            _cfg.validPphases, _cfg.validSphases);

    CatalogPtr relocatedEvCat = relocateEventSingleStep(
            evToRelocateCat, eventWorkingDir, false, false, _cfg.hypodd.step1CtrlFile, _cfg.step1Clustering.minWeight,
            _cfg.step1Clustering.minESdist, _cfg.step1Clustering.maxESdist, _cfg.step1Clustering.minEStoIEratio,
            _cfg.step1Clustering.minDTperEvt, _cfg.step1Clustering.maxDTperEvt, _cfg.step1Clustering.minNumNeigh,
            _cfg.step1Clustering.maxNumNeigh, _cfg.step1Clustering.numEllipsoids, _cfg.step1Clustering.maxEllipsoidSize
    );

    if ( relocatedEvCat )
    {
        SEISCOMP_INFO("Step 1 relocation successful");
        SEISCOMP_INFO("%s", relocationReport(relocatedEvCat).c_str() );

        evToRelocateCat = relocatedEvCat;
    }
    else
    {
        SEISCOMP_ERROR("Failed to perform step 1 origin relocation");

        if (  _cfg.artificialPhases.enable )
        {
            //
            // Alternative Step 1: since step 1 failed we don't have a refined location and thus we won't
            // get good theoretical picks estimate. So let's try to refine location using cross correlation
            // before trying to create artificial phases
            //
            SEISCOMP_INFO("Performing alternative step 1: initial location refinement (with cross correlation)");

            eventWorkingDir = (boost::filesystem::path(subFolder)/"step1-xcorr").string();

            relocatedEvCat  = relocateEventSingleStep(
                    evToRelocateCat, eventWorkingDir, true, false,  _cfg.hypodd.step2CtrlFile, _cfg.step2Clustering.minWeight,
                    _cfg.step2Clustering.minESdist, _cfg.step2Clustering.maxESdist, _cfg.step2Clustering.minEStoIEratio,
                    _cfg.step2Clustering.minDTperEvt, _cfg.step2Clustering.maxDTperEvt, _cfg.step2Clustering.minNumNeigh,
                    _cfg.step2Clustering.maxNumNeigh, _cfg.step2Clustering.numEllipsoids, _cfg.step2Clustering.maxEllipsoidSize
            );

            if ( relocatedEvCat )
            {
                SEISCOMP_INFO("Alternative step 1 relocation successful");
                SEISCOMP_INFO("%s", relocationReport(relocatedEvCat).c_str() );

                evToRelocateCat = relocatedEvCat;
            }
            else
            {
                SEISCOMP_ERROR("Failed to perform alternative step 1 origin relocation");
            }
        }
    }

    //
    // Step 2: relocate the refined location this time with cross correlation
    //
    SEISCOMP_INFO("Performing step 2: relocation with cross correlation (with theoretical phases)");

    eventWorkingDir = (boost::filesystem::path(subFolder)/"step2").string();

    CatalogPtr relocatedEvWithXcorr = relocateEventSingleStep(
            evToRelocateCat, eventWorkingDir, true, _cfg.artificialPhases.enable,
            _cfg.hypodd.step2CtrlFile, _cfg.step2Clustering.minWeight,
            _cfg.step2Clustering.minESdist, _cfg.step2Clustering.maxESdist, _cfg.step2Clustering.minEStoIEratio,
            _cfg.step2Clustering.minDTperEvt, _cfg.step2Clustering.maxDTperEvt, _cfg.step2Clustering.minNumNeigh,
            _cfg.step2Clustering.maxNumNeigh, _cfg.step2Clustering.numEllipsoids, _cfg.step2Clustering.maxEllipsoidSize
    );

    if ( relocatedEvWithXcorr )
    {
        SEISCOMP_INFO("Step 2 relocation successful");
        SEISCOMP_INFO("%s", relocationReport(relocatedEvWithXcorr).c_str() );
    }
    else
    {
        SEISCOMP_ERROR("Failed to perform step 2 origin relocation");
    }

    if ( ! relocatedEvWithXcorr )
        throw runtime_error("Failed origin relocation");

    if ( _workingDirCleanup ) boost::filesystem::remove_all(subFolder);

    return relocatedEvWithXcorr;
}



CatalogPtr 
HypoDD::relocateEventSingleStep(const CatalogCPtr& evToRelocateCat,
                                const string& workingDir,
                                bool doXcorr,
                                bool computeTheoreticalPhases,
                                string hypoddCtrlFile,
                                double minPhaseWeight,
                                double minESdist,
                                double maxESdist,
                                double minEStoIEratio,
                                int minDTperEvt,
                                int maxDTperEvt,
                                int minNumNeigh,
                                int maxNumNeigh,
                                int numEllipsoids,
                                double maxEllipsoidSize)
{
    if ( !Util::createPath(workingDir) )
    {
        string msg = "Unable to create working directory: " + workingDir;
        throw runtime_error(msg);
    }

    // extract event to relocate
    Catalog::Event evToRelocate = evToRelocateCat->getEvents().begin()->second;

    //
    // Select neighbouring events
    //
    int numNeighbours;
    bool keepUnmatchedPhases = doXcorr; //useful for detecting missed picks

    CatalogPtr neighbourCat = selectNeighbouringEvents(
        _ddbgc, evToRelocate, evToRelocateCat, minPhaseWeight, minESdist,  maxESdist,
        minEStoIEratio, minDTperEvt,  maxDTperEvt, minNumNeigh, maxNumNeigh,
        numEllipsoids, maxEllipsoidSize, keepUnmatchedPhases, &numNeighbours
    );

    // add event to the neighbours
    unsigned evToRelocateNewId = neighbourCat->add(evToRelocate.id, *evToRelocateCat, false);

    // add numNeighbours information computed by selectNeighbouringEvents
    evToRelocate = neighbourCat->getEvents().at(evToRelocateNewId);
    evToRelocate.relocInfo.numNeighbours = numNeighbours;
    neighbourCat->updateEvent(evToRelocate);

    // write catalog for debugging purpose
    if ( ! _workingDirCleanup )
    {
        neighbourCat->writeToFile(
            (boost::filesystem::path(workingDir)/"starting-event.csv").string(),
            (boost::filesystem::path(workingDir)/"starting-phase.csv").string(),
            (boost::filesystem::path(workingDir)/"starting-station.csv").string());
    }

    // Create station.dat for hypodd
    string stationFile = (boost::filesystem::path(workingDir)/"station.dat").string();
    createStationDatFile(neighbourCat, stationFile);

    // Create event.dat for hypodd
    string eventFile = (boost::filesystem::path(workingDir)/"event.dat").string();
    createEventDatFile(neighbourCat, eventFile); 

    // Create cross correlated differential travel times file (dt.cc) for hypodd
    string dtccFile = (boost::filesystem::path(workingDir)/"dt.cc").string();
    if ( doXcorr )
    {
        // Perform cross correlation, which also detects picks around theoretical
        // arrival times. The catalog will be updated with those theoretical phases 
        const XCorrCache xcorr = buildXCorrCache(neighbourCat, evToRelocateNewId, computeTheoreticalPhases);
        createDtCc(neighbourCat, evToRelocateNewId, dtccFile, xcorr);
    }
    else
    {
        // Create an empty cross correlated differential travel times file (dt.cc) for hypodd
        ofstream(dtccFile).close();
    }

    // Create differential travel times file (dt.ct) for hypodd
    string dtctFile = (boost::filesystem::path(workingDir)/"dt.ct").string();
    createDtCt(neighbourCat, evToRelocateNewId, dtctFile);

    // run hypodd
    // input : dt.cc dt.ct event.sel station.sel hypoDD.inp
    // output : hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
    runHypodd(workingDir, dtccFile, dtctFile, eventFile, stationFile, hypoddCtrlFile);

    // Load the relocated origin from Hypodd
    string ddrelocFile = (boost::filesystem::path(workingDir)/"hypoDD.reloc").string();
    string ddresidualFile = (boost::filesystem::path(workingDir)/"hypoDD.res").string();
    CatalogPtr relocatedCatalog = loadRelocatedCatalog(neighbourCat, ddrelocFile, ddresidualFile);
    CatalogPtr relocatedEvCat = relocatedCatalog->extractEvent(evToRelocateNewId, true);

    // sometimes hypoDD.reloc file is there but it doesn't contain the relocated event 
    Catalog::Event firstAndOnlyEv = relocatedEvCat->getEvents().begin()->second;
    if ( ! firstAndOnlyEv.relocInfo.isRelocated )
        relocatedEvCat.reset();

    // write catalog for debugging purpose
    if ( ! _workingDirCleanup )
    {
        relocatedCatalog->writeToFile(
            (boost::filesystem::path(workingDir)/"relocated-event.csv").string(),
            (boost::filesystem::path(workingDir)/"relocated-phase.csv").string(),
            (boost::filesystem::path(workingDir)/"relocated-station.csv").string());
    }

    return relocatedEvCat;
}



string HypoDD::relocationReport(const CatalogCPtr& relocatedEv)
{
    Catalog::Event event = relocatedEv->getEvents().begin()->second;
    if ( ! event.relocInfo.isRelocated )
        return "Event not relocated";

    return stringify("Neighboring events %d. "
                     "Cross-correlated P phases %d, S phases %d. Rms residual %.3f [sec]. "
                     "Catalog P phases %d, S phases %d. Rms residual %.2f [sec]. "
                     "Error [km]: East-west %.3f, north-south %.3f, depth %.3f",
                      event.relocInfo.numNeighbours,
                      event.relocInfo.numCCp, event.relocInfo.numCCs, event.relocInfo.rmsResidualCC,
                      event.relocInfo.numCTp, event.relocInfo.numCTs, event.relocInfo.rmsResidualCT,
                      event.relocInfo.lonUncertainty, event.relocInfo.latUncertainty,
                      event.relocInfo.depthUncertainty);
 
}


CatalogPtr HypoDD::selectNeighbouringEvents(const CatalogCPtr& catalog,
                                            const Catalog::Event& refEv,
                                            const CatalogCPtr& refEvCatalog,
                                            double minPhaseWeight,
                                            double minESdist,
                                            double maxESdist,
                                            double minEStoIEratio,
                                            int minDTperEvt,
                                            int maxDTperEvt,
                                            int minNumNeigh,
                                            int maxNumNeigh,
                                            int numEllipsoids,
                                            double maxEllipsoidSize,
                                            bool keepUnmatched,
                                            int* numNeigh) const
{
    SEISCOMP_INFO("Selecting Neighbouring Events for event %s lat %.6f lon %.6f depth %.4f mag %.2f time %s",
                   string(refEv).c_str(), refEv.latitude, refEv.longitude, refEv.depth,
                   refEv.magnitude, refEv.time.iso().c_str());

    // Optimization: make code faster but the result will be the same
    if ( maxNumNeigh <= 0 )
    {
        SEISCOMP_INFO("Disabling ellipsoid algorithm since maxNumNeigh is not set");
        numEllipsoids = 0;
    }

    const map<string,Catalog::Station>& stations = catalog->getStations();
    const map<unsigned,Catalog::Event>& events = catalog->getEvents();
    const multimap<unsigned,Catalog::Phase>& phases = catalog->getPhases();

    /*
     * Build ellipsoids
     *
     * From Waldhauser 2009: to assure a spatially homogeneous subsampling, reference
     * events are selected within each of five concentric, vertically elongated
     * ellipsoidal layers of increasing thickness. Each layer has 8 quadrants.
     */
    vector<HddEllipsoidPtr> ellipsoids;
    double verticalSize = maxEllipsoidSize * 2; // horizontal to vertical axis length 
    for (int i = 0; i < (numEllipsoids-1); i++)
    {
        ellipsoids.push_back( new HddEllipsoid(verticalSize, refEv.latitude, refEv.longitude, refEv.depth) );
        verticalSize /= 2;
    }
    ellipsoids.push_back( new HddEllipsoid(0, verticalSize, refEv.latitude, refEv.longitude, refEv.depth) );

    //
    // sort catalog events by distance and drop the ones further than the outmost ellipsoid
    //
    map<unsigned,double> distanceByEvent; // eventid, distance
    map<unsigned,double> azimuthByEvent;  // eventid, azimuth

    for (const auto& kv : events )
    {
        const Catalog::Event& event = kv.second;

        if (event == refEv)
            continue;

        // drop event if outside the outmost ellipsod boundaries
        HddEllipsoidPtr outmostEllip = ellipsoids[0];
        if ( ! outmostEllip->getOuterEllipsoid().isInside(event.latitude, event.longitude, event.depth) )
            continue;

        // compute distance between current event and reference origin
        double azimuth;
        double distance = computeDistance(refEv, event, &azimuth);

        // keep a list of added events
        distanceByEvent[event.id] = distance;
        azimuthByEvent[event.id]  = azimuth;
    }

    //
    // From the events within distance select the ones who respect the constraints
    //
    multimap<double,CatalogPtr> selectedEvents; // distance, event
    map<unsigned,int> dtCountByEvent;           // eventid, dtCount
    set<string> includedStations;
    set<string> excludedStations;

    for (const auto& kv : distanceByEvent)
    {
        const Catalog::Event& event = events.at(kv.first);
        const double eventDistance = kv.second;

        // if the constraints are met evCat will be added to selectedEvents
        CatalogPtr evCat = catalog->extractEvent(event.id, true);

        // keep track of station distance
        multimap<double, pair<string,string> > stationByDistance; // distance, <stationid,phaseType>
        multimap<double, pair<string,string> > unmatchedPhases; // distance, <stationid,phaseType>

        // Check enough phases (> minDTperEvt) ?
        auto eqlrng = phases.equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
            const Catalog::Phase& phase = it->second;

            // check pick weight
            if (phase.procInfo.weight < minPhaseWeight)
            {
                evCat->removePhase(event.id, phase.stationId, phase.procInfo.type);
                continue;
            }

            // check events/station distance
            const Catalog::Station& station = stations.at(phase.stationId);

            if ( excludedStations.find(station.id) != excludedStations.end() )
            {
                evCat->removePhase(event.id, phase.stationId, phase.procInfo.type);
                continue;
            }

            if ( includedStations.find(station.id) == includedStations.end() )
            {
                // compute distance between station and reference event
                double stationDistance = computeDistance(refEv, station);

                // check this station distance is ok
                if ( ( maxESdist > 0 && stationDistance > maxESdist ) ||  // too far away ?
                     ( stationDistance < minESdist ) )                    // too close ?
                {
                    excludedStations.insert(station.id);
                    evCat->removePhase(event.id, phase.stationId, phase.procInfo.type);
                    continue;
                }

                if ( (stationDistance/eventDistance) < minEStoIEratio ) // ratio too small ?
                {
                    // since this is dependents on the current event we cannot save it into excludedStations 
                    evCat->removePhase(event.id, phase.stationId, phase.procInfo.type);
                    continue;
                }

                // this station is ok for refEv
                includedStations.insert(station.id);
            }

            // compute distance between current event and station
            double stationDistance = computeDistance(event, station);

            // check this station distance is ok
            if ( ( maxESdist > 0 && stationDistance > maxESdist ) ||      // too far away ?
                 ( stationDistance < minESdist )                 ||       // too close ?
                 ( (stationDistance / eventDistance) < minEStoIEratio ) ) // ratio too small ?
            {
                evCat->removePhase(event.id, phase.stationId, phase.procInfo.type);
                continue;
            }

            // now find corresponding phase in reference event phases
            bool peer_found = false;
            auto itRef = refEvCatalog->searchPhase(refEv.id, phase.stationId, phase.procInfo.type);
            if ( itRef != refEvCatalog->getPhases().end() )
            {
                const Catalog::Phase& refPhase = itRef->second;
                if (refPhase.procInfo.weight >= minPhaseWeight)
                    peer_found = true;
            }

            if ( ! peer_found )
            {
                // don't delete this phase, because it respects the constraints and if we add
                // theoretical picks to refEv those phases might become useful for xcorr
                if ( keepUnmatched )
                {
                    unmatchedPhases.emplace(stationDistance, pair<string,string>(phase.stationId, phase.procInfo.type));
                }
                else
                {
                    evCat->removePhase(event.id, phase.stationId, phase.procInfo.type);
                }
                continue;
            }

            stationByDistance.emplace(stationDistance, pair<string,string>(phase.stationId, phase.procInfo.type));
        }

        int dtCount = stationByDistance.size();

        // if not enought phases skip event
        if ( dtCount < minDTperEvt )
        {
            continue;
        }

        // if maxDTperEvt is set then make sure to stay within limits
        if ( maxDTperEvt > 0 )
        {
            if ( dtCount > maxDTperEvt )
            {
                // remove phases belonging to further stations from matched phases
                auto first = std::next(stationByDistance.begin(), maxDTperEvt);
                std::for_each(first, stationByDistance.end(),
                        [evCat, event](const pair<double,pair<string,string>>& kv) { // kv == <distance, <stationid,phaseType> >
                            evCat->removePhase(event.id, kv.second.first, kv.second.second); }
                );
                dtCount = maxDTperEvt;
            }

            if ( (dtCount + unmatchedPhases.size()) > maxDTperEvt )
            {
                // remove phases belonging to further stations from unmatched phases
                auto first = std::next(unmatchedPhases.begin(), maxDTperEvt - stationByDistance.size());
                std::for_each(first, unmatchedPhases.end(),
                        [evCat, event](const pair<double,pair<string,string>>& kv) { // kv == <distance, <stationid,phaseType> >
                            evCat->removePhase(event.id, kv.second.first, kv.second.second); }
                );
                dtCount = maxDTperEvt;
            }
        }

        // add this event to the selected ones
        selectedEvents.emplace(eventDistance, evCat);
        dtCountByEvent.emplace(event.id, dtCount);
    }

    //
    // Finally build the catalog of neighboring events using the elipsoids
    // algorithm or simply the nearest neighbour method
    //
    CatalogPtr neighboringEventCat = new Catalog();
    int numNeighbors = 0;

    if ( numEllipsoids <= 0 )
    {
        //
        // if numEllipsoids is 0 then disable the ellipsoid algorithm and simply select events
        // on the nearest neighbor basis
        // Since selectedEvents is sorted by distance we get closer events first
        //
        for (auto kv : selectedEvents)
        {
            CatalogPtr evCat = kv.second;
            const Catalog::Event& ev = evCat->getEvents().begin()->second;

            // add this event to the catalog
            neighboringEventCat->add(ev.id, *evCat, true);
            numNeighbors++;
            SEISCOMP_INFO("Chose neighbour dtCount %2d distance %5.2f azimuth %3.f depth-diff %6.3f depth %5.3f event %s ",
                            dtCountByEvent[ev.id], distanceByEvent[ev.id], azimuthByEvent[ev.id], refEv.depth-ev.depth, ev.depth, string(ev).c_str() );
            if ( maxNumNeigh > 0 && numNeighbors >= maxNumNeigh) break;
        }
    }
    else
    {
        //
        // Apply Waldauser's concentric ellipsoids algorithm
        //
        vector<int> quadrants = {1,2,3,4,5,6,7,8};
        bool workToDo = true;

        while ( workToDo )
        {
            for(int elpsNum = ellipsoids.size() -1; elpsNum >= 0;  elpsNum--)
            {
                for (int quadrant : quadrants)
                {
                    // if we don't have events or we have already selected maxNumNeigh neighbors exit
                    if ( selectedEvents.empty() || (maxNumNeigh > 0 && numNeighbors >= maxNumNeigh) )
                    {
                        workToDo = false;
                        break;
                    }

                    // since selectedEvents is sorted by distance we get closer events first
                    for (auto it = selectedEvents.begin(); it != selectedEvents.end(); it++)
                    {
                        CatalogPtr evCat = it->second;
                        const Catalog::Event& ev = evCat->getEvents().begin()->second;

                        bool found = ellipsoids[elpsNum]->isInside(ev.latitude, ev.longitude, ev.depth, quadrant);

                        if ( found )
                        {
                            // add this event to the catalog
                            neighboringEventCat->add(ev.id, *evCat, true);
                            numNeighbors++;
                            selectedEvents.erase(it);
                            SEISCOMP_INFO("Chose neighbour ellipsoid %2d quadrant %d dtCount %2d distance %5.2f azimuth %3.f depth-diff %6.3f depth %5.3f event %s ",
                                           elpsNum, quadrant, dtCountByEvent[ev.id], distanceByEvent[ev.id], azimuthByEvent[ev.id], refEv.depth-ev.depth,
                                           ev.depth, string(ev).c_str() ); 
                            break;
                        }
                    }
                }
            }
        }
    }

    if ( numNeigh ) *numNeigh = numNeighbors;

    // Check if enough neighbors were found
    if ( numNeighbors < minNumNeigh )
    {
        string msg = stringify("Skipping event %s, insufficient number of neighbors (%d)",
                               string(refEv).c_str(), numNeighbors);
        SEISCOMP_DEBUG("%s", msg.c_str());
        throw runtime_error(msg);
    }

    return neighboringEventCat;
}



map<unsigned,CatalogPtr> 
HypoDD::selectNeighbouringEventsCatalog(const CatalogCPtr& catalog,
                                        double minPhaseWeight,
                                        double minESdist,
                                        double maxESdist,
                                        double minEStoIEratio,
                                        int minDTperEvt,
                                        int maxDTperEvt,
                                        int minNumNeigh,
                                        int maxNumNeigh,
                                        int numEllipsoids,
                                        double maxEllipsoidSize,
                                        bool keepUnmatched) const
{
    SEISCOMP_INFO("Selecting Catalog Neighbouring Events ");

    // neighbours for each event
    map<unsigned,CatalogPtr> neighboursByEvent;

    // for each event find the neighbours
    for (const auto& kv : catalog->getEvents() )
    {
        Catalog::Event event = kv.second;
        int numNeighbours;

        CatalogPtr neighbourCat; 
        try {
            neighbourCat = selectNeighbouringEvents(
                catalog, event, catalog, minPhaseWeight, minESdist, maxESdist,
                minEStoIEratio, minDTperEvt, maxDTperEvt,  minNumNeigh, maxNumNeigh,
                numEllipsoids, maxEllipsoidSize, keepUnmatched, &numNeighbours
            );
        } catch ( ... ) { }

        if ( neighbourCat )
        {
            // add event to neighbour catalog list
            neighbourCat->add(event.id,*catalog, true);
            // add numNeighbours information
            event.relocInfo.numNeighbours = numNeighbours;
            neighbourCat->updateEvent(event);
            // update the list of events with their respective neighbours
            neighboursByEvent[event.id] = neighbourCat;
        }
    }

    // We don't want to report the same pairs multiple times
    // when creating dt.cc and dt.ct (e.g. pair eventXX-eventYY
    // is the same as pair eventYY-eventXX), we'll remove the
    // pairs that appeared in previous catalogs from the
    // following catalogs
    std::multimap<unsigned,unsigned> existingPairs;

    for (auto kv : neighboursByEvent )
    {
        unsigned currEventId = kv.first;
        CatalogPtr& neighbourCat  = kv.second;

        // remove from currrent catalogl the existing pairs
        auto eqlrng = existingPairs.equal_range(currEventId);
        for (auto existingPair = eqlrng.first; existingPair != eqlrng.second; existingPair++)
        {
            neighbourCat->removeEvent(existingPair->second);
        }

        // remove current pairs from following catalogs
        for ( auto const& kv : neighbourCat->getEvents() )
        {
            if ( kv.first != currEventId )
                existingPairs.emplace(kv.first, currEventId);
        }
    }

    return neighboursByEvent;
}



void
HypoDD::addMissingEventPhases(const CatalogCPtr& searchCatalog,
                              const Catalog::Event& refEv,
                              CatalogPtr& refEvCatalog)
{
    std::vector<Catalog::Phase> newPhases = findMissingEventPhases(searchCatalog, refEv, refEvCatalog);

    for (Catalog::Phase& ph : newPhases)
    {
        refEvCatalog->removePhase(ph.eventId, ph.stationId, ph.procInfo.type);
        refEvCatalog->addPhase(ph, false, false);
        const Catalog::Station& station = searchCatalog->getStations().at(ph.stationId);
        refEvCatalog->addStation(station, true);
    }
}



std::vector<Catalog::Phase>
HypoDD::findMissingEventPhases(const CatalogCPtr& searchCatalog,
                               const Catalog::Event& refEv,
                               const CatalogPtr& refEvCatalog)
{
    //
    // find stations for which the refEv doesn't have phases
    //
    vector<MissingStationPhase> missingPhases = getMissingPhases(searchCatalog, refEv, refEvCatalog);

    //
    // for each missed phase try to detect it
    //
    std::vector<Catalog::Phase> newPhases;
    for ( const MissingStationPhase& pair : missingPhases )
    {
        const Catalog::Station& station = searchCatalog->getStations().at(pair.first);
        const string phaseType = pair.second;

        //
        // loop through each other event and select the ones who have a manually picked phase for the missing station
        //
        vector<HypoDD::PhasePeer> peers = findPhasePeers(station, phaseType, searchCatalog);

        if ( peers.size() <= 0 )
        {
            continue;
        }

        // compute velocity using existing background catalog phases
        double phaseVelocity = 0;
        for (const PhasePeer& peer  : peers)
        {
            const Catalog::Event& event = peer.first;
            const Catalog::Phase& phase = peer.second;
            double travelTime = (phase.time - event.time).length();
            double stationDistance = computeDistance(event, station);
            double vel = stationDistance / travelTime;
            phaseVelocity += vel;
        }
        phaseVelocity /= peers.size();

        Catalog::Phase refEvNewPhase = createThoreticalPhase(station, phaseType, refEv, peers, phaseVelocity);

        newPhases.push_back(refEvNewPhase);
    }

    return newPhases;
}



vector<HypoDD::MissingStationPhase>
HypoDD::getMissingPhases(const CatalogCPtr& searchCatalog,
                         const Catalog::Event& refEv,
                         const CatalogPtr& refEvCatalog) const
{
    const auto& refEvPhases = refEvCatalog->getPhases().equal_range(refEv.id);

    //
    // loop through stations and find those for which the refEv doesn't have phases
    //
    vector<MissingStationPhase> missingPhases;
    for (const auto& kv : searchCatalog->getStations() )
    {
        const Catalog::Station& station = kv.second;

        bool foundP = false, foundS = false;
        for (auto it = refEvPhases.first; it != refEvPhases.second; ++it)
        {
            const Catalog::Phase& phase = it->second;
            if ( station.networkCode == phase.networkCode &&
                 station.stationCode == phase.stationCode)
            {
                if ( phase.procInfo.type == "P" ) foundP = true;
                if ( phase.procInfo.type == "S" ) foundS = true;
            }
            if ( foundP and foundS ) break;
        }
        if ( ! foundP || ! foundS )
        {
            if ( ! foundP )
                missingPhases.push_back( MissingStationPhase(station.id,"P") );
            if ( ! foundS )
                missingPhases.push_back( MissingStationPhase(station.id,"S") );
        }
    }

    return missingPhases;
}



vector<HypoDD::PhasePeer>
HypoDD::findPhasePeers(const Catalog::Station& station, const std::string& phaseType,
                       const CatalogCPtr& searchCatalog) const
{
    //
    // loop through each other event and select the manual phases for the station we
    // are interested in
    //
    vector<PhasePeer> phasePeers;

    for (const auto& kv : searchCatalog->getEvents() )
    {
        const Catalog::Event& event = kv.second; 
        const auto& phases = searchCatalog->getPhases().equal_range(event.id);

        for (auto it = phases.first; it != phases.second; ++it)
        {
            const Catalog::Phase& phase = it->second;

            if ( station.networkCode == phase.networkCode &&
                 station.stationCode == phase.stationCode &&
                 phaseType           == phase.procInfo.type )
            {
                if ( phase.isManual )
                {
                    phasePeers.push_back( PhasePeer(event, phase) );
                }
                break;
            }
        }
    }

    return phasePeers;
}



Catalog::Phase
HypoDD::createThoreticalPhase(const Catalog::Station& station,
                              const string& phaseType,
                              const Catalog::Event& refEv,
                              const vector<HypoDD::PhasePeer>& peers,
                              double phaseVelocity)
{
    const auto xcorrCfg = _cfg.xcorr.at(phaseType);

    // detect locationCode/channelCode   
    struct {
        string locationCode;
        string channelCode;
        Core::Time time;
    } streamInfo = {"", "", Core::Time()};

    for (const PhasePeer& peer : peers)
    {
        //const Catalog::Event& event = peer.first;
        const Catalog::Phase& phase = peer.second;
        // get the closest in time to refEv stream information
        if ( (refEv.time - phase.time).abs() < (refEv.time - streamInfo.time).abs() )
            streamInfo = {phase.locationCode, phase.channelCode, phase.time};
    }

    // initialize the new phase
    Catalog::Phase refEvNewPhase;

    refEvNewPhase.eventId      = refEv.id;
    refEvNewPhase.stationId    = station.id;
    refEvNewPhase.networkCode  = station.networkCode;
    refEvNewPhase.stationCode  = station.stationCode;
    refEvNewPhase.locationCode = streamInfo.locationCode;
    refEvNewPhase.channelCode  = WfMngr::getBandAndInstrumentCodes(streamInfo.channelCode) + xcorrCfg.components[0];
    refEvNewPhase.isManual     = false;
    refEvNewPhase.procInfo.type = phaseType;

    // use phase velocity to compute phase time
    double stationDistance = computeDistance(refEv, station);
    refEvNewPhase.time = refEv.time + Core::TimeSpan(stationDistance / phaseVelocity);

    refEvNewPhase.lowerUncertainty = Catalog::DEFAULT_AUTOMATIC_PICK_UNCERTAINTY;
    refEvNewPhase.upperUncertainty = Catalog::DEFAULT_AUTOMATIC_PICK_UNCERTAINTY;
    refEvNewPhase.procInfo.weight = computePickWeight(refEvNewPhase);
    refEvNewPhase.procInfo.source = Catalog::Phase::Source::THEORETICAL;
    refEvNewPhase.type = phaseType + "t";

    return refEvNewPhase;
}


HypoDD::XCorrCache
HypoDD::buildXCorrCache(std::map<unsigned,CatalogPtr>& neighbourCats,
                       bool computeTheoreticalPhases)
{
    XCorrCache xcorr;
    _counters = {0};
    _wf->resetCounters();

    for (auto& kv : neighbourCats)
    {
        unsigned refEvId = kv.first;
        CatalogPtr& neighbourCat = kv.second;

        const Catalog::Event& refEv = neighbourCat->getEvents().at(refEvId);

        // Compute theoretical phases for stations that have no picks. The cross correlation will
        // be used to detect and fix pick time
        if ( computeTheoreticalPhases )
        {
            addMissingEventPhases(neighbourCat, refEv, neighbourCat);
        }

        buildXcorrDiffTTimePairs(neighbourCat, refEv, xcorr);

        // Update theoretical and automatic phase pick time and uncertainties based on
        // cross-correlation results
        // drop theoretical phases without good cross correlation results
        if ( computeTheoreticalPhases )
        {
            fixPhases(neighbourCat, refEv, xcorr);
        } 
    }

    printCounters();

    return xcorr;
}


HypoDD::XCorrCache
HypoDD::buildXCorrCache(CatalogPtr& catalog, unsigned evToRelocateId,
                        bool computeTheoreticalPhases)
{
    map<unsigned,CatalogPtr> neighbourCats = {{evToRelocateId,catalog}};
    return buildXCorrCache(neighbourCats, computeTheoreticalPhases);
}


/*
 * Compute and store to XCorrCache cross-correlated differential travel times
 * for pairs of earthquake
 */
void HypoDD::buildXcorrDiffTTimePairs(CatalogPtr& catalog,
                                      const Catalog::Event& refEv,
                                      XCorrCache& xcorr)
{
    SEISCOMP_INFO("Computing cross-correlation differential travel times for event %s",
                  string(refEv).c_str() );

    // xcorr settings depending on the phase type
    // (NO snr check for theoretical phases because the pick time is likely wrong and fixed later)
    map<Catalog::Phase::Source, PhaseXCorrCfg> phCfg = {
        {Catalog::Phase::Source::CATALOG,      {_useCatalogDiskCache, &_wfCache, true }},
        {Catalog::Phase::Source::RT_EVENT,     {false,                nullptr,   true }},
        {Catalog::Phase::Source::THEORETICAL,  {false,                nullptr,   false}}
    };

    // keep track of refEv distant to stations
    multimap<double,string> stationByDistance; // <distance, stationid>
    set<string> computedStations;

    //
    // loop through reference event phases
    //
    auto eqlrngRef = catalog->getPhases().equal_range(refEv.id);
    for (auto itRef = eqlrngRef.first; itRef != eqlrngRef.second; ++itRef)
    {
        const Phase& refPhase = itRef->second;

        // save events/station distance
        if ( computedStations.find(refPhase.stationId) == computedStations.end() )
        {
            const Station& station = catalog->getStations().at(refPhase.stationId);
            double stationDistance = computeDistance(refEv, station);
            stationByDistance.emplace(stationDistance, refPhase.stationId);
            computedStations.insert(refPhase.stationId);
        }

        //
        // loop through catalog events and cross correlate phase pairs
        //
        for (const auto& kv : catalog->getEvents() )
        {
            const Event& event = kv.second;

            if (event == refEv)
                continue;

            auto it = catalog->searchPhase(event.id, refPhase.stationId, refPhase.procInfo.type);

            if ( it != catalog->getPhases().end() )
            {
                const Phase& phase = it->second;

                double coeff, lag, dtcc, weight;
                if ( xcorrPhases(refEv, refPhase, phCfg.at(refPhase.procInfo.source),
                                 event, phase, phCfg.at(phase.procInfo.source),
                                 coeff, lag, dtcc, weight) )
                {
                    //
                    // Store good xcorr results
                    //
                    auto& entry = xcorr.getForUpdate(refEv.id, refPhase.stationId, refPhase.procInfo.type);
                    entry.update(event, phase, coeff, lag, dtcc, weight);
                }
            }
        }
    }

    // finalize statistics
    xcorr.computeStats();

    // Print some useful information
    for (const auto& kv : stationByDistance)
    {
        const double stationDistance = kv.first;
        const Station& station = catalog->getStations().at(kv.second);

        bool goodPXcorr = xcorr.has(refEv.id, station.id, Phase::Type::P);
        bool goodSXcorr = xcorr.has(refEv.id, station.id, Phase::Type::S);

        if ( ! goodPXcorr && ! goodSXcorr )
        {
            SEISCOMP_INFO("xcorr: event %5s sta %4s %5s dist %.2f [km] - no good cross correlations pairs",
                          string(refEv).c_str(), station.networkCode.c_str(),
                          station.stationCode.c_str(), stationDistance);
        }
        else
        {
            if ( goodPXcorr )
            {
                const auto& pdata = xcorr.get(refEv.id, station.id, Phase::Type::P);
                SEISCOMP_INFO("xcorr: event %5s sta %4s %5s dist %.2f [km] - "
                          "%d P phases, mean coeff %.2f lag %.2f (events: %s)",
                          string(refEv).c_str(), station.networkCode.c_str(),
                          station.stationCode.c_str(), stationDistance,
                          pdata.ccCount, pdata.mean_coeff, pdata.mean_lag,
                          pdata.peersStr.c_str());
            }
            if ( goodSXcorr )
            {
                const auto& sdata = xcorr.get(refEv.id, station.id, Phase::Type::S);
                SEISCOMP_INFO("xcorr: event %5s sta %4s %5s dist %.2f [km] - "
                          "%d S phases, mean coeff %.2f lag %.2f (events: %s)",
                          string(refEv).c_str(), station.networkCode.c_str(),
                          station.stationCode.c_str(), stationDistance,
                          sdata.ccCount, sdata.mean_coeff, sdata.mean_lag,
                          sdata.peersStr.c_str() );
            }
        }
    }
}

/*
 * Update theoretical and automatic phase pick time and uncertainties based on
 * cross-correlation results
 * Also remove theoretical phases without good cross correlation results
 */
void HypoDD::fixPhases(CatalogPtr& catalog, const Catalog::Event& refEv, XCorrCache& xcorr)
{
    unsigned newP = 0, newS = 0;

    std::vector<Catalog::Phase> phasesToBeRemoved;
    std::vector<Catalog::Phase> newPhases;

    auto eqlrng = catalog->getPhases().equal_range(refEv.id);
    for (auto it = eqlrng.first; it != eqlrng.second; it++)
    {
        const Catalog::Phase& phase = it->second;
        bool goodXcorr = xcorr.has(refEv.id, phase.stationId, phase.procInfo.type);

        // nothing to do if we dont't have good xcorr results of if the phase is manual
        if ( ! goodXcorr || phase.isManual )
        {
            // remove thoretical phases wihtout good xcorr results
            if ( phase.procInfo.source == Catalog::Phase::Source::THEORETICAL )
                phasesToBeRemoved.push_back(phase);
            continue;
        }

        const auto& pdata = xcorr.get(refEv.id, phase.stationId, phase.procInfo.type);

        //
        // Set new phase time and uncertainty
        //
        Catalog::Phase newPhase(phase);
        newPhase.time  -= Core::TimeSpan(pdata.mean_lag);
        newPhase.lowerUncertainty = pdata.mean_lag - pdata.min_lag;
        newPhase.upperUncertainty = pdata.max_lag - pdata.mean_lag;
        newPhase.procInfo.weight = computePickWeight(newPhase);
        newPhase.procInfo.source = Catalog::Phase::Source::XCORR;
        newPhase.type = newPhase.procInfo.type + "x";

        if ( phase.procInfo.source == Catalog::Phase::Source::THEORETICAL )
        {
            newP += newPhase.procInfo.type == "P" ? 1 : 0;
            newS += newPhase.procInfo.type == "S" ? 1 : 0;
        }

        // remove the old phase since the new one will be added
        phasesToBeRemoved.push_back(phase);
        newPhases.push_back(newPhase);
    }

    //
    // Replace automatic/theoretical phases with xcorr detected ones
    //
    for (const Catalog::Phase& ph : phasesToBeRemoved)
    {
        catalog->removePhase(ph);
    }

    for (Catalog::Phase& ph : newPhases)
    {
        catalog->addPhase(ph, false, false);
    }

    const auto& refEvPhases = catalog->getPhases().equal_range(refEv.id);
    SEISCOMP_INFO("Event %s total phases %lu, created %u new phases (%u P and %u S)",
                   string(refEv).c_str(), std::distance(refEvPhases.first, refEvPhases.second),
                   (newP+newS), newP, newS );
}


void HypoDD::printCounters()
{
    unsigned performed        = _counters.xcorr_performed,
             performed_s      = _counters.xcorr_performed_s,
             performed_p      = performed - performed_s,
             perf_theo        = _counters.xcorr_performed_theo,
             perf_theo_s      = _counters.xcorr_performed_s_theo,
             perf_theo_p      = perf_theo - perf_theo_s,
             good_cc          = _counters.xcorr_good_cc,
             good_cc_s        = _counters.xcorr_good_cc_s,
             good_cc_p        = good_cc - good_cc_s,
             good_cc_theo     = _counters.xcorr_good_cc_theo,
             good_cc_s_theo   = _counters.xcorr_good_cc_s_theo,
             good_cc_p_theo   = good_cc_theo - good_cc_s_theo;

    unsigned snr_low, wf_no_avail, wf_cached, wf_downloaded;
    _wf->getCounters(snr_low, wf_no_avail, wf_cached, wf_downloaded);

    SEISCOMP_INFO("Cross correlation performed %u, "
                  "phases with Signal to Noise ratio too low %u, "
                  "phases not available %u (waveforms downloaded %u, "
                  "waveforms loaded from disk cache %u)",
                  performed, snr_low, wf_no_avail, wf_downloaded, wf_cached);

    SEISCOMP_INFO("Total xcorr %u (P %.f%%, S %.f%%) success %.f%% (%u/%u). Successful P %.f%% (%u/%u). Successful S %.f%% (%u/%u)",
                  performed, (performed_p*100./performed), (performed_s*100./performed),
                  (good_cc*100./performed), good_cc, performed,
                  (good_cc_p*100./performed_p), good_cc_p, performed_p,
                  (good_cc_s*100./performed_s), good_cc_s, performed_s);

    if ( perf_theo > 0 )
    {
         unsigned perf_real        = performed - perf_theo,
                  perf_real_s      = performed_s - perf_theo_s,
                  perf_real_p      = perf_real - perf_real_s,
                  good_cc_real     = good_cc - good_cc_theo,
                  good_cc_s_real   = good_cc_s - good_cc_s_theo,
                  good_cc_p_real   = good_cc_real - good_cc_s_real;

         SEISCOMP_INFO("xcorr on actual picks %u/%u (P %.f%%, S %.f%%) success %.f%% (%u/%u). Successful P %.f%% (%u/%u). Successful S %.f%% (%u/%u)", 
                      perf_real, performed, (perf_real_p*100./perf_real), (perf_real_s*100./perf_real),
                      (good_cc_real*100./perf_real), good_cc_real, perf_real,
                      (good_cc_p_real*100./perf_real_p), good_cc_p_real, perf_real_p,
                      (good_cc_s_real*100./perf_real_s), good_cc_s_real, perf_real_s);

         SEISCOMP_INFO("xcorr on theoretical picks %u/%u (P %.f%%, S %.f%%) success %.f%% (%u/%u). Successful P %.f%% (%u/%u). Successful S %.f%% (%u/%u)", 
                      perf_theo, performed, (perf_theo_p*100./perf_theo), (perf_theo_s*100./perf_theo),
                      (good_cc_theo*100./perf_theo), good_cc_theo, perf_theo,
                      (good_cc_p_theo*100./perf_theo_p), good_cc_p_theo, perf_theo_p,
                      (good_cc_s_theo*100./perf_theo_s), good_cc_s_theo, perf_theo_s);
    }
}


Core::TimeWindow
HypoDD::xcorrTimeWindowLong(const Catalog::Phase& phase) const
{
    const auto xcorrCfg = _cfg.xcorr.at(phase.procInfo.type);
    Core::TimeWindow tw = xcorrTimeWindowShort(phase);
    tw.setStartTime( tw.startTime() - Core::TimeSpan(xcorrCfg.maxDelay) );
    tw.setEndTime(   tw.endTime()   + Core::TimeSpan(xcorrCfg.maxDelay) );
    return tw;
}


Core::TimeWindow
HypoDD::xcorrTimeWindowShort(const Catalog::Phase& phase) const
{
    const auto xcorrCfg = _cfg.xcorr.at(phase.procInfo.type);
    double shortDuration = xcorrCfg.endOffset - xcorrCfg.startOffset;
    Core::TimeSpan shortTimeCorrection = Core::TimeSpan(xcorrCfg.startOffset);
    return Core::TimeWindow(phase.time + shortTimeCorrection, shortDuration);
}


bool
HypoDD::xcorrPhases(const Catalog::Event& event1, const Catalog::Phase& phase1, PhaseXCorrCfg& phCfg1,
                    const Catalog::Event& event2, const Catalog::Phase& phase2, PhaseXCorrCfg& phCfg2,
                    double& coeffOut, double& lagOut, double& diffTimeOut, double& weightOut)
{
    if ( phase1.procInfo.type != phase2.procInfo.type )
    {
        SEISCOMP_ERROR("Internal logic error: trying to xcorr different phases (%s and %s)",
                       string(phase1).c_str(), string(phase2).c_str());
        return false;
    }

    auto xcorrCfg = _cfg.xcorr.at(phase1.procInfo.type);
    bool performed = false;
    bool goodCoeff = false;

    //
    // Make sure we are using the same channels in cross correlation
    //
    const string channelCodeRoot1 = WfMngr::getBandAndInstrumentCodes(phase1.channelCode);
    const string channelCodeRoot2 = WfMngr::getBandAndInstrumentCodes(phase2.channelCode);

    string commonChRoot;

    if ( channelCodeRoot1 == channelCodeRoot2 )
    {
        commonChRoot = channelCodeRoot1;
    }
    else 
    {
        //
        // if the channel codes don't match then look for a possible match
        //
        DataModel::ThreeComponents dummy;
        DataModel::SensorLocation *loc = WfMngr::findSensorLocation(phase1.networkCode, phase1.stationCode, phase1.locationCode, phase1.time);
        if ( loc && getThreeComponents(dummy, loc, channelCodeRoot2.c_str(), phase1.time) )
        {
            // phase 1 has the same channels of phase 2
            commonChRoot = channelCodeRoot2;
        }
        else
        {
            loc = WfMngr::findSensorLocation(phase2.networkCode, phase2.stationCode, phase2.locationCode, phase2.time);
            if ( loc && getThreeComponents(dummy, loc, channelCodeRoot1.c_str(), phase2.time) )
            {
                // phase 2 has the same channels of phase 1
                commonChRoot = channelCodeRoot1;
            }
        }
    }

    if ( commonChRoot.empty() )
    {
         SEISCOMP_DEBUG("Cannot find common channels to cross correlate %s and %s (%s and %s)",
                        channelCodeRoot1.c_str(), channelCodeRoot2.c_str(), 
                        string(phase1).c_str(), string(phase2).c_str());
    }

    //
    // perform xcorr on all registered component until we get a good correlation coefficient
    //
    for (const string& component : xcorrCfg.components )
    {
        Catalog::Phase tmpPh1 = phase1;
        Catalog::Phase tmpPh2 = phase2;

        // overwrite phases' component for the xcorr
        if ( commonChRoot.empty() )
        {
            tmpPh1.channelCode = channelCodeRoot1 + component;
            tmpPh2.channelCode = channelCodeRoot2 + component;
        }
        else
        {
            tmpPh1.channelCode = commonChRoot + component;
            tmpPh2.channelCode = commonChRoot + component;
        }

        performed = _xcorrPhases(event1, tmpPh1, phCfg1, event2, tmpPh2, phCfg2,
                                 coeffOut, lagOut, diffTimeOut, weightOut);

        goodCoeff = ( performed && std::abs(coeffOut) >= xcorrCfg.minCoef );

        // if the xcorr was successfull and the coeffiecnt is good then stopi here
        if ( goodCoeff )
        {
            break;
        }
    }

    //
    // Deal with counters
    //
    bool isS = ( phase1.procInfo.type == "S" );
    bool isTheoretical = ( phase1.procInfo.source == Catalog::Phase::Source::XCORR       ||
                           phase2.procInfo.source == Catalog::Phase::Source::XCORR       ||
                           phase1.procInfo.source == Catalog::Phase::Source::THEORETICAL ||
                           phase2.procInfo.source == Catalog::Phase::Source::THEORETICAL );

    if ( performed ) 
    {
        _counters.xcorr_performed++;
        if ( isTheoretical ) _counters.xcorr_performed_theo++;
        if ( isS )
        {
            _counters.xcorr_performed_s++;
            if ( isTheoretical ) _counters.xcorr_performed_s_theo++;
        }

        if ( goodCoeff )
        {
            _counters.xcorr_good_cc++;
            if ( isTheoretical ) _counters.xcorr_good_cc_theo++;
            if ( isS )
            {
                _counters.xcorr_good_cc_s++;
                if ( isTheoretical ) _counters.xcorr_good_cc_s_theo++;
            }
        }
    }

    return goodCoeff;
}


bool
HypoDD::_xcorrPhases(const Catalog::Event& event1, const Catalog::Phase& phase1, PhaseXCorrCfg& phCfg1,
                     const Catalog::Event& event2, const Catalog::Phase& phase2, PhaseXCorrCfg& phCfg2,
                     double& coeffOut, double& lagOut, double& diffTimeOut, double& weightOut)
{
    coeffOut = lagOut = diffTimeOut = weightOut = 0;

    auto xcorrCfg = _cfg.xcorr.at(phase1.procInfo.type);

    Core::TimeWindow tw1 = xcorrTimeWindowLong(phase1);
    Core::TimeWindow tw2 = xcorrTimeWindowLong(phase2);

    // load the long trace 1, because we want to cache the long version. Then we'll trim it.
    GenericRecordCPtr tr1 = _wf->getWaveform(tw1, event1, phase1, phCfg1.cache, phCfg1.useDiskCache, phCfg1.allowSnrCheck);
    if ( !tr1 )
    {
        return false;
    }

    // load the long trace 2, because we want to cache the long version. Then we'll trim it
    GenericRecordCPtr tr2 = _wf->getWaveform(tw2, event2, phase2, phCfg2.cache, phCfg2.useDiskCache, phCfg2.allowSnrCheck);
    if ( !tr2 )
    {
        return false;
    }

    // trust the manual pick on phase 2: keet trace2 short and xcorr it with
    // a larger trace1 window
    double xcorr_coeff = 0, xcorr_lag = 0;
 
    if ( phase2.isManual || (! phase1.isManual && ! phase2.isManual) )
    {
        // trim tr2 to shorter length, we want to cross correlate the short with the long one
        GenericRecordPtr tr2Short = new GenericRecord(*tr2);
        Core::TimeWindow tw2Short = xcorrTimeWindowShort(phase2);
        if ( ! WfMngr::trim(*tr2Short, tw2Short) )
        {
            SEISCOMP_DEBUG("Cannot trim phase2 waveform, skipping cross correlation "
                             "for phase pair phase1='%s', phase2='%s'",
                             string(phase1).c_str(), string(phase2).c_str());
            return false;
        }

        if ( ! xcorr(tr1, tr2Short, xcorrCfg.maxDelay, true, xcorr_lag, xcorr_coeff) )
        {
            return false;
        }
    }

    // trust the manual pick on phase 1: keet trace1 short and xcorr it with
    // a larger trace2 window
    double xcorr_coeff2 = 0, xcorr_lag2 = 0;

    if ( phase1.isManual || (! phase1.isManual && ! phase2.isManual) )
    {
        // trim tr1 to shorter length, we want to cross correlate the short with the long one
        GenericRecordPtr tr1Short = new GenericRecord(*tr1);
        Core::TimeWindow tw1Short = xcorrTimeWindowShort(phase1);
        if ( ! WfMngr::trim(*tr1Short, tw1Short) )
        {
            SEISCOMP_DEBUG("Cannot trim phase1 waveform, skipping cross correlation "
                             "for phase pair phase1='%s', phase2='%s'",
                             string(phase1).c_str(), string(phase2).c_str());
            return false;
        }

        if ( ! xcorr(tr1Short, tr2, xcorrCfg.maxDelay, true, xcorr_lag2, xcorr_coeff2) )
        {
            return false;
        }
    }

    if ( std::abs(xcorr_coeff2) > std::abs(xcorr_coeff) )
    {
        // swap
        xcorr_coeff = xcorr_coeff2;
        xcorr_lag = xcorr_lag2;
    }

    // compute differential travel time and weight of measurement
    double travel_time1 = phase1.time - event1.time;
    double travel_time2 = phase2.time - event2.time;
    coeffOut  = xcorr_coeff;
    lagOut    = xcorr_lag;
    diffTimeOut = travel_time1 - travel_time2 - xcorr_lag;
    weightOut = xcorr_coeff * xcorr_coeff;

    return true;
}


/*
 * Calculate the correlation series (tr1 and tr2 are already demeaned)
 *
 * delayOut is the shift in seconds (positive or negative) between tr1 and tr2 middle points
 * to get the highest correlation coefficient (coeffOut) between the 2 traces
 * A delayOut of 0 is when tr1 and tr2 middle points are aligned
 */
bool
HypoDD::xcorr(const GenericRecordCPtr& tr1, const GenericRecordCPtr& tr2, double maxDelay,
              bool qualityCheck, double& delayOut, double& coeffOut) const
{
    coeffOut = std::nan("");

    if (tr1->samplingFrequency() != tr2->samplingFrequency())
    {
        SEISCOMP_INFO("Cannot cross correlate traces with different sampling freq (%f!=%f)",
                      tr1->samplingFrequency(), tr2->samplingFrequency());
        return false;
    }

    const double freq = tr1->samplingFrequency();
    const int maxDelaySmps = maxDelay * freq; // secs to samples

    // check longest/shortest trace
    const bool swap = tr1->data()->size() > tr2->data()->size();
    GenericRecordCPtr trShorter = swap ? tr2 : tr1;
    GenericRecordCPtr trLonger  = swap ? tr1 : tr2; 

    const double *smpsS = DoubleArray::ConstCast(trShorter->data())->typedData();
    const double *smpsL = DoubleArray::ConstCast(trLonger->data())->typedData();
    const int smpsSsize = trShorter->data()->size();
    const int smpsLsize = trLonger->data()->size();

    //
    // for later quality check: save local maxima/minima
    //
    struct LocalMaxima {
        bool notDecreasing = false;
        double prevCoeff = -1;
        vector<double> values;
        void update(double coeff)
        {
            if ( ! std::isfinite(coeff) )
                return;

            if (coeff < prevCoeff && notDecreasing )
                values.push_back(prevCoeff);
            notDecreasing = coeff >= prevCoeff;
            prevCoeff = coeff;
        }
    };
    LocalMaxima localMaxs, localMins;

    for (int delay = -maxDelaySmps; delay < maxDelaySmps; delay++)
    {
        double numer = 0, denomL = 0, denomS = 0;
        for (int idxS = 0; idxS < smpsSsize; idxS++)
        {
            int idxL = idxS + (smpsLsize-smpsSsize)/2 + delay;
            if (idxL < 0 || idxL >= smpsLsize)
                continue;
            numer  += smpsS[idxS] * smpsL[idxL];
            denomL += smpsL[idxL] * smpsL[idxL];
            denomS += smpsS[idxS] * smpsS[idxS];
        }
        const double denom =  std::sqrt(denomS * denomL);
        const double coeff = numer / denom;
        if ( std::abs(coeff) > std::abs(coeffOut) || ! std::isfinite(coeffOut) )
        {
            coeffOut = coeff;
            delayOut = delay / freq; // samples to secs
        }

        // for later quality check
        localMaxs.update(coeff);
        localMins.update(-coeff);
    }

    if ( swap )
    {
        delayOut = -delayOut;
    }

    /*
     * To avoid errors introduced by cycle skipping the differential time measurement is
     * only accepted, if all side lobe maxima CCslm of the cross-correlation function 
     * fulfill the following condition:
     *
     *                CCslm < CCmax - ( (1.0-CCmax) / 2.0 )
     *
     * where CCmax corresponds to the global maximum of the cross-correlation function.
     * By discarding measurements with local maxima CCslm close to the global maximum CC,
     * the number of potential blunders due to cycle skipping is significantly reduced.
     *
     * See Diehl et al. (2017): The induced earthquake sequence related to the St. Gallen
     * deep geothermal project: Fault reactivation and fluid interactions imaged by
     * microseismicity
     * */
    if ( qualityCheck && std::isfinite(coeffOut))
    {
        double threshold = std::abs(coeffOut) - ( (1.0 - std::abs(coeffOut)) / 2.0 );
        int numMax = 0;
        vector<double> localMs = coeffOut > 0 ? localMaxs.values : localMins.values;
        for (double CCslm : localMs)
        {
            if (std::isfinite(CCslm) && CCslm >= threshold) numMax++;
            if (numMax > 1)
            {
                coeffOut = std::nan("");
                break;
            }
        }
    }

    if ( ! std::isfinite(coeffOut) )
    {
        coeffOut = 0;
        delayOut = 0.;
    }

    return true;
}



/*
 * Create absolute travel times difference file for pairs of earthquakes.
 */ 
void
HypoDD::createDtCt(map<unsigned,CatalogPtr>& neighbourCats, const string& dtctFile) const
{
    SEISCOMP_INFO("Creating differential travel time file %s", dtctFile.c_str());

    ofstream outStream(dtctFile);
    if ( !outStream.is_open() )
        throw runtime_error("Cannot create file " + dtctFile);

    for (const auto& kv : neighbourCats)
        writeAbsTTimePairs(kv.second, kv.first, outStream);
}


void HypoDD::createDtCt(CatalogPtr& catalog, unsigned evToRelocateId, const string& dtctFile) const
{
    map<unsigned,CatalogPtr> neighbourCats = { {evToRelocateId, catalog} };
    createDtCt(neighbourCats, dtctFile);
}


/* 
 * Create absolute travel times file (dt.ct) for hypodd
 *
 * Each event pair is listed by a header line (in free format)
 * #, ID1, ID2
 * followed by nobs lines of observations (in free format):
 * STA, TT1, TT2, WGHT, PHA
 * 
 */
void HypoDD::writeAbsTTimePairs(const CatalogCPtr& catalog,
                                unsigned evToRelocateId,
                                ofstream& outStream) const
{
    auto search = catalog->getEvents().find(evToRelocateId);
    if (search == catalog->getEvents().end())
    {
        string msg = stringify("Cannot find event id %u in the catalog.", evToRelocateId);
        throw runtime_error(msg);
    }
    const Catalog::Event& refEv = search->second;

    // loop through catalog events
    for (const auto& kv : catalog->getEvents() )
    {
        const Catalog::Event& event = kv.second;

        if (event == refEv)
            continue;

        int dtCount = 0;
        stringstream evStream;
        evStream << stringify("# %10u %10u", refEv.id, event.id) << endl;

        // loop through event phases
        auto eqlrng = catalog->getPhases().equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
            const Catalog::Phase& phase = it->second;

            // fetch phase pair from reference event and compute absolute travel time difference
            auto itRef = catalog->searchPhase(refEv.id, phase.stationId, phase.procInfo.type);

            if ( itRef != catalog->getPhases().end() )
            { 
                const Catalog::Phase& refPhase = itRef->second;

                double ref_travel_time = refPhase.time - refEv.time;
                if (ref_travel_time < 0)
                {
                    SEISCOMP_DEBUG("Ignoring phase '%s' with negative travel time (event '%s')",
                                   string(refPhase).c_str(), string(refEv).c_str());
                    continue;
                }
                double travel_time = phase.time - event.time;
                if (travel_time < 0)
                {
                    SEISCOMP_DEBUG("Ignoring phase '%s' with negative travel time (event '%s')",
                                   string(phase).c_str(), string(event).c_str());
                    continue;
                }

                // get common observation weight for pair (FIXME: take the lower one? average?)
                double weight = (refPhase.procInfo.weight + phase.procInfo.weight) / 2.0;

                evStream << stringify("%-12s %.6f %.6f %.2f %s",
                                      refPhase.stationId.c_str(), ref_travel_time,
                                      travel_time, weight, refPhase.procInfo.type.c_str());
                evStream << endl;
                dtCount++;
            }
        }
        if (dtCount > 0)
            outStream << evStream.str();
    }
}



/*
 * Create differential travel times file from cross correlation
 * for pairs of earthquakes.
 */
void
HypoDD::createDtCc(map<unsigned,CatalogPtr>& neighbourCats,
                   const string& dtccFile,
                   const XCorrCache& xcorr)
{
    SEISCOMP_INFO("Creating Cross correlation differential travel time file %s", dtccFile.c_str());

    ofstream outStream(dtccFile);
    if ( !outStream.is_open() )
        throw runtime_error("Cannot create file " + dtccFile);

    for (auto& kv : neighbourCats)
    {
        unsigned evToRelocateId = kv.first;
        CatalogPtr& neighbourCat = kv.second;
        writeXcorrDiffTTimePairs(neighbourCat, evToRelocateId, xcorr, outStream);
    }
}


void HypoDD::createDtCc(CatalogPtr& catalog,
                        unsigned evToRelocateId,
                        const string& dtccFile,
                        const XCorrCache& xcorr)
{
    map<unsigned,CatalogPtr> neighbourCats = {{evToRelocateId,catalog}};
    createDtCc(neighbourCats, dtccFile, xcorr);
}

/*
 * Compute and store to file differential travel times from cross
 * correlation for pairs of earthquakes.
 *
 * Each event pair is listed by a header line (in free format)
 * #, ID1, ID2, OTC
 * followed by lines with observations (in free format):
 * STA, DT, WGHT, PHA
 *
 */
void HypoDD::writeXcorrDiffTTimePairs(CatalogPtr& catalog,
                                      unsigned refEvId,
                                      const XCorrCache& xcorr,
                                      ofstream& outStream)
{
    auto search = catalog->getEvents().find(refEvId);
    if (search == catalog->getEvents().end())
    {
        string msg = stringify("Cannot find event id %u in the catalog.", refEvId);
        throw runtime_error(msg);
    }
    const Catalog::Event& refEv = search->second;

    // loop through neighbouring events
    //
    for (const auto& kv : catalog->getEvents() )
    {
        const Catalog::Event& event = kv.second;

        if (event == refEv)
            continue;

        int dtCount = 0;
        stringstream evStream;
        evStream << stringify("# %10u %10u       0.0", refEv.id, event.id) << endl;

        // loop through event phases
        auto eqlrng = catalog->getPhases().equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
            const Catalog::Phase& phase = it->second;

            // fetch xcorr results for this event pair at station/phase
            if ( xcorr.has(refEv.id, event.id, phase.stationId, phase.procInfo.type) )
            {
                const auto& data = xcorr.get(refEv.id, event.id, phase.stationId, phase.procInfo.type);

                evStream << stringify("%-12s %.6f %.4f %s", phase.stationId.c_str(),
                                      data.dtcc, data.weight, phase.procInfo.type.c_str());
                evStream << endl;
                dtCount++;
            }
        }
        if ( dtCount > 0 )
            outStream << evStream.str();
    }
}



/*
 * Reads the event pairs matched in dt.ct which are selected by ph2dt and
 * calculate cross correlated differential travel_times for every pair.
 * input dt.ct 
 * output dt.cc
 */
set<unsigned>
HypoDD::createDtCcPh2dt(const CatalogCPtr& catalog, const string& dtctFile, const string& dtccFile)
{
    SEISCOMP_INFO("Creating Cross correlation differential travel time file %s from input file %s",
                   dtccFile.c_str(), dtctFile.c_str());

    if ( ! Util::fileExists(dtctFile) )
        throw runtime_error("Unable to perform cross correlation, cannot find file: " + dtctFile);

    ofstream outStream(dtccFile);
    if ( !outStream.is_open() )
        throw runtime_error("Cannot create file " + dtccFile);

    _counters = {0};
    _wf->resetCounters();

    set<unsigned> selectedEvents;

    const std::map<unsigned,Catalog::Event>& events = catalog->getEvents();
    const Catalog::Event *ev1 = nullptr, *ev2 = nullptr;
    int dtCount = 0;
    stringstream evStream;

    // read file one line a time
    ifstream in(dtctFile);
    while (!in.eof())
    {
        string row;
        std::getline(in, row);
        if (in.bad() || in.fail())
            break;

        // split line on space
        static const std::regex regex(R"([\s]+)", std::regex::optimize);
        std::sregex_token_iterator it{row.begin(), row.end(), regex, -1};
        std::vector<std::string> fields{it, {}};

        // remove the first empty element if the line start with spaces
        if ( !fields.empty() && fields[0] == "")
            fields.erase(fields.begin());

        // check beginning of a new event pair line (# ID1 ID2)
        if (fields[0] == "#" && fields.size() == 3)
        {
            unsigned evId1 = std::stoul(fields[1]);
            unsigned evId2 = std::stoul(fields[2]);
            auto search1 = events.find(evId1);
            auto search2 = events.find(evId2);
            if (search1 == events.end() || search2 == events.end())
            {
                string msg = stringify("Internal logic error: file %s contains events ids (%s or %s) "
                                       "that are not present in the input catalog.",
                                       dtctFile.c_str(), string(*ev1).c_str(), string(*ev2).c_str());
                throw runtime_error(msg.c_str());
            }
            ev1 = &search1->second;
            ev2 = &search2->second;

            selectedEvents.insert(evId1);
            selectedEvents.insert(evId2);

            // write the pairs has been built up to now
            if (dtCount > 0 )
                outStream << evStream.str();
            evStream.str("");
            evStream.clear();
            dtCount = 0;

            evStream << stringify("# %10u %10u       0.0", ev1->id, ev2->id) << endl;
        }
        // observation line (STA, TT1, TT2, WGHT, PHA)
        else if(ev1 != nullptr && ev2 != nullptr && fields.size() == 5)
        {
            string stationId = fields[0];
            string phaseType = fields[4];

            // loop through event 1 phases
            auto eqlrng = catalog->getPhases().equal_range(ev1->id);
            for (auto it = eqlrng.first; it != eqlrng.second; ++it)
            {
                const Catalog::Phase& phase1 = it->second;
                if (phase1.stationId == stationId &&
                    phase1.procInfo.type == phaseType )
                {
                    // loop through event 2 phases
                    eqlrng = catalog->getPhases().equal_range(ev2->id);
                    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
                    {
                        const Catalog::Phase& phase2 = it->second;
                        if (phase2.stationId == stationId &&
                            phase2.procInfo.type == phaseType )
                        {
                            PhaseXCorrCfg phCfg {_useCatalogDiskCache, &_wfCache, true};
                            double coeff, lag, dtcc, weight;

                            if ( xcorrPhases(*ev1, phase1, phCfg, *ev2, phase2, phCfg,
                                             coeff, lag, dtcc, weight) )
                            {
                                evStream << stringify("%-12s %.6f %.4f %s", stationId.c_str(),
                                                      dtcc, weight, phaseType.c_str());
                                evStream << endl;
                                dtCount++;
                            }
                            break;
                        }
                    }
                    break;
                }
            }
        }
        else
        {
            ev1 = ev2 = nullptr;
            SEISCOMP_WARNING("Skipping unrecognized line from '%s' (line='%s')",
                           dtctFile.c_str(), row.c_str());
        }
    }

    if (dtCount > 0 )
        outStream << evStream.str();

    printCounters();

    return selectedEvents;
}


/*
 *  Write the station.dat input file for ph2dt and hypodd
 *  One station per line:
 *  STA, LAT, LON, ELV, MODID
 *
 *  E.g.
 NCAAS 38.4301 -121.11   12
 NCABA 38.8793 -121.067  25
 NCABJ 39.1658 -121.193  35
 NCABR 39.1381 -121.48   14
 *
 */
void HypoDD::createStationDatFile(const CatalogCPtr& catalog, const string& staFileName) const
{
    SEISCOMP_INFO("Creating station file %s", staFileName.c_str());

    ofstream outStream(staFileName);
    if ( !outStream.is_open() ) {
        string msg = "Cannot create file " + staFileName;
        throw runtime_error(msg);
    }

    for (const auto& kv :  catalog->getStations() )
    {
        const Catalog::Station& station = kv.second;
        outStream << stringify("%-12s %12.6f %12.6f %12.f",
                              station.id.c_str(), station.latitude,
                              station.longitude, station.elevation);
        outStream << endl;
    }
}


/* Write the phase.dat input file for ph2dt
 * ph2dt accepts hypocenter, followed by its travel time data in the following format:
 * #, YR, MO, DY, HR, MN, SC, LAT, LON, DEP, MAG, EH, EZ, RMS, ID
 * followed by nobs lines of observations:
 * STA, TT, WGHT, PHA
 * e.g.
 *
 #  1985  1 24  2 19 58.71  37.8832 -122.2415    9.80 1.40 0.2 0.5 0.0    38542
 NCCSP       2.850  -1.000   P
 NCCSP       2.910   0.016   P
 NCCBW       3.430  -1.000   P
 NCCBW       3.480   0.031   P
 #  1996 11  9  7  8 36.70  37.8810 -122.2457    9.14 1.10 0.6 0.6 0.0   484120
 NCCVP       1.860  -1.000   P
 NCCSP       2.770   0.250   P
 NCCMC       2.810   0.125   P
 NCCBW       3.360   0.250   P
 *
 */
void HypoDD::createPhaseDatFile(const CatalogCPtr& catalog, const string& phaseFileName) const
{
    SEISCOMP_INFO("Creating phase file %s", phaseFileName.c_str());

    ofstream outStream(phaseFileName);
    if ( !outStream.is_open() ) {
        string msg = "Cannot create file " + phaseFileName;
        throw runtime_error(msg);
    }

    for (const auto& kv :  catalog->getEvents() )
    {
        const Catalog::Event& event = kv.second;

        int year, month, day, hour, min, sec, usec;
        if ( ! event.time.get(&year, &month, &day, &hour, &min, &sec, &usec) )
        {
            SEISCOMP_WARNING("Cannot convert origin time for event '%s'", string(event).c_str());
            continue;
        }

        outStream << stringify("# %d %d %d %d %d %.2f %.6f %.6f %.3f %.2f %.4f %.4f %.4f %u",
                              year, month, day, hour, min, sec + double(usec)/1.e6,
                              event.latitude,event.longitude,event.depth,
                              event.magnitude, event.horiz_err, event.vert_err,
                              event.rms, event.id);
        outStream << endl;

        auto eqlrng = catalog->getPhases().equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
            const Catalog::Phase& phase = it->second;

            double travel_time = phase.time - event.time;
            if (travel_time < 0)
            {
                SEISCOMP_DEBUG("Ignoring phase '%s' with negative travel time (event '%s')",
                               string(phase).c_str(), string(event).c_str());
                continue; 
            }

            outStream << stringify("%-12s %12.6f %5.2f %4s",
                                  phase.stationId.c_str(), travel_time,
                                  phase.procInfo.weight, phase.procInfo.type.c_str());
            outStream << endl;
        }
    }
}


/* Write the event.dat input file for hypodd
 * One event per line:
 * DATE, TIME, LAT, LON, DEP, MAG, EH, EV, RMS, ID
 * e.g.
 *
19850124   2195871   37.8832  -122.2415      9.800   1.4    0.15    0.51   0.02      38542
19911126  14274555   37.8738  -122.2432      9.950   1.4    0.22    0.53   0.09     238298
19861019  20503808   37.8802  -122.2405      9.370   1.6    0.16    0.50   0.05      86036
19850814  18015544   37.8828  -122.2497      7.940   2.0    0.13    0.45   0.06      52942
19850527    430907   37.8778  -122.2412      9.050   1.1    0.26    0.75   0.03      48565
19850402   5571645   37.8825  -122.2420      9.440   1.9    0.12    0.30   0.04      45165
 *
 */
void HypoDD::createEventDatFile(const CatalogCPtr& catalog, const string& eventFileName) const
{
    SEISCOMP_INFO("Creating event file %s", eventFileName.c_str());

    ofstream outStream(eventFileName);
    if ( !outStream.is_open() )
    {
        string msg = "Cannot create file " + eventFileName;
        throw runtime_error(msg);
    }

    for (const auto& kv :  catalog->getEvents() )
    {
        const Catalog::Event& event = kv.second;

        int year, month, day, hour, min, sec, usec;
        if ( ! event.time.get(&year, &month, &day, &hour, &min, &sec, &usec) )
        {
            SEISCOMP_WARNING("Cannot convert origin time for event '%s'", string(event).c_str());
            continue;
        }

        outStream << stringify("%d%02d%02d  %02d%02d%04d %.6f %.6f %.3f %.2f %.4f %.4f %.4f %u",
                              year, month, day, hour, min, int(sec * 1e2 + usec / 1e4),
                              event.latitude, event.longitude, event.depth,
                              event.magnitude, event.horiz_err, event.vert_err,
                              event.rms, event.id);
        outStream << endl;
    }
}



/*
 * run ph2dt
 * input files: ph2dt.inp station.dat phase.dat
 * output files: station.sel event.sel event.dat dt.ct
 */
void HypoDD::runPh2dt(const string& workingDir, const string& stationFile, const string& phaseFile) const
{
    SEISCOMP_INFO("Running ph2dt...");

    if ( !Util::fileExists(stationFile) )
        throw runtime_error("Unable to run ph2dt, file doesn't exist: " + stationFile);

    if ( !Util::fileExists(phaseFile) )
        throw runtime_error("Unable to run ph2dt, file doesn't exist: " + phaseFile);

    if ( !Util::fileExists(_cfg.ph2dt.ctrlFile) )
        throw runtime_error("Unable to run ph2dt, control file doesn't exist: " + _cfg.ph2dt.ctrlFile);

    // copy control file while replacing input/output file names
    map<int,string> linesToReplace = {
        {1, boost::filesystem::path(stationFile).filename().string()},// requires boost 1.60 boost::filesystem::path(stationFile).lexically_relative(workingDir).string()},
        {2, boost::filesystem::path(phaseFile).filename().string()},  // requires boost 1.60 boost::filesystem::path(phaseFile).lexically_relative(workingDir).string()},
    };
    copyFileAndReplaceLines(_cfg.ph2dt.ctrlFile,
                            (boost::filesystem::path(workingDir)/"ph2dt.inp").string(),
                            linesToReplace);

    // run ph2dt (use /bin/sh to get stdout/strerr redirection)
    string cmd = stringify("%s %s >ph2dt.out 2>&1",
                           _cfg.ph2dt.exec.c_str(), "ph2dt.inp");
    ::startExternalProcess({"/bin/sh", "-c", cmd}, true, workingDir);
}



/*
 * run hypodd executable
 * input files: dt.cc dt.ct event.sel station.sel hypoDD.inp
 * output files: hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
 */
void HypoDD::runHypodd(const string& workingDir, const string& dtccFile, const string& dtctFile,
                       const string& eventFile, const string& stationFile, const std::string& ctrlFile) const
{
    SEISCOMP_INFO("Running hypodd...");

    if ( !Util::fileExists(dtccFile) )
        throw runtime_error("Unable to run hypodd, file doesn't exist: " + dtccFile);

    if ( !Util::fileExists(dtctFile) )
        throw runtime_error("Unable to run hypodd, file doesn't exist: " + dtctFile);

    if ( !Util::fileExists(eventFile) )
        throw runtime_error("Unable to run hypodd, file doesn't exist: " + eventFile);

    if ( !Util::fileExists(stationFile) )
        throw runtime_error("Unable to run hypodd, file doesn't exist: " + stationFile);

    if ( !Util::fileExists(ctrlFile) )
        throw runtime_error("Unable to run hypodd, control file doesn't exist: " + ctrlFile);

    // check if hypodd.inp is for version 2.1
    ifstream ctrlFileStrm(ctrlFile);
    if ( ! ctrlFileStrm.is_open() )
    {
        string msg = stringify("Cannot open hypodd control file %s", ctrlFile.c_str());
        throw runtime_error(msg);
    }

    int lineOffset = 0;
    string line;
    if( std::getline(ctrlFileStrm, line) && line == "hypoDD_2")
        lineOffset = 1;

    // copy control file while replacing input/output file names
    map<int,string> linesToReplace = {
        {lineOffset + 1, boost::filesystem::path(dtccFile   ).filename().string()}, // requires boost 1.60 boost::filesystem::path(dtccFile  ).lexically_relative(workingDir).string() },
        {lineOffset + 2, boost::filesystem::path(dtctFile   ).filename().string()}, // requires boost 1.60 boost::filesystem::path(dtctFile   ).lexically_relative(workingDir).string() },
        {lineOffset + 3, boost::filesystem::path(eventFile  ).filename().string()}, // requires boost 1.60 boost::filesystem::path(eventFile  ).lexically_relative(workingDir).string() },
        {lineOffset + 4, boost::filesystem::path(stationFile).filename().string()}, // requires boost 1.60 boost::filesystem::path(stationFile).lexically_relative(workingDir).string() },
        {lineOffset + 5, "hypoDD.loc"},
        {lineOffset + 6, "hypoDD.reloc"},
        {lineOffset + 7, "hypoDD.sta"},
        {lineOffset + 8, "hypoDD.res"},
        {lineOffset + 9, "hypoDD.src"}
    };
    copyFileAndReplaceLines(ctrlFile, (boost::filesystem::path(workingDir)/"hypoDD.inp").string(),
                            linesToReplace);

    // run Hypodd (use /bin/sh to get stdout/strerr redirection)
    string cmd = stringify("%s %s >hypoDD.out 2>&1", _cfg.hypodd.exec.c_str(), "hypoDD.inp");
    ::startExternalProcess({"/bin/sh", "-c", cmd}, true, workingDir);
}


/*
 * load a catalog from hypodd output file
 * input: hypoDD.reloc
 *
 * One event per line (written in fixed, but may be read in free format):
 * ID, LAT, LON, DEPTH, X, Y, Z, EX, EY, EZ, YR, MO, DY, HR, MI, SC, MAG, NCCP, NCCS, NCTP,
NCTS, RCC, RCT, CID
 *
 */
CatalogPtr HypoDD::loadRelocatedCatalog(const CatalogCPtr& originalCatalog,
                                        const std::string& ddrelocFile,
                                        const std::string& ddresidualFile) const
{
    SEISCOMP_INFO("Loading catalog relocated by hypodd...");

    if ( !Util::fileExists(ddrelocFile) )
        throw runtime_error("Cannot load hypodd relocated catalog file: " + ddrelocFile);

    map<string,Catalog::Station> stations = originalCatalog->getStations();
    map<unsigned,Catalog::Event> events = originalCatalog->getEvents();
    multimap<unsigned,Catalog::Phase> phases = originalCatalog->getPhases();

    // read relocation file one line a time
    ifstream in(ddrelocFile);
    while (!in.eof())
    {
        string row;
        std::getline(in, row);
        if (in.bad() || in.fail())
            break;

        // split line on space
        static const std::regex regex(R"([\s]+)", std::regex::optimize);
        std::sregex_token_iterator it{row.begin(), row.end(), regex, -1};
        std::vector<std::string> fields{it, {}};

        // remove the first empty element if the line start with spaces
        if ( !fields.empty() && fields[0] == "")
            fields.erase(fields.begin());

        if (fields.size() != 24)
        {
            SEISCOMP_WARNING("Skipping unrecognized line from '%s' (line='%s')",
                             ddrelocFile.c_str(), row.c_str());
            continue;
        }

        // load corresponding event and update information
        unsigned eventId = std::stoul(fields[0]);
        auto search = events.find(eventId);
        if (search == events.end())
        {
            // skip events that are not part of the passed catalog
            continue;
        }
        Catalog::Event& event = search->second;
        event.latitude  = std::stod(fields[1]);
        event.longitude = std::stod(fields[2]);
        event.depth     = std::stod(fields[3]);

        int year  = std::stoi(fields[10]);
        int month = std::stoi(fields[11]);
        int day   = std::stoi(fields[12]);
        int hour  = std::stoi(fields[13]);
        int min   = std::stoi(fields[14]);
        double seconds = std::stod(fields[15]);
        int sec  = int(seconds);
        int usec = (seconds - sec) * 1.e6;

        event.time = Core::Time(year, month, day, hour, min, sec, usec);

        event.relocInfo.isRelocated = true;
        event.relocInfo.lonUncertainty   = std::stod(fields[7])/1000.;
        event.relocInfo.latUncertainty   = std::stod(fields[8])/1000.;
        event.relocInfo.depthUncertainty = std::stod(fields[9])/1000.;
        event.relocInfo.numCCp = std::stoi(fields[17]);
        event.relocInfo.numCCs = std::stoi(fields[18]);
        event.relocInfo.numCTp = std::stoi(fields[19]);
        event.relocInfo.numCTs = std::stoi(fields[20]);
        event.relocInfo.rmsResidualCC = std::stod(fields[21]);
        event.relocInfo.rmsResidualCT = std::stod(fields[22]);
        if  ( (event.relocInfo.numCTp +  event.relocInfo.numCTs) > 0 &&
              (event.relocInfo.numCCp +  event.relocInfo.numCCs) > 0   )
            event.rms = (event.relocInfo.rmsResidualCC + event.relocInfo.rmsResidualCT) / 2.;
        else if ( (event.relocInfo.numCTp +  event.relocInfo.numCTs) > 0)
            event.rms = event.relocInfo.rmsResidualCT;
        else if ( (event.relocInfo.numCCp +  event.relocInfo.numCCs) > 0)
            event.rms = event.relocInfo.rmsResidualCC;
        else
            event.rms = 0;

    }

    // read residual file one line a time to fetch residuals and final weights
    if ( ! ddresidualFile.empty() )
    {
        struct residual {
            double residuals = 0;
            double weights = 0;
            int count = 0;
        };
        map<string,struct residual> resInfos;

        // 1=ccP; 2=ccS; 3=ctP; 4=ctS
        map<string, string> dataTypeMap = { {"1","P"}, {"2","S"}, {"3","P"}, {"4","S"} };

        ifstream in(ddresidualFile);
        while (!in.eof())
        {
            string row;
            std::getline(in, row);
            if (in.bad() || in.fail())
                break;

            // split line on space
            static const std::regex regex(R"([\s]+)", std::regex::optimize);
            std::sregex_token_iterator it{row.begin(), row.end(), regex, -1};
            std::vector<std::string> fields{it, {}};

            // remove the first empty element if the line start with spaces
            if ( !fields.empty() && fields[0] == "")
                fields.erase(fields.begin());

            if (fields.size() != 9)
            {
                SEISCOMP_WARNING("Skipping unrecognized line from '%s' (line='%s')",
                                 ddresidualFile.c_str(), row.c_str());
                continue;
            }

            string stationId = fields[0];
            unsigned ev1Id = std::stoul(fields[2]);
            unsigned ev2Id = std::stoul(fields[3]);
            string dataType = dataTypeMap[ fields[4] ]; // 1=ccP; 2=ccS; 3=ctP; 4=ctS
            double residual = std::stod(fields[6]) / 1000.; //ms -> s
            double finalWeight = std::stod(fields[7]);

            string key1 = to_string(ev1Id) + "+" + stationId + "+" + dataType;
            struct residual& info1 = resInfos[key1];
            info1.residuals += residual;
            info1.weights += finalWeight;
            info1.count++;

            string key2 = to_string(ev2Id) + "+" + stationId + "+" + dataType;
            struct residual& info2 = resInfos[key2];
            info2.residuals += residual;
            info2.weights += finalWeight;
            info2.count++;
        }

        for (auto& pair : phases)
        {
            Catalog::Phase &phase = pair.second;
            string key = to_string(phase.eventId) + "+" + phase.stationId + "+" + phase.procInfo.type;
            if ( resInfos.find(key) != resInfos.end() )
            {
                struct residual& info = resInfos[key];
                phase.relocInfo.isRelocated = true;
                phase.relocInfo.residual = info.residuals / info.count;
                phase.relocInfo.finalWeight = info.weights / info.count;
            }
        }
    }

    return new Catalog(stations, events, phases);
}


} // HDD
} // Seiscomp


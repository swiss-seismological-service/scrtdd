/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 *   You can redistribute and/or modify this program under the             *
 *   terms of the "SED Public License for Seiscomp Contributions"          *
 *                                                                         *
 *   You should have received a copy of the "SED Public License for        *
 *   Seiscomp Contributions" with this. If not, you can find it at         *
 *   http://www.seismo.ethz.ch/static/seiscomp_contrib/license.txt         *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   "SED Public License for Seiscomp Contributions" for more details      *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#include "hypodd.h"
#include "csvreader.h"

#include <seiscomp3/core/datetime.h>
#include <seiscomp3/core/strings.h>
#include <seiscomp3/core/typedarray.h>
#include <seiscomp3/io/recordinput.h>
#include <seiscomp3/io/records/mseedrecord.h>
#include <seiscomp3/processing/operator/transformation.h>
#include <seiscomp3/processing/operator/ncomps.h>
#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/filter/butterworth.h>
#include <seiscomp3/client/inventory.h>
#include <seiscomp3/datamodel/station.h>
#include <seiscomp3/datamodel/event.h>
#include <seiscomp3/datamodel/pick.h>
#include <seiscomp3/datamodel/origin.h>
#include <seiscomp3/datamodel/magnitude.h>
#include <seiscomp3/datamodel/amplitude.h>
#include <seiscomp3/utils/files.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
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

    SEISCOMP_INFO("Executing command: %s", cmdline.c_str());

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



template <class Map, class Val> typename Map::iterator
searchByValue(Map & SearchMap, const Val & SearchVal)
{
    typename Map::iterator iRet = SearchMap.end();
    for (typename Map::iterator iTer = SearchMap.begin(); iTer != SearchMap.end(); iTer ++)
    {
        if (iTer->second == SearchVal)
        {
            iRet = iTer;
            break;
        }
    }
    return iRet;
}



template <class Map, class Val> typename Map::const_iterator
searchByValue(const Map & SearchMap, const Val & SearchVal)
{
    typename Map::const_iterator iRet = SearchMap.end();
    for (typename Map::const_iterator iTer = SearchMap.begin(); iTer != SearchMap.end(); iTer ++)
    {
        if (iTer->second == SearchVal)
        {
            iRet = iTer;
            break;
        }
    }
    return iRet;
}



template <class T>
T nextPowerOf2(T a, T min=1, T max=1<<31)
{
    int b = min;
    while (b < a)
    {
        b <<= 1;
        if (b > max)
            return -1;
    }
    return b;
}



void writeTrace(GenericRecordPtr trace, string file)
{
    std::ofstream ofs(file);
    IO::MSeedRecord msRec(*trace);
    int reclen = msRec.data()->size()*msRec.data()->bytes() + 64;
    reclen = nextPowerOf2<int>(reclen, 128, 1048576); // MINRECLEN 128, MAXRECLEN 1048576
    if (reclen > 0)
    {
        msRec.setOutputRecordLength(reclen);
        msRec.write(ofs);
    }
}



GenericRecordPtr readTrace(string file)
{
    std::ifstream ifs(file);
    IO::MSeedRecord msRec(Array::DOUBLE, Record::Hint::DATA_ONLY);
    msRec.read(ifs);
    GenericRecordPtr trace = new GenericRecord(msRec);
    trace->setData(msRec.data()->clone()); // copy data too
    return trace;
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



DataModel::Station* findStation(const string& netCode,
                                const string& stationCode,
                                const Core::Time& atTime)
{
    DataModel::Inventory *inv = Client::Inventory::Instance()->inventory();
    if ( ! inv )
    {
        SEISCOMP_ERROR("Inventory not available");
        return nullptr;
    }

    for (size_t n = 0; n < inv->networkCount(); ++n )
    {
        DataModel::Network *net = inv->network(n);
        if ( net->start() > atTime ) continue;
        try { if ( net->end() < atTime ) continue; }
        catch ( ... ) {}

        if ( !Core::wildcmp(netCode, net->code()) ) continue;

        for ( size_t s = 0; s < net->stationCount(); ++s )
        {
            DataModel::Station *sta = net->station(s);
            if ( sta->start() > atTime ) continue;
            try { if ( sta->end() < atTime ) continue; }
            catch ( ... ) {}

            if ( !Core::wildcmp(stationCode, sta->code()) ) continue;

            return sta;
        }
    }

    return nullptr;
}



DataModel::SensorLocation *findSensorLocation(DataModel::Station *station,
                                              const std::string &code,
                                              const Core::Time &atTime)
{
    if ( station )
    {
        for ( size_t i = 0; i < station->sensorLocationCount(); ++i )
        {
            DataModel::SensorLocation *loc = station->sensorLocation(i);

            try {
                if ( loc->end() <= atTime ) continue;
            }
            catch ( Seiscomp::Core::ValueException& ) {}

            if ( loc->start() > atTime ) continue;

            if ( loc->code() == code )
                return loc;
        }
    }
    return nullptr;
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


}


namespace Seiscomp {
namespace HDD {



DataModel::PublicObject* DataSource::getObject(const Core::RTTI& classType,
                                                const std::string& publicID)
{
    DataModel::PublicObject* ret = nullptr;

    if ( _eventParameters && ! ret)
    {
        if (classType == DataModel::Pick::TypeInfo() )
            ret = _eventParameters->findPick(publicID);
        else if (classType == DataModel::Amplitude::TypeInfo() )
            ret = _eventParameters->findAmplitude(publicID);
        else if (classType == DataModel::Origin::TypeInfo() )
            ret = _eventParameters->findOrigin(publicID);
        else if (classType == DataModel::Event::TypeInfo() )
            ret = _eventParameters->findEvent(publicID);
    }

    if ( _cache && ! ret)
    {
        ret = _cache->find(classType, publicID);
    }

   return ret;
}



void DataSource::loadArrivals(DataModel::Origin* org)
{
    bool found = false;

    if ( _eventParameters && ! found)
    {
        DataModel::Origin* epOrg = _eventParameters->findOrigin(org->publicID());
        if ( epOrg )
        {
            found = true;
            for (size_t i = 0; i < epOrg->arrivalCount(); i++)
                org->add(DataModel::Arrival::Cast(epOrg->arrival(i)->clone()));
        }
    }

    if ( _query && ! found )
    {
        found = true;
        _query->loadArrivals(org);
    }
}



DataModel::Event* DataSource::getParentEvent(const std::string& originID)
{
    DataModel::Event* ret = nullptr;

    if ( _eventParameters && ! ret )
    {
        for (size_t i = 0; i <  _eventParameters->eventCount() && !ret; i++)
        {
            DataModel::Event* ev = _eventParameters->event(i);
            for (size_t j = 0; j < ev->originReferenceCount() && !ret; j++)
            {
                DataModel::OriginReference* orgRef = ev->originReference(j);
                if (orgRef->originID() == originID)
                    ret = ev;
            }
        }
    }

    if ( _query && ! ret)
    {
        ret = _query->getEvent(originID);
    }

    return ret;
}



Catalog::Catalog() : Catalog(map<string,Station>(),
                             map<unsigned,Event>(),
                             multimap<unsigned,Phase>())
{
}



Catalog::Catalog(const map<string,Station>& stations,
                 const map<unsigned,Event>& events,
                 const multimap<unsigned,Phase>& phases)
{
    _stations = stations;
    _events = events;
    _phases = phases;
}



Catalog::Catalog(const string& stationFile, const string& eventFile, const string& phaFile)
{
    if ( !Util::fileExists(stationFile) )
    {
        string msg = "File " + stationFile + " does not exist";
        throw runtime_error(msg);
    }

    if ( !Util::fileExists(eventFile) )
    {
        string msg = "File " + eventFile + " does not exist";
        throw runtime_error(msg);
    }

    if ( !Util::fileExists(phaFile) )
    {
        string msg = "File " + phaFile + " does not exist";
        throw runtime_error(msg);
    }

    vector<map<string,string> > stations = CSV::readWithHeader(stationFile);

    for (const auto& row : stations )
    {
        Station sta;
        sta.id = row.at("id");
        sta.latitude = std::stod(row.at("latitude"));
        sta.longitude = std::stod(row.at("longitude"));
        sta.elevation = std::stod(row.at("elevation"));
        sta.networkCode = row.at("networkCode");
        sta.stationCode = row.at("stationCode");
        _stations[sta.id] = sta;
}

    vector<map<string,string> >events = CSV::readWithHeader(eventFile);

    for (const auto& row : events )
    {
        Event ev;
        ev.id          = std::stoul(row.at("id"));
        ev.time        = Core::Time::FromString(row.at("isotime").c_str(), "%FT%T.%fZ"); //iso format
        ev.latitude    = std::stod(row.at("latitude"));
        ev.longitude   = std::stod(row.at("longitude"));
        ev.depth       = std::stod(row.at("depth"));
        ev.magnitude   = std::stod(row.at("magnitude"));
        ev.horiz_err   = std::stod(row.at("horizontal_err"));
        ev.vert_err    = std::stod(row.at("vertical_err"));
        ev.rms         = std::stod(row.at("rms"));
        _events[ev.id] = ev;
    }

    vector<map<string,string> >phases = CSV::readWithHeader(phaFile);

    for (const auto& row : phases )
    {
        Phase ph;
        ph.eventId    = std::stoul(row.at("eventId"));
        ph.stationId  = row.at("stationId");
        ph.time        = Core::Time::FromString(row.at("isotime").c_str(), "%FT%T.%fZ"); //iso format
        ph.weight      = std::stod(row.at("weight"));
        ph.type        = row.at("type");
        ph.networkCode   = row.at("networkCode");
        ph.stationCode   = row.at("stationCode");
        ph.locationCode  = row.at("locationCode");
        ph.channelCode   = row.at("channelCode");
        _phases.emplace(ph.eventId, ph);
    }
}



void Catalog::add(const std::vector<DataModel::Origin*>& origins,
                  DataSource& dataSrc)
{
    for(DataModel::Origin* org : origins)
    {
        if ( org->arrivalCount() == 0)
            dataSrc.loadArrivals(org);

        // Add event
        Event ev;
        ev.id          = 0;
        ev.time        = org->time().value();
        ev.latitude    = org->latitude();
        ev.longitude   = org->longitude();
        ev.depth       = org->depth(); // km
        ev.horiz_err   = 0;
        ev.vert_err    = 0;
        ev.rms         = org->quality().standardError();
        DataModel::MagnitudePtr mag;
        // try to fetch preferred magnitude stored in the event
        DataModel::EventPtr parentEvent = dataSrc.getParentEvent(org->publicID());
        if ( parentEvent )
        {
            mag = dataSrc.get<DataModel::Magnitude>(parentEvent->preferredMagnitudeID());
        }
        if ( mag )
        {
            ev.magnitude = mag->magnitude();
        }
        else
        {
            SEISCOMP_WARNING("Origin %s: cannot load preferred magnitude from parent event, set it to 0",
                             org->publicID().c_str());
            ev.magnitude = 0.;
        }

        SEISCOMP_DEBUG("Adding origin '%s' to the catalog", org->publicID().c_str());

        // make sure we didn't ger duplicated events
        if ( this->searchEvent(ev) != _events.end())
        {
            SEISCOMP_WARNING("Skipping duplicated origin %s (%s)",
                             org->publicID().c_str(), string(ev).c_str());
            continue;
        }

        this->addEvent(ev, false);
        Event& newEvent = searchByValue(this->_events, ev)->second;  // fetch the id

        // Add Phases
        for ( size_t i = 0; i < org->arrivalCount(); ++i )
        {
            DataModel::Arrival *orgArr = org->arrival(i);
            const DataModel::Phase& orgPh = orgArr->phase();

            DataModel::PickPtr pick = dataSrc.get<DataModel::Pick>(orgArr->pickID());
            if ( !pick )
            {
                SEISCOMP_ERROR("Cannot load pick '%s' (origin %s)",
                               orgArr->pickID().c_str(), org->publicID().c_str());
                continue;
            }

            // find the station
            Station sta;
            sta.networkCode = pick->waveformID().networkCode();
            sta.stationCode = pick->waveformID().stationCode();

            // add station if not already there
            if (searchStation(sta) == _stations.end())
            {
                DataModel::Station* orgArrStation = findStation(sta.networkCode, sta.stationCode,
                                                                pick->time());
                if ( !orgArrStation )
                {
                    SEISCOMP_ERROR("Cannot load station %s.%s information for arrival '%s' (origin '%s')."
                                   "All picks associated with this station will not be used.",
                                    sta.networkCode.c_str(), sta.stationCode.c_str(),
                                   orgArr->pickID().c_str(), org->publicID().c_str());
                    continue;
                }
                sta.latitude = orgArrStation->latitude();
                sta.longitude = orgArrStation->longitude();
                sta.elevation = orgArrStation->elevation(); // meter
                this->addStation(sta, false);
            }
            // the station has to be there at this point
            sta = this->searchStation(sta)->second;

            Phase ph;
            ph.eventId    = newEvent.id;
            ph.stationId  = sta.id;
            ph.time       = pick->time().value();
            try {
                ph.weight = orgArr->weight();
            } catch ( Core::ValueException& ) {
                ph.weight = 1.;
                SEISCOMP_DEBUG("Pick '%s' (origin %s) has no weight set, use default weight of 1.",
                              orgArr->pickID().c_str(), org->publicID().c_str());
            }

            ph.type        = orgPh.code();
            ph.networkCode =  pick->waveformID().networkCode();
            ph.stationCode =  pick->waveformID().stationCode();
            ph.locationCode =  pick->waveformID().locationCode();
            ph.channelCode =  pick->waveformID().channelCode();
            this->addPhase(ph, false);
        }
    }
}



void Catalog::add(const std::vector<std::string>& ids, DataSource& dataSrc)
{
    vector<DataModel::Origin*> origins;

    for(const string& id : ids)
    {
        DataModel::OriginPtr org = dataSrc.get<DataModel::Origin>(id);
        if ( !org )
        {
            DataModel::EventPtr ev = dataSrc.get<DataModel::Event>(id);
            if ( ev )
            {
                org = dataSrc.get<DataModel::Origin>(ev->preferredOriginID());
            }
        }
        if ( !org )
        {
            SEISCOMP_ERROR("Cannot find origin/event with id %s", id.c_str());
            continue;
        }
        origins.push_back(org.get());
    }

    add(origins, dataSrc);
}



void Catalog::add(const std::string& idFile,  DataSource& dataSrc)
{
    if ( !Util::fileExists(idFile) )
    {
        string msg = "File " + idFile + " does not exist";
        throw runtime_error(msg);
    }

    vector<string> ids;
    vector< map<string,string> > rows = CSV::readWithHeader(idFile);

    for(const auto& row : rows)
    {
        const string& id = row.at("seiscompId");
        ids.push_back(id);
    }

    add(ids, dataSrc);
}



CatalogPtr Catalog::merge(const CatalogCPtr& other) const
{
    CatalogPtr mergedCatalog = new Catalog(getStations(), getEvents(), getPhases());

    for (const auto& kv :  other->getEvents() )
    {
        const Catalog::Event& event = kv.second;
        mergedCatalog->copyEvent(event, other, false);
    }

    return mergedCatalog;
}



bool Catalog::copyEvent(const Catalog::Event& event, const CatalogCPtr& other, bool keepEvId)
{
    if ( keepEvId )
    {
        if  (_events.find(event.id) != _events.end() )
            throw runtime_error("Cannot add event, internal logic error");
        _events[event.id] = event;
    }
    else
    {
        if ( ! addEvent(event, true) )
            throw runtime_error("Cannot add event, internal logic error");
    }

    unsigned newEventId = searchEvent(event)->first;

    auto eqlrng = other->getPhases().equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
        Catalog::Phase phase = it->second;

        auto search = other->getStations().find(phase.stationId);
        if (search == other->getStations().end())
        {
            string msg = stringify("Malformed catalog: cannot find station '%s' "
                                    " referenced by phase '%s'",
                                   phase.stationId.c_str(), string(phase).c_str());
            throw runtime_error(msg);
        }
        const Catalog::Station& station = search->second;
        addStation(station, true);
        string newStationId = searchStation(station)->first;

        phase.eventId = newEventId;
        phase.stationId = newStationId;

        addPhase(phase, true);
    }

    return true;
}



void Catalog::removeEvent(const Event& event)
{
    map<unsigned,Catalog::Event>::const_iterator it = searchEvent(event);
    if ( it != _events.end() )
    {
        _events.erase(it);
    }
    auto eqlrng = _phases.equal_range(event.id);
    _phases.erase(eqlrng.first, eqlrng.second);
}


void Catalog::removeEvent(unsigned eventId)
{
    map<unsigned,Catalog::Event>::const_iterator it = _events.find(eventId);
    if ( it != _events.end() )
    {
        _events.erase(it);
    }
    auto eqlrng = _phases.equal_range(eventId);
    _phases.erase(eqlrng.first, eqlrng.second);
}


map<string,Catalog::Station>::const_iterator Catalog::searchStation(const Station& station) const
{
    return searchByValue(_stations, station);
}


map<unsigned,Catalog::Event>::const_iterator Catalog::searchEvent(const Event& event) const
{
    return searchByValue(_events, event);
}


map<unsigned,Catalog::Phase>::const_iterator Catalog::searchPhase(const Phase& phase) const
{
    return searchByValue(_phases, phase);
}


bool Catalog::addStation(const Station& station, bool checkDuplicate)
{
    if (checkDuplicate && searchStation(station) != _stations.end())
    {
        return false;
    }
    Station newStation = station;
    newStation.id = newStation.networkCode + newStation.stationCode;
    _stations[newStation.id] = newStation;
    return true;
}


bool Catalog::addEvent(const Event& event, bool checkDuplicate)
{
    decltype(_events)::key_type maxKey = _events.empty() ? 0 : _events.rbegin()->first;
    if (checkDuplicate && searchEvent(event) != _events.end())
    {
        return false;
    }
    Event newEvent = event;
    newEvent.id = maxKey + 1;
    _events[newEvent.id] = newEvent;
    return true;
}


bool Catalog::addPhase(const Phase& phase, bool checkDuplicate)
{
    if (checkDuplicate && searchPhase(phase) != _phases.end())
    {
        return false;
    }
    _phases.emplace(phase.eventId, phase);
    return true;
}

void Catalog::writeToFile(string eventFile, string phaseFile, string stationFile) const
{
    stringstream evStreamNoReloc;
    stringstream evStreamReloc;

    evStreamNoReloc << "id,isotime,latitude,longitude,depth,magnitude,horizontal_err,vertical_err,rms"; 
    evStreamReloc << evStreamNoReloc.str() << ",relocated,lonUncertainty,latUncertainty,depthUncertainty,numCCp,numCCs,numCTp,numCTs,residualCC,residualCT" << endl;
    evStreamNoReloc << endl;

    bool relocInfo = false;
    for (const auto& kv : _events )
    {
        const Catalog::Event& ev = kv.second;

        stringstream evStream;
        evStream << stringify("%u,%s,%.6f,%.6f,%.4f,%.2f,%.4f,%.4f,%.4f",
                              ev.id,ev.time.iso().c_str(),
                              ev.latitude,ev.longitude,ev.depth,ev.magnitude,
                              ev.horiz_err,ev.vert_err,ev.rms);

        evStreamNoReloc << evStream.str() << endl;
        evStreamReloc   << evStream.str();

        if ( ! ev.relocInfo.isRelocated )
        {
            evStreamReloc << ",false,,,,,,,,,";
        }
        else
        {
            relocInfo = true;
            evStreamReloc <<  stringify(",true,%.6f,%.6f,%.6f,%d,%d,%d,%d,%.4f,%.4f",
                          ev.relocInfo.lonUncertainty, ev.relocInfo.latUncertainty,
                          ev.relocInfo.depthUncertainty, ev.relocInfo.numCCp,
                          ev.relocInfo.numCCs, ev.relocInfo.numCTp, ev.relocInfo.numCTs,
                          ev.relocInfo.rmsResidualCC, ev.relocInfo.rmsResidualCT);
        }
        evStreamReloc << endl;
    }

    ofstream evStream(eventFile);
    evStream << ( relocInfo ? evStreamReloc.str() : evStreamNoReloc.str() );

    ofstream phStream(phaseFile);
    phStream << "eventId,stationId,isotime,weight,type,networkCode,stationCode,locationCode,channelCode";
    if (relocInfo)
        phStream << ",usedInReloc,residual,finalWeight";
    phStream << endl;
    for (const auto& kv : _phases )
    {
        const Catalog::Phase& ph = kv.second;
        phStream << stringify("%u,%s,%s,%.2f,%s,%s,%s,%s,%s",
                              ph.eventId, ph.stationId.c_str(), ph.time.iso().c_str(),
                              ph.weight, ph.type.c_str(), ph.networkCode.c_str(),
                              ph.stationCode.c_str(), ph.locationCode.c_str(), ph.channelCode.c_str());

        if (relocInfo)
        {
            if ( ! ph.relocInfo.isRelocated )
            {
                phStream << ",false,,";
            }
            else
            {
                phStream << stringify(",true,%.4f,%.2f",
                                      ph.relocInfo.residual, ph.relocInfo.finalWeight);
            }
        }
        phStream << endl;
    }

    ofstream staStream(stationFile);
    staStream << "id,latitude,longitude,elevation,networkCode,stationCode" << endl;
    for (const auto& kv : _stations )
    {
        const Catalog::Station& sta = kv.second;
        staStream << stringify("%s,%.6f,%.6f,%.1f,%s,%s",
                               sta.id.c_str(), sta.latitude, sta.longitude, sta.elevation,
                               sta.networkCode.c_str(), sta.stationCode.c_str())
                  << endl;
    }
    
}


HypoDD::HypoDD(const CatalogCPtr& catalog, const Config& cfg, const string& workingDir)
{
    _cfg = cfg;
    _workingDir = workingDir;
    setCatalog(catalog);

    if ( !Util::pathExists(_workingDir) )
    {
        if ( !Util::createPath(_workingDir) )
        {
            string msg = "Unable to create working directory: " + _workingDir;
            throw runtime_error(msg);
        }
    }

    _cacheDir = (boost::filesystem::path(_workingDir)/"wfcache").string();
    if ( !Util::pathExists(_cacheDir) )
    {
        if ( !Util::createPath(_cacheDir) )
        {
            string msg = "Unable to create cache directory: " + _cacheDir;
            throw runtime_error(msg);
        }
    }
}



HypoDD::~HypoDD()
{
    if ( _workingDirCleanup )
    {
        // delete all expect the cache directory
        for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(_workingDir), {}))
        {       
            if ( ! boost::filesystem::equivalent(entry, _cacheDir) )
                boost::filesystem::remove_all(entry);
        }
    }
}



void HypoDD::setCatalog(const CatalogCPtr& catalog)
{
    _ddbgc = filterOutPhases(catalog, _cfg.validPphases, _cfg.validSphases);
}



// Creates dir name from event. This id has the following format:
// OriginTime_Lat_Lon_CreationDate
// eg 20111210115715_46343_007519_20111210115740
string HypoDD::generateWorkingSubDir(const Catalog::Event& ev) const
{
    char buf[20];

    string id;
    id = ev.time.toString("%Y%m%d%H%M%S");

    id += "_";

    // Latitude
    sprintf(buf, "%05d", int(ev.latitude*1000));
    id += buf;

    id += "_";

    // Longitude
    sprintf(buf, "%06d", int(ev.longitude*1000));
    id += buf;

    id += "_";

    Core::Time t = Core::Time::GMT();

    id += t.toString("%Y%m%d%H%M%S");

    return id;
}



void HypoDD::preloadData()
{

    for (const auto& kv : _ddbgc->getEvents() )
    {
        const Catalog::Event& event = kv.second;
        auto eqlrng = _ddbgc->getPhases().equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
            const Catalog::Phase& phase = it->second;
            auto xcorrCfg = _cfg.xcorr[phase.type];

            double duration = (xcorrCfg.endOffset - xcorrCfg.startOffset) + xcorrCfg.maxDelay * 2;
            Core::TimeSpan timeCorrection = Core::TimeSpan(xcorrCfg.startOffset) - Core::TimeSpan(xcorrCfg.maxDelay);
            Core::TimeWindow tw = Core::TimeWindow(phase.time + timeCorrection, duration);

            getWaveform(tw, event, phase, _wfCache, _useCatalogDiskCache);
        }
    }
}


/*
 * Build a catalog with requested phases only and for the same event/station pair
 * make sure to have only one phase. If multiple phases are found, keep the first
 * one arrived
 */
CatalogPtr HypoDD::filterOutPhases(const CatalogCPtr& catalog,
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
                        SEISCOMP_DEBUG("Preferring phase '%s' over '%s'",
                                       string(phase).c_str(), string(existingPhase).c_str());
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
                        SEISCOMP_DEBUG("Preferring phase '%s' over '%s'",
                                       string(phase).c_str(), string(existingPhase).c_str());
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
        phase.type = "P";
        filteredPhases.emplace(phase.eventId, phase);
    }
    for (auto& it : filteredS)
    {
        Catalog::Phase& phase = it.second;
        phase.type = "S";
        filteredPhases.emplace(phase.eventId, phase);
    }

    return new Catalog(catalog->getStations(), catalog->getEvents(), filteredPhases);
}


CatalogPtr HypoDD::relocateCatalog(bool force, bool usePh2dt)
{
    SEISCOMP_INFO("Starting HypoDD relocator in multiple events mode");

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
    _ddbgc->writeToFile((boost::filesystem::path(catalogWorkingDir)/"event.csv").string(),
                        (boost::filesystem::path(catalogWorkingDir)/"phase.csv").string(),
                        (boost::filesystem::path(catalogWorkingDir)/"station.csv").string() );

    // Create station.dat for hypodd (if not already generated)
    string stationFile = (boost::filesystem::path(catalogWorkingDir)/"station.dat").string();
    if ( force || ! Util::fileExists(stationFile) )
    {
        createStationDatFile(stationFile, _ddbgc);
    }

    string eventFile = (boost::filesystem::path(catalogWorkingDir)/"event.dat").string();
    string dtctFile = (boost::filesystem::path(catalogWorkingDir)/"dt.ct").string();
    string dtccFile = (boost::filesystem::path(catalogWorkingDir)/"dt.cc").string(); 

    if ( ! usePh2dt )
    {
        // Create event.dat for hypodd (if not already generated)
        if ( force || ! Util::fileExists(eventFile) )
        {
            createEventDatFile(eventFile, _ddbgc);
        }

        // calculate absolute travel times from catalog phases
        // Create dt.ct (if not already generated)
        if ( force || ! Util::fileExists(dtctFile) )
        {
            createDtCtCatalog(_ddbgc, dtctFile);
        }

        // calculate cross correlated differential travel times
        // Create dt.cc (if not already generated)
        if ( force || ! Util::fileExists(dtccFile) )
        {
            createDtCcCatalog(_ddbgc, dtccFile);
        }
    }
    else
    {
        // Create phase.dat for ph2dt (if not already generated)
        string phaseFile = (boost::filesystem::path(catalogWorkingDir)/"phase.dat").string();
        if ( force || ! Util::fileExists(phaseFile) )
        {
            createPhaseDatFile(phaseFile, _ddbgc);
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
            createDtCcPh2dt(dtctFile, dtccFile);
        }
    }

    // run hypodd
    // input : dt.cc dt.ct event.sel station.sel hypoDD.inp
    // output : hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
    string ddrelocFile = (boost::filesystem::path(catalogWorkingDir)/"hypoDD.reloc").string();
    string ddresidualFile = (boost::filesystem::path(catalogWorkingDir)/"hypoDD.res").string();
    if ( force || ! Util::fileExists(ddrelocFile) || ! Util::fileExists(ddresidualFile) )
    {
        runHypodd(catalogWorkingDir, dtccFile, dtctFile, eventFile, stationFile, _cfg.hypodd.step2CtrlFile);
    }

    // load a catalog from hypodd output file
    // input: hypoDD.reloc
    CatalogPtr relocatedCatalog = loadRelocatedCatalog(_ddbgc, ddrelocFile, ddresidualFile);

    // write catalog for debugging purpose
    relocatedCatalog->writeToFile((boost::filesystem::path(catalogWorkingDir)/"event-reloc.csv").string(),
                                  (boost::filesystem::path(catalogWorkingDir)/"phase-reloc.csv").string(),
                                  (boost::filesystem::path(catalogWorkingDir)/"station-reloc.csv").string());
    return relocatedCatalog;
}



CatalogPtr HypoDD::relocateSingleEvent(const CatalogCPtr& singleEvent)
{
    SEISCOMP_INFO("Starting HypoDD relocator in single event mode");

    const CatalogCPtr evToRelocateCat = filterOutPhases(singleEvent, _cfg.validPphases, _cfg.validSphases);

    // there must be only one event in the catalog, the origin to relocate
    const Catalog::Event& evToRelocate = evToRelocateCat->getEvents().begin()->second;

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

    string eventWorkingDir = (boost::filesystem::path(subFolder)/"step1").string();
    if ( !Util::createPath(eventWorkingDir) )
    {
        string msg = "Unable to create working directory: " + eventWorkingDir;
        throw runtime_error(msg);
    }

    // write catalog for debugging purpose
    _ddbgc->writeToFile((boost::filesystem::path(subFolder)/"event.csv").string(),
                        (boost::filesystem::path(subFolder)/"phase.csv").string(),
                        (boost::filesystem::path(subFolder)/"station.csv").string() );

    // Select neighbouring events
    CatalogPtr neighbourCat = selectNeighbouringEvents(
        _ddbgc->merge(evToRelocateCat), evToRelocate, _cfg.dtct.minWeight, _cfg.dtct.minESdist,
        _cfg.dtct.maxESdist, _cfg.dtct.minEStoIEratio, _cfg.dtct.maxIEdist, _cfg.dtct.minDTperEvt,
        _cfg.dtct.minNumNeigh, _cfg.dtct.maxNumNeigh
    );
    // add event
    neighbourCat->copyEvent(evToRelocate, evToRelocateCat, false);
    // extract the new id of the event
    unsigned evToRelocateNewId = neighbourCat->searchEvent(evToRelocate)->first;

    // Create station.dat for hypodd
    string stationFile = (boost::filesystem::path(eventWorkingDir)/"station.dat").string();
    createStationDatFile(stationFile, neighbourCat);

    // Create event.dat for hypodd
    string eventFile = (boost::filesystem::path(eventWorkingDir)/"event.dat").string();
    createEventDatFile(eventFile, neighbourCat);

    // Create differential travel times file (dt.ct) for hypodd
    string dtctFile = (boost::filesystem::path(eventWorkingDir)/"dt.ct").string();
    createDtCtSingleEvent(neighbourCat, evToRelocateNewId, dtctFile);
 
    // Create an empty cross correlated differential travel times file (dt.cc) for hypodd
    string dtccFile = (boost::filesystem::path(eventWorkingDir)/"dt.cc").string();
    ofstream(dtccFile).close();

    // run hypodd
    // input : dt.cc dt.ct event.sel station.sel hypoDD.inp
    // output : hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
    runHypodd(eventWorkingDir, dtccFile, dtctFile, eventFile, stationFile, _cfg.hypodd.step1CtrlFile);

    // Load the relocated origin from Hypodd
    string ddrelocFile = (boost::filesystem::path(eventWorkingDir)/"hypoDD.reloc").string();
    string ddresidualFile = (boost::filesystem::path(eventWorkingDir)/"hypoDD.res").string();
    CatalogPtr relocatedCatalog = loadRelocatedCatalog(neighbourCat, ddrelocFile, ddresidualFile);
    CatalogPtr relocatedEvCat = extractEvent(relocatedCatalog, evToRelocateNewId);
    // write catalog for debugging purpose
    relocatedCatalog->writeToFile((boost::filesystem::path(eventWorkingDir)/"event-reloc.csv").string(),
                                  (boost::filesystem::path(eventWorkingDir)/"phase-reloc.csv").string(),
                                  (boost::filesystem::path(eventWorkingDir)/"station-reloc.csv").string());

    if ( _workingDirCleanup )
    {
        boost::filesystem::remove_all(eventWorkingDir);
    }

    //
    // Step 2: relocate the refined location this time with cross correlation
    //

    CatalogPtr relocatedEvWithXcorr;

    try {

        eventWorkingDir = (boost::filesystem::path(subFolder)/"step2").string();
        if ( !Util::createPath(eventWorkingDir) )
        {
            string msg = "Unable to create working directory: " + eventWorkingDir;
            throw runtime_error(msg);
        }

        // Select neighbouring events from the relocated origin
        const Catalog::Event& relocatedEv = relocatedEvCat->getEvents().begin()->second;
        neighbourCat = selectNeighbouringEvents(
            _ddbgc->merge(relocatedEvCat), relocatedEv, _cfg.dtcc.minWeight, _cfg.dtcc.minESdist,
            _cfg.dtcc.maxESdist, _cfg.dtcc.minEStoIEratio, _cfg.dtcc.maxIEdist, _cfg.dtcc.minDTperEvt,
            _cfg.dtcc.minNumNeigh, _cfg.dtcc.maxNumNeigh
        );
        // add event
        neighbourCat->copyEvent(relocatedEv, relocatedEvCat, false);
        // extract the new id of the event
        unsigned refinedLocNewId = neighbourCat->searchEvent(relocatedEv)->first;

        // Create station.dat for hypodd
        stationFile = (boost::filesystem::path(eventWorkingDir)/"station.dat").string();
        createStationDatFile(stationFile, neighbourCat);

        // Create event.dat for hypodd
        eventFile = (boost::filesystem::path(eventWorkingDir)/"event.dat").string();
        createEventDatFile(eventFile, neighbourCat);

        // Create differential travel times file (dt.ct) for hypodd
        dtctFile = (boost::filesystem::path(eventWorkingDir)/"dt.ct").string();
        createDtCtSingleEvent(neighbourCat, refinedLocNewId, dtctFile);

        // Create cross correlated differential travel times file (dt.cc) for hypodd
        dtccFile = (boost::filesystem::path(eventWorkingDir)/"dt.cc").string();
        createDtCcSingleEvent(neighbourCat, refinedLocNewId, dtccFile);

        // run hypodd
        // input : dt.cc dt.ct event.sel station.sel hypoDD.inp
        // output : hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
        runHypodd(eventWorkingDir, dtccFile, dtctFile, eventFile, stationFile, _cfg.hypodd.step2CtrlFile);

        // Load the relocated origin from Hypodd
        ddrelocFile = (boost::filesystem::path(eventWorkingDir)/"hypoDD.reloc").string();
        ddresidualFile = (boost::filesystem::path(eventWorkingDir)/"hypoDD.res").string();
        relocatedCatalog = loadRelocatedCatalog(neighbourCat, ddrelocFile, ddresidualFile);
        relocatedEvWithXcorr = extractEvent(relocatedCatalog, refinedLocNewId);
        // write catalog for debugging purpose
        relocatedCatalog->writeToFile((boost::filesystem::path(eventWorkingDir)/"event-reloc.csv").string(),
                                      (boost::filesystem::path(eventWorkingDir)/"phase-reloc.csv").string(),
                                      (boost::filesystem::path(eventWorkingDir)/"station-reloc.csv").string());

        if ( _workingDirCleanup )
        {
            boost::filesystem::remove_all(eventWorkingDir);
        }

    } catch ( exception &e ) {
        SEISCOMP_ERROR("%s", e.what());
        SEISCOMP_ERROR("It was not possible to use cross correlation when relocating origin");
    }

    return relocatedEvWithXcorr ? relocatedEvWithXcorr : relocatedEvCat;
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
void HypoDD::createStationDatFile(const string& staFileName, const CatalogCPtr& catalog) const
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
void HypoDD::createPhaseDatFile(const string& phaseFileName, const CatalogCPtr& catalog) const
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
                SEISCOMP_ERROR("Ignoring phase '%s' with negative travel time (event '%s')",
                               string(phase).c_str(), string(event).c_str());
                continue; 
            }

            outStream << stringify("%-12s %12.6f %5.2f %4s",
                                  phase.stationId.c_str(), travel_time,
                                  phase.weight, phase.type.c_str());
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
void HypoDD::createEventDatFile(const string& eventFileName, const CatalogCPtr& catalog) const
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
    map<int,string> linesToReplace = { {1,stationFile}, {2,phaseFile} };
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
        {lineOffset + 1, dtccFile},
        {lineOffset + 2, dtctFile},
        {lineOffset + 3, eventFile},
        {lineOffset + 4, stationFile},
        {lineOffset + 5, "hypoDD.loc"},
        {lineOffset + 6, "hypoDD.reloc"},
        {lineOffset + 7, "hypoDD.sta"},
        {lineOffset + 8, "hypoDD.res"},
        {lineOffset + 9, "hypoDD.src"}
    };
    copyFileAndReplaceLines(ctrlFile, (boost::filesystem::path(workingDir)/"hypodd.inp").string(),
                            linesToReplace);

    // run Hypodd (use /bin/sh to get stdout/strerr redirection)
    string cmd = stringify("%s %s >hypodd.out 2>&1", _cfg.hypodd.exec.c_str(), "hypodd.inp");
    ::startExternalProcess({"/bin/sh", "-c", cmd}, true, workingDir);
}




CatalogPtr HypoDD::selectNeighbouringEvents(const CatalogCPtr& catalog,
                                            const Catalog::Event& refEv,
                                            double minPhaseWeight,
                                            double minESdist,
                                            double maxESdist,
                                            double minEStoIEratio,
                                            double maxIEdist,
                                            int minDTperEvt,
                                            int minNumNeigh,
                                            int maxNumNeigh) const
{
    SEISCOMP_DEBUG("Selecting Neighbouring Events for event %s", string(refEv).c_str());

    map<double,unsigned> eventByDistance; // distance, eventid
    map<unsigned,double> distanceByEvent; // eventid, distance
    map<unsigned,double> azimuthByEvent;  // eventid, azimuth

    // loop through every event in the catalog and select the ones within maxIEdist distance
    for (const auto& kv : catalog->getEvents() )
    {
        const Catalog::Event& event = kv.second;

        if (event == refEv)
            continue;

        // compute distance between current event and reference origin
        double azimuth;
        double distance = computeDistance(refEv.latitude, refEv.longitude, refEv.depth,
                                          event.latitude, event.longitude, event.depth,
                                          &azimuth);
        // too far away ?
        if ( maxIEdist > 0 && distance > maxIEdist )
            continue;

        // keep a list of added events sorted by distance
        eventByDistance[distance] = event.id;
        distanceByEvent[event.id] = distance;
        azimuthByEvent[event.id]  = azimuth;
    }

    // From the events within distance select the ones who respect the constraints
    list<unsigned> selectedEvents;
    set<string> includedStations;
    set<string> excludedStations;

    for (const auto& kv : eventByDistance)
    {
        const Catalog::Event& event = catalog->getEvents().at(kv.second);

        //  enought phases (> minDTperEvt) ?
        int dtCount = 0;
        auto eqlrng = catalog->getPhases().equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
            const Catalog::Phase& phase = it->second;

            // check pick weight
            if (phase.weight < minPhaseWeight)
                continue;

            // check events/station distance
            auto search = catalog->getStations().find(phase.stationId);
            if (search == catalog->getStations().end())
            {
                string msg = stringify("Malformed catalog: cannot find station '%s' "
                                   " referenced by phase '%s' for event '%u'",
                                   phase.stationId.c_str(), string(phase).c_str(),  event.id);
                throw runtime_error(msg);
            }

            const Catalog::Station& station = search->second;

            if ( excludedStations.find(station.id) != excludedStations.end() )
                continue;

            if ( includedStations.find(station.id) == includedStations.end() )
            {
                // compute distance between station and reference event
                double stationDistance = computeDistance(refEv.latitude, refEv.longitude, refEv.depth,
                                                         station.latitude, station.longitude,
                                                         -(station.elevation/1000.));

                // check this station distance is ok
                if ( ( maxESdist > 0 && stationDistance > maxESdist ) ||                  // too far away ?
                     ( stationDistance < minESdist )                 ||                  // too close ?
                     ( (stationDistance/distanceByEvent[event.id]) < minEStoIEratio ) ) // ratio too small ?
                {
                    excludedStations.insert(station.id);
                    continue;
                }
                includedStations.insert(station.id);
            }

            // compute distance between current event and station
            double stationDistance = computeDistance(event.latitude, event.longitude, event.depth,
                                                     station.latitude, station.longitude,
                                                     -(station.elevation/1000.));

            // check this station distance is ok
            if ( ( maxESdist > 0 && stationDistance > maxESdist ) ||                  // too far away ?
                 ( stationDistance < minESdist )                 ||                  // too close ?
                 ( (stationDistance / distanceByEvent[event.id]) < minEStoIEratio ) ) // ratio too small ?
            {
                continue;
            }

            // this phase is ok, now find corresponding phase in reference event phases
            auto eqlrng2 = catalog->getPhases().equal_range(refEv.id);
            for (auto it2 = eqlrng2.first; it2 != eqlrng2.second; ++it2)
            {
                const Catalog::Phase& refPhase = it2->second;
                if (phase.stationId == refPhase.stationId && 
                    phase.type == refPhase.type)
                {
                    if (refPhase.weight >= minPhaseWeight)
                        dtCount++;
                    break;
                }
            }
        }

        // not enought phases ?
        if ( dtCount < minDTperEvt )
        {
            continue;
        } 

        // add this event to the selected ones
        selectedEvents.push_back(event.id);
        SEISCOMP_DEBUG("Selecting possible event %u distance %.1f azimuth %.1f",
                       event.id, kv.first, azimuthByEvent[event.id]);
    }

    // Finally build the catalog of selected events
    // We prefer closer events over further ones but at the same time
    // we want to have a broad azimuth coverage
    CatalogPtr filteredCatalog = new Catalog();

    int numEvents = 0;
    int nextBin = 0;
    int numBins = 1;

    while ( !selectedEvents.empty() )
    {
        // select only maxNumNeigh
        if ( maxNumNeigh > 0 && numEvents >= maxNumNeigh )
            break;

        if ( ++nextBin == numBins )
        {
            numBins = std::min(8, int(selectedEvents.size()));
            if ( maxNumNeigh > 0 )
                numBins = std::min(numBins, maxNumNeigh-numEvents);
            nextBin = 0;
        }

        // Split numBins in 2 sets: one set for events above the reference event and
        // the other set below the reference event. This is to reduce depth errors
        int halfBins =  numBins / 2;
        int azimuthBins = nextBin < halfBins ? halfBins : numBins - halfBins;
        int expectedAzBin = nextBin < halfBins ? nextBin : nextBin - halfBins;
        bool expectedAbove = nextBin < halfBins;

        for (auto it = selectedEvents.begin(); it != selectedEvents.end(); it++)
        {
            const Catalog::Event& event = catalog->getEvents().at( *it );

            double azimuthBinSize = 360. / azimuthBins;
            int eventAzBin = int(azimuthByEvent[event.id] / azimuthBinSize) % azimuthBins;

            bool isAboveRefEv = event.depth < refEv.depth;
            bool depthOk = (isAboveRefEv == expectedAbove) || (numBins == 1);

            if ( eventAzBin == expectedAzBin && depthOk )
            {
                // add this event to the catalog
                filteredCatalog->copyEvent(event, catalog, true);
                numEvents++;
                selectedEvents.erase(it);
                SEISCOMP_DEBUG("Chose neighbour event %u distance %.1f azimuth %.1f depth %.3f "
                               "(azimuthBins %d eventAzBin %d deltaDepth %.3f)",
                               event.id, distanceByEvent[event.id], azimuthByEvent[event.id], event.depth,
                               azimuthBins, eventAzBin, (refEv.depth-event.depth));
                break;
            }
        }
    }

    // Check if enough neighbors were found
    if ( numEvents < minNumNeigh )
    {
        string msg = stringify("Skipping event %s, insufficient number of neighbors (%d)",
                               string(refEv).c_str(), numEvents);
        SEISCOMP_DEBUG("%s", msg.c_str());
        throw runtime_error(msg);
    }

    return filteredCatalog;
}



map<unsigned,CatalogPtr> 
HypoDD::selectNeighbouringEventsCatalog(const CatalogCPtr& catalog,
                                        double minPhaseWeight,
                                        double minESdist,
                                        double maxESdist,
                                        double minEStoIEratio,
                                        double maxIEdist,
                                        int minDTperEvt,
                                        int minNumNeigh,
                                        int maxNumNeigh) const
{
    SEISCOMP_INFO("Selecting Catalog Neighbouring Events ");

    // build the list of neighbours for each event
    map<unsigned,CatalogPtr> neighbourCats;
    CatalogPtr tmpCatalog = new Catalog(*catalog);
    bool redo;
    do {
        neighbourCats.clear();
        vector<unsigned> removedEvents;
        
        // for each event find the neighbours
        for (const auto& kv : tmpCatalog->getEvents() )
        {
            const Catalog::Event& event = kv.second;

            CatalogPtr neighbourCat; 
            try {
                neighbourCat = selectNeighbouringEvents(
                    tmpCatalog, event, minPhaseWeight, minESdist,
                    maxESdist, minEStoIEratio, maxIEdist,
                    minDTperEvt, minNumNeigh, maxNumNeigh
                );
            } catch ( ... ) { }

            if ( neighbourCat )
            {
                // add event to neighbour catalog list and update the list
                neighbourCat->copyEvent(event, catalog, true);
                neighbourCats[event.id] = neighbourCat;
            }
            else
            {
                removedEvents.push_back( event.id );
            }
        }

        // remove event because in next iterations we don't want this
        for (unsigned eventIdtoRemove : removedEvents)
        {
            tmpCatalog->removeEvent(eventIdtoRemove);
        }
        
        // check if the removed events were used as neighbor of any event
        // if so restart the loop
        redo = false;
        for (auto kv : neighbourCats )
        {
            CatalogPtr& currCat  = kv.second;
            for (unsigned eventIdtoRemove : removedEvents)
            {
                auto search = currCat->getEvents().find(eventIdtoRemove);
                if (search != currCat->getEvents().end())
                {
                    SEISCOMP_DEBUG("Restart Neighbouring Events selection ");
                    redo = true;
                    break;
                }
            }
            if ( redo ) break;
        }
        
    } while( redo );

    // We don't want to report the same pairs multiple times
    // when creating dt.cc and dt.ct (e.g. pair eventXX-eventYY
    // is the same as pair eventYY-eventXX), we'll remove the
    // pairs that appeared in previous catalogs from the
    // following catalogs
    std::multimap<unsigned,unsigned> existingPairs;

    for (auto kv : neighbourCats )
    {
        unsigned currEventId = kv.first;
        CatalogPtr& currCat  = kv.second;

        // remove from currrent catalogl the existing pairs
        auto eqlrng = existingPairs.equal_range(currEventId);
        for (auto existingPair = eqlrng.first; existingPair != eqlrng.second; existingPair++)
        {
            currCat->removeEvent(existingPair->second);
        }

        // remove current pairs from following catalogs
        for ( auto const& kv : currCat->getEvents() )
        {
            if ( kv.first != currEventId )
                existingPairs.emplace(kv.first, currEventId);
        }
    }

    return neighbourCats;
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
            string msg = stringify("Malformed catalog: cannot find relocated event %u"
                                   " in the original catalog", eventId);
            throw runtime_error(msg);
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
        event.rms = (event.relocInfo.rmsResidualCC + event.relocInfo.rmsResidualCT) / 2.;
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
            string key = to_string(phase.eventId) + "+" + phase.stationId + "+" + phase.type;
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


CatalogPtr HypoDD::extractEvent(const CatalogCPtr& catalog, unsigned eventId) const
{
    CatalogPtr eventToExtract = new Catalog();

    auto search = catalog->getEvents().find(eventId);
    if (search == catalog->getEvents().end())
    {
        string msg = stringify("Cannot find event id %u in the catalog.", eventId);
        throw runtime_error(msg);
    }

    const Catalog::Event& event = search->second;
    eventToExtract->addEvent(event, false);

    unsigned newEventId = eventToExtract->searchEvent(event)->first;

    auto eqlrng = catalog->getPhases().equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
        Catalog::Phase phase = it->second;

        auto search = catalog->getStations().find(phase.stationId);
        if (search == catalog->getStations().end())
        {
            string msg = stringify("Malformed catalog: cannot find station '%s' "
                                   " referenced by phase '%s' for event '%s'",
                                   phase.stationId.c_str(), string(phase).c_str(), string(event).c_str());
            throw runtime_error(msg);
        }
        const Catalog::Station& station = search->second;
        eventToExtract->addStation(station, true);
        string newStationId = eventToExtract->searchStation(station)->first;

        phase.eventId = newEventId;
        phase.stationId = newStationId;

        eventToExtract->addPhase(phase, true);
    } 

    return eventToExtract;
}



/* 
 * Create absolute travel times file (dt.ct) for hypodd
 * This is for full catalog mode
 */
void HypoDD::createDtCtCatalog(const CatalogCPtr& catalog,
                               const string& dtctFile) const
{
    SEISCOMP_INFO("Creating differential travel time file %s", dtctFile.c_str());

    map<unsigned,CatalogPtr> neighbourCats = selectNeighbouringEventsCatalog(
        catalog, _cfg.dtct.minWeight, _cfg.dtct.minESdist,
        _cfg.dtct.maxESdist, _cfg.dtct.minEStoIEratio, _cfg.dtct.maxIEdist,
        _cfg.dtct.minDTperEvt, _cfg.dtct.minNumNeigh, _cfg.dtct.maxNumNeigh
    );

    ofstream outStream(dtctFile);
    if ( !outStream.is_open() )
        throw runtime_error("Cannot create file " + dtctFile);

    for (const auto& kv : neighbourCats)
        buildAbsTTimePairs(kv.second, kv.first, outStream);
}



/* 
 * Create absolute travel times file (dt.ct) for hypodd
 * This is for single event mode
 */
void HypoDD::createDtCtSingleEvent(const CatalogCPtr& catalog,
                                   unsigned evToRelocateId,
                                   const string& dtctFile) const
{
    SEISCOMP_INFO("Creating differential travel time file %s", dtctFile.c_str());

    ofstream outStream(dtctFile);
    if ( !outStream.is_open() )
        throw runtime_error("Cannot create file " + dtctFile);

    buildAbsTTimePairs(catalog, evToRelocateId, outStream);
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
void HypoDD::buildAbsTTimePairs(const CatalogCPtr& catalog,
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

            // loop through reference event phases
            auto eqlrng2 = catalog->getPhases().equal_range(refEv.id);
            for (auto it2 = eqlrng2.first; it2 != eqlrng2.second; ++it2)
            {
                const Catalog::Phase& refPhase = it2->second;

                if (phase.stationId == refPhase.stationId && 
                    phase.type == refPhase.type)
                {
                    double ref_travel_time = refPhase.time - refEv.time;
                    if (ref_travel_time < 0)
                    {
                        SEISCOMP_WARNING("Ignoring phase '%s' with negative travel time (event '%s')",
                                       string(refPhase).c_str(), string(refEv).c_str());
                        break;
                    }
                    double travel_time = phase.time - event.time;
                    if (travel_time < 0)
                    {
                        SEISCOMP_WARNING("Ignoring phase '%s' with negative travel time (event '%s')",
                                       string(phase).c_str(), string(event).c_str());
                        break;
                    }

                    // get common observation weight for pair (FIXME: take the lower one? average?)
                    double weight = (refPhase.weight + phase.weight) / 2.0;

                    evStream << stringify("%-12s %.6f %.6f %.2f %s",
                                          refPhase.stationId.c_str(), ref_travel_time,
                                          travel_time, weight, refPhase.type.c_str());
                    evStream << endl;
                    dtCount++;
                    break;
                }
            }
        }
        if (dtCount > 0)
            outStream << evStream.str();
    }
}



/*
 * Compute and store to file differential travel times from cross
 * correlation for pairs of earthquakes.
 * This is for full catalog mode
 */
void HypoDD::createDtCcCatalog(const CatalogCPtr& catalog,
                               const string& dtccFile)
{
    SEISCOMP_INFO("Creating Cross correlation differential travel time file %s", dtccFile.c_str());

    map<unsigned,CatalogPtr> neighbourCats = selectNeighbouringEventsCatalog(
        catalog, _cfg.dtcc.minWeight, _cfg.dtcc.minESdist,
        _cfg.dtcc.maxESdist, _cfg.dtcc.minEStoIEratio, _cfg.dtcc.maxIEdist,
        _cfg.dtcc.minDTperEvt, _cfg.dtcc.minNumNeigh, _cfg.dtcc.maxNumNeigh
    );

    ofstream outStream(dtccFile);
    if ( !outStream.is_open() )
        throw runtime_error("Cannot create file " + dtccFile);

    for (const auto& kv : neighbourCats)
    {
        buildXcorrDiffTTimePairs(kv.second, kv.first, outStream,
                                 _wfCache, _useCatalogDiskCache,
                                 _wfCache, _useCatalogDiskCache);
    }
}


/*
 * Compute and store to file differential travel times from cross
 * correlation for pairs of earthquakes.
 * This is for single event mode
 */
void HypoDD::createDtCcSingleEvent(const CatalogCPtr& catalog,
                                   unsigned evToRelocateId,
                                   const string& dtccFile)
{
    SEISCOMP_INFO("Creating Cross correlation differential travel time file %s", dtccFile.c_str());

    ofstream outStream(dtccFile);
    if ( !outStream.is_open() )
        throw runtime_error("Cannot create file " + dtccFile);

    map<string, GenericRecordPtr> tmpCahe;

    buildXcorrDiffTTimePairs(catalog, evToRelocateId, outStream,
                             _wfCache, _useCatalogDiskCache,
                             (_useSingleEvDiskCache ? _wfCache : tmpCahe),
                              _useSingleEvDiskCache);
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
void HypoDD::buildXcorrDiffTTimePairs(const CatalogCPtr& catalog,
                                      unsigned evToRelocateId,
                                      ofstream& outStream,
                                      std::map<std::string,GenericRecordPtr>& catalogCache,
                                      bool useDiskCacheCatalog,
                                      std::map<std::string,GenericRecordPtr>& refEvCache,
                                      bool useDiskCacheRefEv)
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
        evStream << stringify("# %10u %10u       0.0", refEv.id, event.id) << endl;

        // loop through event phases
        auto eqlrng = catalog->getPhases().equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
            const Catalog::Phase& phase = it->second;

            // loop through reference event phases
            auto eqlrng2 = catalog->getPhases().equal_range(refEv.id);
            for (auto it2 = eqlrng2.first; it2 != eqlrng2.second; ++it2)
            {
                const Catalog::Phase& refPhase = it2->second;

                if (phase.stationId == refPhase.stationId && 
                    phase.type == refPhase.type)
                {
                    double dtcc, weight;
                    if ( xcorr(refEv, refPhase, event, phase, dtcc, weight,
                               refEvCache, useDiskCacheRefEv,
                               catalogCache, useDiskCacheCatalog) )
                    {
                        evStream << stringify("%-12s %.6f %.4f %s", refPhase.stationId.c_str(),
                                              dtcc, weight,  refPhase.type.c_str());
                        evStream << endl;
                        dtCount++;
                    }
                    break;
                }
            }
        }
        if (dtCount > 0)
            outStream << evStream.str();
    }
}



/*
 * Reads the event pairs matched in dt.ct which are selected by ph2dt and
 * calculate cross correlated differential travel_times for every pair.
 * input dt.ct 
 * output dt.cc
 */
void HypoDD::createDtCcPh2dt(const string& dtctFile, const string& dtccFile)
{
    SEISCOMP_INFO("Creating Cross correlation differential travel time file %s", dtccFile.c_str());

    if ( !Util::fileExists(dtctFile) )
        throw runtime_error("Unable to perform cross correlation, cannot find file: " + dtctFile);

    ofstream outStream(dtccFile);
    if ( !outStream.is_open() )
        throw runtime_error("Cannot create file " + dtccFile);

    const std::map<unsigned,Catalog::Event>& events = _ddbgc->getEvents();
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
                string msg = stringify("Relocated catalog contains events ids (%s or %s) "
                                       "that are not present in the original catalog.",
                                       string(*ev1).c_str(), string(*ev2).c_str());
                throw runtime_error(msg.c_str());
            }
            ev1 = &search1->second;
            ev2 = &search2->second;

            // write the pairs has been built up to now
            if (dtCount > 0 )
                outStream << evStream.str();
            evStream = stringstream();
            dtCount = 0;

            evStream << stringify("# %10u %10u       0.0", ev1->id, ev2->id) << endl;
        }
        // observation line (STA, TT1, TT2, WGHT, PHA)
        else if(ev1 != nullptr && ev2 != nullptr && fields.size() == 5)
        {
            string stationId = fields[0];
            string phaseType = fields[4];

            // loop through event 1 phases
            auto eqlrng = _ddbgc->getPhases().equal_range(ev1->id);
            for (auto it = eqlrng.first; it != eqlrng.second; ++it)
            {
                const Catalog::Phase& phase1 = it->second;
                if (phase1.stationId == stationId &&
                    phase1.type == phaseType )
                {
                    // loop through event 2 phases
                    eqlrng = _ddbgc->getPhases().equal_range(ev2->id);
                    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
                    {
                        const Catalog::Phase& phase2 = it->second;
                        if (phase2.stationId == stationId &&
                            phase2.type == phaseType )
                        {
                            double dtcc, weight;
                            if ( xcorr( *ev1, phase1, *ev2, phase2, dtcc, weight,
                                       _wfCache, _useCatalogDiskCache,
                                       _wfCache, _useCatalogDiskCache) )
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
}



string 
HypoDD::waveformId(const Catalog::Phase& ph, const Core::TimeWindow& tw) const
{
    return waveformId(ph.networkCode, ph.stationCode, ph.locationCode, ph.channelCode, tw);
}

string
HypoDD::waveformId(const string& networkCode, const string& stationCode,
                   const string& locationCode, const string& channelCode,
                   const Core::TimeWindow& tw) const
{
    return stringify("%s.%s.%s.%s.%s.%s",
                     networkCode.c_str(), stationCode.c_str(),
                     locationCode.c_str(), channelCode.c_str(),
                     tw.startTime().iso().c_str(),
                     tw.endTime().iso().c_str());
}


bool
HypoDD::xcorr(const Catalog::Event& event1, const Catalog::Phase& phase1,
              const Catalog::Event& event2, const Catalog::Phase& phase2,
              double& dtccOut, double& weightOut,
              std::map<std::string,GenericRecordPtr>& cache1,  bool useDiskCache1,
              std::map<std::string,GenericRecordPtr>& cache2,  bool useDiskCache2)
{
    dtccOut = 0;
    weightOut = 0;
    auto xcorrCfg = _cfg.xcorr[phase1.type];

    SEISCOMP_DEBUG("Calculating cross correlation for phase pair phase1='%s', phase2='%s'",
                   string(phase1).c_str(), string(phase2).c_str());

    // compute start time and duration for the two traces
    double shortDuration = xcorrCfg.endOffset - xcorrCfg.startOffset;
    Core::TimeSpan shortTimeCorrection = Core::TimeSpan(xcorrCfg.startOffset);
    double longDuration = shortDuration + xcorrCfg.maxDelay * 2;
    Core::TimeSpan longTimeCorrection = shortTimeCorrection - Core::TimeSpan(xcorrCfg.maxDelay);

    Core::TimeWindow tw1 = Core::TimeWindow(phase1.time + longTimeCorrection, longDuration);
    Core::TimeWindow tw2 = Core::TimeWindow(phase2.time + longTimeCorrection, longDuration); 

    // Check if we have already excluded those traces (couldn't load, snr too high, etc..)
    const string wfId1 = waveformId(phase1, tw1);
    const string wfId2 = waveformId(phase2, tw2);
    if ( _excludedWfs.count(wfId1) != 0 || _excludedWfs.count(wfId2) != 0 )
    {
        return false;
    }

    // load the long trace 1, because we want to cache the long version. Then we'll trim it.
    GenericRecordPtr tr1 = getWaveform(tw1, event1, phase1, cache1, useDiskCache1);
    if ( !tr1 )
    {
        _excludedWfs.insert(wfId1);
        return false;
    }

    // load the long trace 2, because we want to cache the long version. Then we'll trim it
    GenericRecordPtr tr2 = getWaveform(tw2, event2, phase2, cache2, useDiskCache2);
    if ( !tr2 )
    {
        _excludedWfs.insert(wfId2);
        return false;
    }

    // trim tr2 to shorter length, we want to cross correlate the short with the long one
    GenericRecordPtr tr2Short = new GenericRecord(*tr2);
    Core::TimeWindow tw2Short = Core::TimeWindow(phase2.time + shortTimeCorrection, shortDuration);
    if ( !trim(*tr2Short, tw2Short) )
    {
        SEISCOMP_WARNING("Cannot trim phase2 waveform, skipping cross correlation "
                         "for phase pair phase1='%s', phase2='%s'",
                         string(phase1).c_str(), string(phase2).c_str());
        _excludedWfs.insert(wfId2);
        return false;
    }

    double xcorr_coeff, xcorr_dt;
    if ( ! xcorr(tr1, tr2Short, xcorrCfg.maxDelay, true, xcorr_dt, xcorr_coeff) )
    {
        return false;
    }

    // trim tr1 to shorter length, we want to cross correlate the short with the long one
    GenericRecordPtr tr1Short = new GenericRecord(*tr1);
    Core::TimeWindow tw1Short = Core::TimeWindow(phase1.time + shortTimeCorrection, shortDuration);
    if ( !trim(*tr1Short, tw1Short) )
    {
        SEISCOMP_WARNING("Cannot trim phase1 waveform, skipping cross correlation "
                         "for phase pair phase1='%s', phase2='%s'",
                         string(phase1).c_str(), string(phase2).c_str());
        _excludedWfs.insert(wfId1);
        return false;
    }

    double xcorr_coeff2, xcorr_dt2;
    if ( ! xcorr(tr1Short, tr2, xcorrCfg.maxDelay, true, xcorr_dt2, xcorr_coeff2) )
    {
        return false;
    }

    if ( ! std::isfinite(xcorr_coeff) && ! std::isfinite(xcorr_coeff2) )
    {
        return false;
    }

    if ( ! std::isfinite(xcorr_coeff) ||
         (std::isfinite(xcorr_coeff2) && xcorr_coeff2 > xcorr_coeff) )
    {
        xcorr_coeff = xcorr_coeff2;
        xcorr_dt = xcorr_dt2;
    }

    if ( xcorr_coeff < xcorrCfg.minCoef )
    {
        return false;
    }

    // compute differential travel time and weight of measurement
    double travel_time1 = phase1.time - event1.time;
    double travel_time2 = phase2.time - event2.time;
    dtccOut = travel_time1 - travel_time2 - xcorr_dt;
    weightOut = xcorr_coeff * xcorr_coeff;

    return true;
}



// Calculate the correlation series (tr1 and tr2 are already demeaned)
bool
HypoDD::xcorr(const GenericRecordCPtr& tr1, const GenericRecordCPtr& tr2, double maxDelay,
              double qualityCheck, double& delayOut, double& coeffOut) const
{
    delayOut = 0.;
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
    GenericRecordCPtr trLonger = swap ? tr1 : tr2; 

    const double *smpsS = DoubleArray::ConstCast(trShorter->data())->typedData();
    const double *smpsL = DoubleArray::ConstCast(trLonger->data())->typedData();
    const int smpsSsize = trShorter->data()->size();
    const int smpsLsize = trLonger->data()->size();

    // for later quality check: save local maxima
    vector<double> localMaxs;
    bool notDecreasing = false;
    double prevCoeff = -1;

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
        if ( coeff > coeffOut || !std::isfinite(coeffOut) )
        {
            coeffOut = coeff;
            delayOut = delay / freq; // samples to secs
        }

        // for later quality check: save local maxima
        if (coeff < prevCoeff && notDecreasing ) localMaxs.push_back(prevCoeff);
        notDecreasing = coeff >= prevCoeff;
        prevCoeff = coeff;
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
    if ( qualityCheck && std::isfinite(coeffOut) )
    {
        double threshold = coeffOut - ( (1.0 - coeffOut) / 2.0 );
        int numMax = 0;
        for (double CCslm : localMaxs)
        {
            if (std::isfinite(CCslm) && CCslm >= threshold) numMax++;
            if (numMax > 1)
            {
                coeffOut = std::nan("");
                SEISCOMP_DEBUG("Cycle skipping detected when cross correlating traces");
                break;
            }
        }
    }

    return true;
}



double
HypoDD::S2Nratio(const GenericRecordCPtr& tr, const Core::Time& pickTime,
                 double noiseOffsetStart, double noiseOffsetEnd,
                 double signalOffsetStart, double signalOffsetEnd) const
{
    const double *data = DoubleArray::ConstCast(tr->data())->typedData();
    const int size = tr->data()->size();
    const double freq = tr->samplingFrequency();
    const Core::Time dataStartTime = tr->startTime();

    // convert time w.r.t. guiding pick time to sample number
    auto secToSample = [&, freq, size](double sec) { 
        return std::min(std::max(std::round(sec * freq), 0.), size-1.);
    };
    const double pickOffset = (pickTime - dataStartTime).length();
    const int noiseStart  = secToSample(noiseOffsetStart  + pickOffset);
    const int noiseEnd    = secToSample(noiseOffsetEnd    + pickOffset);
    const int signalStart = secToSample(signalOffsetStart + pickOffset);
    const int signalEnd   = secToSample(signalOffsetEnd   + pickOffset);

    if ( (std::min({noiseStart,noiseEnd,signalStart,signalEnd}) < 0)    ||
         (std::max({noiseStart,noiseEnd,signalStart,signalEnd}) >= size) )
    {
        SEISCOMP_ERROR("Cannot compute S2N ratio: noise/signal windows exceed waveform boundaries");
        return -1;
    }

    // Get maximum (absolute) amplitude in noise window:
    double noiseMax = -1.0;
    for (int i = noiseStart; i < noiseEnd; i++)
    {
        noiseMax = std::max(std::abs(data[i]), noiseMax);
    }

    // Get maximum (absolute) amplitude in signal window:
    double signalMax = -1.0;
    for (int i = signalStart; i < signalEnd; i++)
    {
        signalMax = std::max(std::abs(data[i]), signalMax);
    }

    return signalMax/noiseMax;
}



/*
 * Return the waveform from the memory cache if present, otherwise load it
 */
GenericRecordPtr
HypoDD::getWaveform(const Core::TimeWindow& tw,
                    const Catalog::Event& ev,
                    const Catalog::Phase& ph,
                    map<string,GenericRecordPtr>& memCache,
                    bool useDiskCache) const
{
    string wfDesc = stringify("Waveform for Phase '%s' and Time slice from %s length %.2f sec",
                              string(ph).c_str(), tw.startTime().iso().c_str(), tw.length());

    /*
     * first try to load the waveform from the memory cache, if present
     */
    const string wfId = waveformId(ph, tw);

    const auto it = memCache.find(wfId);

    if ( it != memCache.end() )
    {
        return it->second; // waveform cached, just return it 
    }

    /*
     * Load the waveform, possibly perform a projection 123->ZNE or ZNE->ZRT,
     * filter it and finally save the result in the memory cache  for later re-use
     */
    bool projectionRequired = true;

    string channelCodeRoot = ph.channelCode.substr(0, ph.channelCode.size()-1);
    Core::Time refTime = tw.startTime();

    DataModel::ThreeComponents tc;
    DataModel::Station* station = findStation(ph.networkCode, ph.stationCode, refTime);
    DataModel::SensorLocation *loc = findSensorLocation(station, ph.locationCode, refTime);

    bool allComponents = false;

    if ( loc )
    {
        allComponents = getThreeComponents(tc, loc, channelCodeRoot.c_str(), refTime);

        if ( ( tc.comps[ThreeComponents::Vertical] &&
               tc.comps[ThreeComponents::Vertical]->code() == ph.channelCode)        ||
             ( tc.comps[ThreeComponents::FirstHorizontal] &&
               tc.comps[ThreeComponents::FirstHorizontal]->code() == ph.channelCode) ||
             ( tc.comps[ThreeComponents::SecondHorizontal] &&
               tc.comps[ThreeComponents::SecondHorizontal]->code() == ph.channelCode)
           )
        {
            projectionRequired = false;
        }
    }
    else
    {
        SEISCOMP_WARNING("Unable to fetch SensorLocation info (%s)", wfDesc.c_str());
        projectionRequired = false; // let's try to load the waveform anyway
    }

    // Compute the waveform time window to load
    Core::TimeWindow twToLoad = tw;
    if ( _cfg.snr.minSnr > 0 )
    {
        // if the SNR window is bigger than the xcorr window, than extend
        // the waveform time window
        Core::Time winStart = std::min( {
              tw.startTime(),
              ph.time + Core::TimeSpan(_cfg.snr.noiseStart),
              ph.time + Core::TimeSpan(_cfg.snr.signalStart)
        });
        Core::Time winEnd = std::max( {
                tw.endTime(),
                ph.time + Core::TimeSpan(_cfg.snr.noiseEnd),
                ph.time + Core::TimeSpan(_cfg.snr.signalEnd)
        });
        twToLoad = Core::TimeWindow(winStart, winEnd);
    }

    // Load waveform:
    // - if no projection required, just load the requested component
    // - otherwise perform the projection 123->ZNE or ZNE->ZRT
    GenericRecordPtr trace;
    try {

        if ( ! projectionRequired )
        {
            trace = loadWaveform(twToLoad, ph.networkCode, ph.stationCode,
                                 ph.locationCode, ph.channelCode, useDiskCache);
        }
        else 
        {
            if ( ! allComponents )
            {
                SEISCOMP_WARNING("Unable to fetch orientation information (%s)", wfDesc.c_str());
                return nullptr;
            }
            trace = loadProjectWaveform(twToLoad, ev, ph, tc, loc, useDiskCache);
        }

    } catch ( exception &e ) {
        SEISCOMP_WARNING("%s", e.what());
        return nullptr;
    }

    // fitler waveform
    filter(*trace, true, _cfg.wfFilter.filterOrder, _cfg.wfFilter.filterFmin,
           _cfg.wfFilter.filterFmax, _cfg.wfFilter.resampleFreq);

    // check SNR threshold
    if ( _cfg.snr.minSnr > 0 )
    {
        double snr = S2Nratio(trace, ph.time, _cfg.snr.noiseStart, _cfg.snr.noiseEnd,
                              _cfg.snr.signalStart, _cfg.snr.signalEnd);
        if ( snr < _cfg.snr.minSnr ) 
        {
            SEISCOMP_DEBUG("Trace has too low SNR (%.2f), discard it (%s)", snr, wfDesc.c_str());
            return nullptr;
        }
    }

    // Trim waveform in case we loaded more data than requested (to compute SNR)
    if ( twToLoad != tw )
    {
        if ( ! trim(*trace, tw) )
        {
            string msg = stringify("Incomplete trace, not enough data (%s)", wfDesc.c_str());
            throw runtime_error(msg);
        }
    }

    // save waveform into the cache
    memCache[wfId] = trace;

    return trace; 
}



GenericRecordPtr
HypoDD::loadProjectWaveform(const Core::TimeWindow& tw,
                            const Catalog::Event& ev,
                            const Catalog::Phase& ph,
                            const DataModel::ThreeComponents& tc,
                            const DataModel::SensorLocation *loc,
                            bool useDiskCache) const
{
    string wfDesc = stringify("Waveform for Phase '%s' and Time slice from %s length %.2f sec",
                              string(ph).c_str(), tw.startTime().iso().c_str(), tw.length());

    SEISCOMP_DEBUG("Loading the 3 components waveforms (%s %s %s) to perform the projection...",
                   tc.comps[ThreeComponents::Vertical]->code().c_str(),
                   tc.comps[ThreeComponents::FirstHorizontal]->code().c_str(),
                   tc.comps[ThreeComponents::SecondHorizontal]->code().c_str());

    // orientation ZNE
    Math::Matrix3d orientationZNE;
    Math::Vector3d n;
    n.fromAngles(+deg2rad(tc.comps[ThreeComponents::Vertical]->azimuth()),
                 -deg2rad(tc.comps[ThreeComponents::Vertical]->dip())).normalize();
    orientationZNE.setColumn(2, n);
    n.fromAngles(+deg2rad(tc.comps[ThreeComponents::FirstHorizontal]->azimuth()),
                 -deg2rad(tc.comps[ThreeComponents::FirstHorizontal]->dip())).normalize();
    orientationZNE.setColumn(1, n);
    n.fromAngles(+deg2rad(tc.comps[ThreeComponents::SecondHorizontal]->azimuth()),
                 -deg2rad(tc.comps[ThreeComponents::SecondHorizontal]->dip())).normalize();
    orientationZNE.setColumn(0, n);

    // orientation ZRT
    Math::Matrix3d orientationZRT;
    double delta, az, baz;
    Math::Geo::delazi(ev.latitude, ev.longitude, loc->latitude(), loc->longitude(),
                      &delta, &az, &baz);
    orientationZRT.loadRotateZ(deg2rad(baz + 180.0));

    // transformation matrix
    Math::Matrix3d transformation;
    map<string,string> chCodeMap;

    string channelCodeRoot = ph.channelCode.substr(0, ph.channelCode.size()-1);
    char component = *ph.channelCode.rbegin();

    if ( component == 'Z' || component == 'N'  || component == 'E' )
    {
        transformation = orientationZNE;
        chCodeMap[channelCodeRoot + "Z"] = tc.comps[ThreeComponents::Vertical]->code();
        chCodeMap[channelCodeRoot + "N"] = tc.comps[ThreeComponents::FirstHorizontal]->code();
        chCodeMap[channelCodeRoot + "E"] = tc.comps[ThreeComponents::SecondHorizontal]->code();

        SEISCOMP_DEBUG("Performing ZNE projection (channelCode %s -> %s) for %s",
            chCodeMap[ph.channelCode].c_str(), ph.channelCode.c_str(), wfDesc.c_str());
    }
    else if ( component == 'R'  || component == 'T' )
    {
        transformation.mult(orientationZRT, orientationZNE);
        //chCodeMap[channelCodeRoot + "Z"] = tc.comps[ThreeComponents::Vertical]->code();
        chCodeMap[channelCodeRoot + "R"] = tc.comps[ThreeComponents::FirstHorizontal]->code();
        chCodeMap[channelCodeRoot + "T"] = tc.comps[ThreeComponents::SecondHorizontal]->code();

        SEISCOMP_DEBUG("Performing ZRT projection (channelCode %s -> %s) for %s",
            chCodeMap[ph.channelCode].c_str(), ph.channelCode.c_str(), wfDesc.c_str()); 
    }
    else
    {
        string msg = stringify("Unknown channel '%c', cannot load %s", component, wfDesc.c_str());
        throw runtime_error(msg); 
    }

    // Load the available components
    GenericRecordPtr tr1 = loadWaveform(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                                        tc.comps[ThreeComponents::Vertical]->code(), useDiskCache);
    GenericRecordPtr tr2 = loadWaveform(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                                        tc.comps[ThreeComponents::FirstHorizontal]->code(), useDiskCache);
    GenericRecordPtr tr3 = loadWaveform(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                                        tc.comps[ThreeComponents::SecondHorizontal]->code(), useDiskCache);

    // The wrapper will direct 3 codes into the right slots using the
    // Stream configuration class and will finally use the transformation
    // operator. The advantage is that it will apply the configured gain for
    typedef Operator::StreamConfigWrapper<double,3,Operator::Transformation> OpWrapper;

    // Define the final operator class:
    //  1. Send channel codes to right slots
    //  2. Align 3 channels sample wise
    //  3. Transform the resulting 3 component trace with a rotation matrix
    typedef NCompsOperator<double,3,OpWrapper> Rotator;

    Processing::Stream streams[3];
    streams[0].init(tc.comps[0]);
    streams[1].init(tc.comps[1]);
    streams[2].init(tc.comps[2]);
    Rotator op(OpWrapper(streams, Operator::Transformation<double,3>(transformation)));


    class DataStorer
    {
        public:
        DataStorer(const string& channelCode, const Core::TimeWindow& tw, map<string,string> chCodeMap)
        : chMap(chCodeMap[channelCode], channelCode)
        , _seq( new TimeWindowBuffer(tw) )
        {  }

        bool store(const Record *rec)
        {
            if (rec->channelCode() == chMap.first)
                _seq->feed(rec);
            return true;
        }

        const std::pair<string,string> chMap;
        std::shared_ptr<RecordSequence> _seq;
    };

    DataStorer projectedData(ph.channelCode, tw, chCodeMap);

    // The function that will be called after a transformed record was created
    //op.setStoreFunc(boost::bind(&RecordSequence::feed, seq, _1));
    op.setStoreFunc(boost::bind(&DataStorer::store, projectedData, _1));

    op.feed(tr1.get());
    op.feed(tr2.get());
    op.feed(tr3.get());

    std::shared_ptr<RecordSequence> seq = projectedData._seq;

    if ( seq->empty() )
    {
        string msg = stringify("No data after the projection for %s", wfDesc.c_str());
        throw runtime_error(msg);
    }

    GenericRecordPtr trace = new GenericRecord();

    if ( ! merge(*trace, *seq) )
    {
        string msg = stringify("Data records could not be merged into a single trace (%s)", wfDesc.c_str());
        throw runtime_error(msg);
    }

    trace->setChannelCode(ph.channelCode);

    if ( ! trim(*trace, tw) )
    {
        string msg = stringify("Incomplete trace, not enough data (%s)", wfDesc.c_str());
        throw runtime_error(msg);
    }

    return trace;
}



/*
 * Read a waveform from a chached copy on disk if present, otherwise
 * from the configured RecordStream
 */
GenericRecordPtr
HypoDD::loadWaveform(const Core::TimeWindow& tw,
                     const string& networkCode,
                     const string& stationCode,
                     const string& locationCode,
                     const string& channelCode,
                     bool useDiskCache) const
{
    string cacheFile = waveformId(networkCode, stationCode, locationCode, channelCode, tw) + ".mseed";

    cacheFile = (boost::filesystem::path(_cacheDir)/cacheFile).string();

    GenericRecordPtr trace;
    // First try to read trace from disk cache
    if ( useDiskCache && Util::fileExists(cacheFile) )
    {
        try {
            trace = readTrace(cacheFile);
        } catch ( ... ) {
            SEISCOMP_WARNING("Couldn't load cached waveform %s, read it from record stream",
                             cacheFile.c_str());
        }
    }

    // if the trace is not cached then read it from the configured recordStream
    if ( !trace )
    {
        trace = readWaveformFromRecordStream(tw,networkCode, stationCode, locationCode, channelCode);
        // then save the trace to disk for later usage
        if ( useDiskCache )
        {
            try {
                writeTrace(trace, cacheFile);
            } catch ( ... ) {
                SEISCOMP_WARNING("Couldn't write waveform cache to disk %s", cacheFile.c_str());
            }
        }
    }

    return trace;
}


GenericRecordPtr
HypoDD::readWaveformFromRecordStream(const Core::TimeWindow& tw,
                                     const string& networkCode,
                                     const string& stationCode,
                                     const string& locationCode,
                                     const string& channelCode) const
{
    IO::RecordStreamPtr rs = IO::RecordStream::Open( _cfg.dtcc.recordStreamURL.c_str() );
    if ( rs == nullptr )
    {
        string msg = "Cannot open RecordStream: " + _cfg.dtcc.recordStreamURL;
        throw runtime_error(msg);
    }

    rs->setTimeWindow(tw);
    rs->addStream(networkCode, stationCode, locationCode, channelCode);

    // Store each record in a RecordSequence
    IO::RecordInput inp(rs.get(), Array::DOUBLE, Record::DATA_ONLY);
    std::shared_ptr<RecordSequence> seq( new TimeWindowBuffer(tw) );
    RecordPtr rec;
    while ( rec = inp.next() )
    {
        seq->feed(rec.get());
    }

    if ( seq->empty() )
    {
        string msg = stringify("Data could not be loaded (stream %s.%s.%s.%s from %s length %.2f)",
                               networkCode.c_str(), stationCode.c_str(),
                               locationCode.c_str(), channelCode.c_str(),
                               tw.startTime().iso().c_str(), tw.length());
        throw runtime_error(msg);
    }

    GenericRecordPtr trace = new GenericRecord();

    if ( ! merge(*trace, *seq) )
    {
        string msg = stringify("Data records could not be merged into a single trace "
                               "(%s.%s.%s.%s from %s length %.2f)",
                               networkCode.c_str(), stationCode.c_str(),
                               locationCode.c_str(), channelCode.c_str(),
                               tw.startTime().iso().c_str(), tw.length());
        throw runtime_error(msg);
    }

    if ( ! trim(*trace, tw) )
    {
        string msg = stringify("Incomplete trace, not enough data for requested"
                               " time window (%s.%s.%s.%s from %s length %.2f)",
                               networkCode.c_str(), stationCode.c_str(),
                               locationCode.c_str(), channelCode.c_str(),
                               tw.startTime().iso().c_str(), tw.length());
        throw runtime_error(msg);
    }

    return trace;
}


bool HypoDD::merge(GenericRecord &trace, const RecordSequence& seq) const
{
    if ( seq.empty() )
    {
        return false;
    }

    RecordCPtr first = seq.front();
    RecordCPtr last;
    double samplingFrequency = first->samplingFrequency();
    Core::TimeSpan maxAllowedGap, maxAllowedOverlap;

    maxAllowedGap = Core::TimeSpan((double)(0.5 / samplingFrequency));
    maxAllowedOverlap = Core::TimeSpan((double)(-0.5 / samplingFrequency));

    trace.setNetworkCode(first->networkCode());
    trace.setStationCode(first->stationCode());
    trace.setLocationCode(first->locationCode());
    trace.setChannelCode(first->channelCode());

    trace.setStartTime(first->startTime());
    trace.setSamplingFrequency(samplingFrequency);

    Array::DataType datatype = first->data()->dataType();
    ArrayPtr arr = ArrayFactory::Create(datatype, datatype, 0, nullptr);

    for (const RecordCPtr& rec : seq )
    {
        if ( rec->samplingFrequency() != samplingFrequency ) {
            SEISCOMP_DEBUG("%s.%s.%s.%s: record sampling frequencies are not consistent: %f != %f",
                          trace.networkCode().c_str(),
                          trace.stationCode().c_str(),
                          trace.locationCode().c_str(),
                          trace.channelCode().c_str(),
                          samplingFrequency, rec->samplingFrequency());
            return false;
        }

        // Check for gaps and overlaps
        if ( last ) {
            Core::TimeSpan diff = rec->startTime()-last->endTime();
            if ( diff > maxAllowedGap ) {
                SEISCOMP_DEBUG("%s.%s.%s.%s: gap detected of %d.%06ds",
                              trace.networkCode().c_str(),
                              trace.stationCode().c_str(),
                              trace.locationCode().c_str(),
                              trace.channelCode().c_str(),
                              (int)diff.seconds(), (int)diff.microseconds());
                return false;
            }

            if ( diff < maxAllowedOverlap ) {
                SEISCOMP_DEBUG("%s.%s.%s.%s: overlap detected of %fs",
                              trace.networkCode().c_str(),
                              trace.stationCode().c_str(),
                              trace.locationCode().c_str(),
                              trace.channelCode().c_str(),
                              (double)diff);
                return false;
            }
        }

        arr->append( (Array *)(rec->data()));

        last = rec;
    }

    trace.setData(arr.get());

    return true;
}


bool HypoDD::trim(GenericRecord &trace, const Core::TimeWindow& tw) const
{
    int ofs = (int)(double(tw.startTime() - trace.startTime())*trace.samplingFrequency());
    int samples = (int)(tw.length()*trace.samplingFrequency());

    // Not enough data at start of time window
    if ( ofs < 0 )
    {
        SEISCOMP_DEBUG("%s: need %d more samples in past",
                       trace.streamID().c_str(), -ofs);
        return false;
    }

    // Not enough data at end of time window
    if ( ofs+samples > trace.data()->size() )
    {
        SEISCOMP_DEBUG("%s: need %d more samples past the end",
                       trace.streamID().c_str(), trace.data()->size()-samples-ofs);
        return false;
    }

    ArrayPtr sliced = trace.data()->slice(ofs, ofs+samples);

    trace.setStartTime(tw.startTime());
    trace.setData(sliced.get());

    return true;
}


void HypoDD::filter(GenericRecord &trace, bool demeaning,
                    int order, double fmin, double fmax, double resampleFreq) const
{
    DoubleArray *data = DoubleArray::Cast(trace.data());

    if (demeaning)
    {
        double mean = data->mean();
        int cnt = data->size();
        for ( int i = 0; i < cnt; ++i ) (*data)[i] -= mean;
    }

    if ( resampleFreq > 0)
    {
        resample(trace, resampleFreq, true);
    }

    if ( fmin > 0 && fmax > 0 )
    {
        Math::Filtering::IIR::ButterworthHighLowpass<double> bp(order, fmin, fmax);
        bp.setSamplingFrequency(trace.samplingFrequency());
        bp.apply(data->size(), data->typedData());
    }
    else if ( fmin > 0 )
    {
        Math::Filtering::IIR::ButterworthHighpass<double> hp(order, fmin);
        hp.setSamplingFrequency(trace.samplingFrequency());
        hp.apply(data->size(), data->typedData());
    }
    else if ( fmax > 0 )
    {
        Math::Filtering::IIR::ButterworthLowpass<double> lp(order, fmax);
        lp.setSamplingFrequency(trace.samplingFrequency());
        lp.apply(data->size(), data->typedData());
    }
}


void HypoDD::resample(GenericRecord &trace, double sf, bool average) const
{
    if ( sf <= 0 )
        return;

    if ( trace.samplingFrequency() == sf )
        return;

    DoubleArray *data = DoubleArray::Cast(trace.data());
    double step = trace.samplingFrequency() / sf;

    if ( trace.samplingFrequency() < sf ) // upsampling
    {
        double fi = data->size() - 1;
        data->resize( data->size() / step );

        for( int i = data->size()-1; i >= 0; i-- )
        {
            (*data)[i] = (*data)[(int)fi];
            fi -= step;
        }
    }
    else // downsampling
    {
        int w = average?step*0.5 + 0.5:0;
        int i = 0;
        double fi = 0.0;
        int cnt = data->size();

        if ( w <= 0 )
        {
            while ( fi < cnt ) {
                (*data)[i++] = (*data)[(int)fi];
                fi += step;
            }
        }
        else
        {
            while ( fi < cnt )
            {
                int ci = (int)fi;
                double scale = 1.0;
                double v = (*data)[ci];

                for ( int g = 1; g < w; ++g )
                {
                    if ( ci >= g )
                    {
                        v += (*data)[ci-g];
                        scale += 1.0;
                    }

                    if ( ci+g < cnt )
                    {
                        v += (*data)[ci+g];
                        scale += 1.0;
                    }
                }

                v /= scale;

                (*data)[i++] = v;
                fi += step;
            }
        }
        data->resize(i);
    }
    trace.setSamplingFrequency((double)sf);
    trace.dataUpdated();
}


} // HDD
} // Seiscomp

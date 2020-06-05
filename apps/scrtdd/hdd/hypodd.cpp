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
#include "utils.h"

#include <seiscomp3/core/datetime.h>
#include <seiscomp3/core/strings.h>
#include <seiscomp3/core/typedarray.h>
#include <seiscomp3/io/recordinput.h>
#include <seiscomp3/client/inventory.h>
#include <seiscomp3/utils/files.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <boost/filesystem.hpp>
#include <boost/bind.hpp>
#include <boost/range/iterator_range_core.hpp>

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp3/logging/log.h>
#include <seiscomp3/logging/file.h>

using namespace std;
using namespace Seiscomp;
using Seiscomp::Core::stringify;
using DataModel::ThreeComponents;
using Event = HDD::Catalog::Event;
using Phase = HDD::Catalog::Phase;
using Station = HDD::Catalog::Station;
using CacheType = HDD::WfMngr::CacheType;


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

    _tmpCacheDir = (boost::filesystem::path(_workingDir)/"tmpcache").string();
    if ( ! Util::pathExists(_tmpCacheDir) )
    {
        if ( ! Util::createPath(_tmpCacheDir) )
        {
            string msg = "Unable to create cache directory: " + _tmpCacheDir;
            throw runtime_error(msg);
        }
    }

    _wfDebugDir = (boost::filesystem::path(_workingDir)/"wfdebug").string();

    _wf = new WfMngr(_cfg.ddObservations2.recordStreamURL, _cacheDir, _tmpCacheDir, _wfDebugDir);
    _wf->setProcessing(_cfg.wfFilter.filterStr, _cfg.wfFilter.resampleFreq);
    _wf->setSnr(_cfg.snr.minSnr, _cfg.snr.noiseStart, _cfg.snr.noiseEnd, _cfg.snr.signalStart, _cfg.snr.signalEnd);

    setUseCatalogDiskCache(true);
    setWaveformCacheAll(false);
    setWaveformDebug(false);

    _ttt = TravelTimeTableInterface::Create(_cfg.ttt.type.c_str());
    _ttt->setModel(_cfg.ttt.model.c_str());
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
                 ! boost::filesystem::equivalent(entry, _tmpCacheDir) &&
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
    _ddbgc = Catalog::filterPhasesAndSetWeights(_srcCat, Phase::Source::CATALOG,
                                                _cfg.validPphases, _cfg.validSphases);
}


void HypoDD::setWaveformDebug(bool debug)
{
    _waveformDebug = debug;
    if ( _waveformDebug )
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

    _wf->setWaveformDebug(_waveformDebug);
}


// Creates dir name from event. This id has the following format:
// OriginTime_Lat_Lon_CreationDate_Random
// eg 20111210115715_46343_007519_20111210115740_6666
string HypoDD::generateWorkingSubDir(const Event& ev) const
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
        const Event& event = kv.second;
        auto eqlrng = _ddbgc->getPhases().equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
            const Phase& phase = it->second;
            Core::TimeWindow tw = xcorrTimeWindowLong(phase);
            const auto xcorrCfg = _cfg.xcorr.at(phase.procInfo.type);

            for (string component : xcorrCfg.components )
            {
                Phase tmpPh = phase;
                tmpPh.channelCode = WfMngr::getBandAndInstrumentCodes(tmpPh.channelCode) + component;
                _wf->getWaveform(tw, event, tmpPh, &_wfCache, CacheType::PERMANENT, true);
            }

            numPhases++;
            if (  phase.procInfo.type == Phase::Type::S ) numSPhases++;
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


CatalogPtr HypoDD::relocateCatalog()
{
    SEISCOMP_INFO("Starting HypoDD relocator in multiple events mode");

    CatalogPtr catToReloc( new Catalog(*_ddbgc) );

    // Create working directory (used to be a working directory, now it's just logs/debug)
    string catalogWorkingDir = (boost::filesystem::path(_workingDir)/"catalog").string();
    if ( ! Util::pathExists(catalogWorkingDir) )
    {
        if ( !Util::createPath(catalogWorkingDir) )
        {
            string msg = "Unable to create working directory: " + catalogWorkingDir;
            throw runtime_error(msg);
        }
    }

    // prepare file logger
    std::shared_ptr<Logging::Output> processingInfoOutput;
    if ( ! _workingDirCleanup )
    {
        processingInfoOutput = std::shared_ptr<Logging::Output>(new Logging::FileOutput(
            (boost::filesystem::path(catalogWorkingDir)/"info.log").string().c_str()));
        processingInfoOutput->subscribe(Seiscomp::Logging::_SCInfoChannel);
        processingInfoOutput->subscribe(Seiscomp::Logging::_SCWarningChannel);
        processingInfoOutput->subscribe(Seiscomp::Logging::_SCErrorChannel);
    }

    // Find Neighbouring Events in the catalog
    deque< list<NeighboursPtr> > clusters = selectNeighbouringEventsCatalog(
        catToReloc, _cfg.ddObservations2.minWeight,
        _cfg.ddObservations2.minESdist, _cfg.ddObservations2.maxESdist,
        _cfg.ddObservations2.minEStoIEratio, _cfg.ddObservations2.minDTperEvt,
        _cfg.ddObservations2.maxDTperEvt, _cfg.ddObservations2.minNumNeigh,
        _cfg.ddObservations2.maxNumNeigh, _cfg.ddObservations2.numEllipsoids,
        _cfg.ddObservations2.maxEllipsoidSize, true
    );

    SEISCOMP_INFO("Found %lu event clusters", clusters.size() );

    // write catalog for debugging purpose
    if ( ! _workingDirCleanup )
    {
        catToReloc->writeToFile(
            (boost::filesystem::path(catalogWorkingDir)/"starting-event.csv").string(),
            (boost::filesystem::path(catalogWorkingDir)/"starting-phase.csv").string(),
            (boost::filesystem::path(catalogWorkingDir)/"starting-station.csv").string() );
    }

    //
    // Relocate one cluster a time
    //
    CatalogPtr relocatedCatalog( new Catalog() );

    for (unsigned clusterId = 0; clusterId < clusters.size(); clusterId++ )
    {
        const list<NeighboursPtr>& neighCluster = clusters[clusterId];

        SEISCOMP_INFO("Relocating cluster %u (%lu events)", clusterId+1, neighCluster.size() );

        if ( ! _workingDirCleanup )
        {
            CatalogPtr catToDump( new Catalog() );
            for (const NeighboursPtr& n : neighCluster )
                catToDump->add(n->refEvId, *catToReloc, true);
            string prefix = "cluster-" + to_string(clusterId + 1);
            catToDump->writeToFile(
                (boost::filesystem::path(catalogWorkingDir)/(prefix+"-event.csv")).string(),
                (boost::filesystem::path(catalogWorkingDir)/(prefix+"-phase.csv")).string(),
                (boost::filesystem::path(catalogWorkingDir)/(prefix+"-station.csv")).string() ); 
        }

        // Perform cross correlation, which also detects picks around theoretical
        // arrival times. The catalog will be updated with those theoretical phases 
        const XCorrCache xcorr = buildXCorrCache(catToReloc, neighCluster, _useArtificialPhases);

        // The actual relocation
        CatalogPtr relocatedCluster = relocate(catToReloc, neighCluster, false, xcorr);

        relocatedCatalog->add(*relocatedCluster, true);

        if ( ! _workingDirCleanup )
        {
            string prefix = "relocated-cluster-" + to_string(clusterId + 1);
            relocatedCluster->writeToFile(
                (boost::filesystem::path(catalogWorkingDir)/(prefix+"-event.csv")).string(),
                (boost::filesystem::path(catalogWorkingDir)/(prefix+"-phase.csv")).string(),
                (boost::filesystem::path(catalogWorkingDir)/(prefix+"-station.csv")).string() ); 
        }
    }

    // write catalog for debugging purpose
    if ( ! _workingDirCleanup )
    {
        relocatedCatalog->writeToFile(
            (boost::filesystem::path(catalogWorkingDir)/"relocated-event.csv").string(),
            (boost::filesystem::path(catalogWorkingDir)/"relocated-phase.csv").string(),
            (boost::filesystem::path(catalogWorkingDir)/"relocated-station.csv").string());
    }

    if ( _workingDirCleanup ) boost::filesystem::remove_all(catalogWorkingDir);

    return relocatedCatalog;
}



CatalogPtr HypoDD::relocateSingleEvent(const CatalogCPtr& singleEvent)
{
    // there must be only one event in the catalog, the origin to relocate
    Event evToRelocate = singleEvent->getEvents().begin()->second;
    const auto& evToRelocatePhases = singleEvent->getPhases().equal_range(evToRelocate.id);

    SEISCOMP_INFO("Starting HypoDD relocator in single event mode: event %s (%ld phases)",
                  string(evToRelocate).c_str(), std::distance(evToRelocatePhases.first, evToRelocatePhases.second));

    // Create working directory (used to be a working directory, now it's just logs/debug)
    string subFolder = generateWorkingSubDir(evToRelocate);
    subFolder = (boost::filesystem::path(_workingDir)/subFolder).string();
    if ( ! Util::pathExists(subFolder) )
    {
        if ( ! Util::createPath(subFolder) )
        {
            string msg = "Unable to create working directory: " + subFolder;
            throw runtime_error(msg);
        }
    }

    // prepare file logger
    std::shared_ptr<Logging::Output> processingInfoOutput;
    if ( ! _workingDirCleanup )
    {
        processingInfoOutput = std::shared_ptr<Logging::Output>(new Logging::FileOutput(
            (boost::filesystem::path(subFolder)/"info.log").string().c_str()));
        processingInfoOutput->subscribe(Seiscomp::Logging::_SCInfoChannel);
        processingInfoOutput->subscribe(Seiscomp::Logging::_SCWarningChannel);
        processingInfoOutput->subscribe(Seiscomp::Logging::_SCErrorChannel);
    }

    //
    // Step 1: refine location without cross correlation
    //
    SEISCOMP_INFO("Performing step 1: initial location refinement (no cross correlation)");

    string eventWorkingDir = (boost::filesystem::path(subFolder)/"step1").string();

    CatalogPtr evToRelocateCat = Catalog::filterPhasesAndSetWeights(singleEvent,
                                                                    Phase::Source::RT_EVENT,
                                                                    _cfg.validPphases,
                                                                    _cfg.validSphases);

    CatalogPtr relocatedEvCat = relocateEventSingleStep(
            evToRelocateCat, eventWorkingDir, false, false, _cfg.ddObservations1.minWeight,
            _cfg.ddObservations1.minESdist, _cfg.ddObservations1.maxESdist, 
            _cfg.ddObservations1.minEStoIEratio, _cfg.ddObservations1.minDTperEvt, 
            _cfg.ddObservations1.maxDTperEvt, _cfg.ddObservations1.minNumNeigh,
            _cfg.ddObservations1.maxNumNeigh, _cfg.ddObservations1.numEllipsoids,
            _cfg.ddObservations1.maxEllipsoidSize
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
    }

    //
    // Step 2: relocate the refined location this time with cross correlation
    //
    SEISCOMP_INFO("Performing step 2: relocation with cross correlation");

    eventWorkingDir = (boost::filesystem::path(subFolder)/"step2").string();

    CatalogPtr relocatedEvWithXcorr = relocateEventSingleStep(
            evToRelocateCat, eventWorkingDir, true, _useArtificialPhases,
            _cfg.ddObservations2.minWeight, _cfg.ddObservations2.minESdist,
            _cfg.ddObservations2.maxESdist, _cfg.ddObservations2.minEStoIEratio,
            _cfg.ddObservations2.minDTperEvt, _cfg.ddObservations2.maxDTperEvt,
            _cfg.ddObservations2.minNumNeigh, _cfg.ddObservations2.maxNumNeigh,
            _cfg.ddObservations2.numEllipsoids, _cfg.ddObservations2.maxEllipsoidSize
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
    if ( ! Util::createPath(workingDir) )
    {
        string msg = "Unable to create working directory: " + workingDir;
        throw runtime_error(msg);
    }

    CatalogPtr relocatedEvCat;

    try
    {
        // extract event to relocate
        const Event& evToRelocate = evToRelocateCat->getEvents().begin()->second;

        //
        // Select neighbouring events
        //
        bool keepUnmatchedPhases = doXcorr; //useful for detecting missed picks

        NeighboursPtr neighbours = selectNeighbouringEvents(
            _ddbgc, evToRelocate, evToRelocateCat, minPhaseWeight, minESdist,  maxESdist,
            minEStoIEratio, minDTperEvt,  maxDTperEvt, minNumNeigh, maxNumNeigh,
            numEllipsoids, maxEllipsoidSize, keepUnmatchedPhases
        );

        //
        // Prepare catalog to relocate
        //
        CatalogPtr catalog = neighbours->fromNeighbours(_ddbgc);
        unsigned evToRelocateNewId = catalog->add(evToRelocate.id, *evToRelocateCat, false);
        neighbours->refEvId = evToRelocateNewId;

        // write catalog for debugging purpose
        if ( ! _workingDirCleanup )
        {
            catalog->writeToFile(
                (boost::filesystem::path(workingDir)/"starting-event.csv").string(),
                (boost::filesystem::path(workingDir)/"starting-phase.csv").string(),
                (boost::filesystem::path(workingDir)/"starting-station.csv").string());
        }

        XCorrCache xcorr;
        if ( doXcorr )
        {
            // Perform cross correlation, which also detects picks around theoretical
            // arrival times. The catalog will be updated with those theoretical phases 
            xcorr = buildXCorrCache(catalog, {neighbours}, computeTheoreticalPhases);
        }

        // The actual relocation
        relocatedEvCat = relocate(catalog, {neighbours}, true, xcorr);

        // write catalog for debugging purpose
        if ( ! _workingDirCleanup )
        {
            relocatedEvCat->writeToFile(
                (boost::filesystem::path(workingDir)/"relocated-event.csv").string(),
                (boost::filesystem::path(workingDir)/"relocated-phase.csv").string(),
                (boost::filesystem::path(workingDir)/"relocated-station.csv").string());
        }

    } catch ( exception &e ) {
        SEISCOMP_ERROR("%s", e.what());
    }

    return relocatedEvCat;
}


CatalogPtr
HypoDD::relocate(CatalogPtr& catalog,
                 const std::list<NeighboursPtr>& neighCluster,
                 bool keepNeighboursFixed,
                 const XCorrCache& xcorr) const
{
    ObservationParams obsparams;

    for ( unsigned iteration=0; iteration < _cfg.solver.algoIterations; iteration++ )
    {
        //
        // compute parameters for this loop iteration
        //
        auto interpolate = [&](double start, double end) -> double
        {
            if ( _cfg.solver.algoIterations < 2 ) return (start + end) / 2;
            return start + (end - start) * iteration / (_cfg.solver.algoIterations-1);
        };

        double dampingFactor = interpolate(_cfg.solver.dampingFactorStart,
                                           _cfg.solver.dampingFactorEnd);
        double downWeightingByResidual = interpolate(_cfg.solver.downWeightingByResidualStart,
                                                     _cfg.solver.downWeightingByResidualEnd); 
        double meanLonShiftConstraint = interpolate(_cfg.solver.meanShiftConstraintStart[0],
                                                    _cfg.solver.meanShiftConstraintEnd[0]);
        double meanLatShiftConstraint = interpolate(_cfg.solver.meanShiftConstraintStart[1],
                                                    _cfg.solver.meanShiftConstraintEnd[1]); 
        double meanDepthShiftConstraint = interpolate(_cfg.solver.meanShiftConstraintStart[2],
                                                      _cfg.solver.meanShiftConstraintEnd[2]); 
        double meanTTShiftConstraint = interpolate(_cfg.solver.meanShiftConstraintStart[3],
                                                    _cfg.solver.meanShiftConstraintEnd[3]);
        double absTTDiffObsWeight = interpolate(1.0, _cfg.solver.absTTDiffObsWeight);
        double xcorrObsWeight     = interpolate(1.0, _cfg.solver.xcorrObsWeight);

        SEISCOMP_INFO("Solving iteration %u num events %lu (observ. weight TTdiff %.2f xcorr %.2f "
                      "dampingFactor %.2f downWeightingByResidual %.2f "
                      "meanShiftConstrainWeight %.2f,%.2f,%.2f,%.2f)",
                      iteration, neighCluster.size(), absTTDiffObsWeight, xcorrObsWeight,
                      dampingFactor, downWeightingByResidual,
                      meanLonShiftConstraint, meanLatShiftConstraint,
                      meanDepthShiftConstraint, meanTTShiftConstraint); 

        // Create a solver and then add observations
        Solver solver(_cfg.solver.type);

        //
        // Add absolute travel time/xcorr differences to the solver (the observations)
        //
        for (const NeighboursPtr& neighbours : neighCluster)
        {
            addObservations(solver, absTTDiffObsWeight, xcorrObsWeight, catalog,
                            neighbours, keepNeighboursFixed, xcorr, obsparams);
        }
        obsparams.addToSolver(solver);


        //
        // solve the system
        //
        try {
            solver.solve(_cfg.solver.solverIterations, dampingFactor,
                         downWeightingByResidual, meanLonShiftConstraint,
                         meanLatShiftConstraint,  meanDepthShiftConstraint,
                         meanTTShiftConstraint, _cfg.solver.L2normalization);
        } catch ( exception &e ) {
            SEISCOMP_INFO("Cannot solve the double-difference system, stop here (%s)", e.what());
            break;
        }

        // prepare for next iteration
        obsparams = ObservationParams();

        // update event parameters
        catalog = updateRelocatedEvents(solver, catalog, neighCluster, obsparams);
    }

    // build the relocated catalog from the results of relocations
    CatalogPtr relocatedCatalog( new Catalog() );
    for (const NeighboursPtr& neighbours : neighCluster)
    {
        const Event& event = catalog->getEvents().at(neighbours->refEvId);
        if ( event.relocInfo.isRelocated )
        {
            relocatedCatalog->add(neighbours->refEvId, *catalog, true);
        }
    }

    return relocatedCatalog;
}


string HypoDD::relocationReport(const CatalogCPtr& relocatedEv)
{
    Event event = relocatedEv->getEvents().begin()->second;
    if ( ! event.relocInfo.isRelocated )
        return "Event not relocated";

    return stringify("Rms residual %.3f [sec]. Neighboring events %u used P %u used S %u."
                     "Num observations: %u (avg. a priory weight %.2f final weight %.2f). "
                     "Num xcross observations: P phases %u S phases %u. "
                     "Num absolute TT diff observations: P phases %u S phases %u",
                      event.rms, event.relocInfo.numNeighbours, 
                      event.relocInfo.usedP, event.relocInfo.usedS,
                      (event.relocInfo.numCCp+event.relocInfo.numCCs+
                      event.relocInfo.numCTp+event.relocInfo.numCTs),
                      event.relocInfo.meanObsWeight, event.relocInfo.meanFinalObsWeight,
                      event.relocInfo.numCCp, event.relocInfo.numCCs,
                      event.relocInfo.numCTp, event.relocInfo.numCTs);
}


/*
 * Add to the Solver the absolute travel times differences and the
 * differential travel times from cross correlation for pairs of
 * earthquakes
 */ 
void
HypoDD::addObservations(Solver& solver, double absTTDiffObsWeight, double xcorrObsWeight,
                        const CatalogCPtr& catalog, const NeighboursPtr& neighbours,
                        bool keepNeighboursFixed, const XCorrCache& xcorr,
                        ObservationParams& obsparams ) const
{
    // copy event because we'll update it
    const Event& refEv = catalog->getEvents().at(neighbours->refEvId);

    //
    // loop through reference event phases
    //
    auto eqlrng = catalog->getPhases().equal_range(refEv.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
        const Phase& refPhase = it->second;
        const Station& station = catalog->getStations().at(refPhase.stationId);
        char phaseTypeAsChar = static_cast<char>(refPhase.procInfo.type);

        //
        // loop through neighbouring events and look for the matching phase
        // 
        for ( unsigned neighEvId : neighbours->ids )
        {
            const Event& event = catalog->getEvents().at(neighEvId);

            if ( ! neighbours->has(neighEvId, refPhase.stationId, refPhase.procInfo.type) )
                continue;

            const Phase& phase = catalog->searchPhase(event.id, refPhase.stationId,
                                                      refPhase.procInfo.type)->second;

            //
            // compute travel times for both event and refEvent
            //
            double ref_travel_time = (refPhase.time - refEv.time).length();
            if (ref_travel_time < 0)
            {
                SEISCOMP_DEBUG("Ignoring phase %s with negative travel time",
                               string(refPhase).c_str());
                continue;
            }

            double travel_time = (phase.time - event.time).length();
            if (travel_time < 0)
            {
                SEISCOMP_DEBUG("Ignoring phase %s with negative travel time",
                               string(phase).c_str());
                continue;
            }

            try {
                obsparams.add(_ttt, refEv, station, phaseTypeAsChar);
                obsparams.add(_ttt, event, station, phaseTypeAsChar);
            } catch ( exception &e ) {
                SEISCOMP_DEBUG("Skipping observation (ev %u-%u sta %s phase %c): %s",
                    refEv.id, event.id, station.id.c_str(), phaseTypeAsChar, e.what());
                continue;
            }

            //
            // Conpute absolute trave time differences to the solver
            //
            double diffTime = ref_travel_time - travel_time;
            double weight   =  _cfg.solver.usePickUncertainty
                            ? (refPhase.procInfo.weight + phase.procInfo.weight) / 2.0
                            : 1.0;
            bool isXcorr = false;

            //
            // Check if we have xcorr results for current event/refEvent pair at station/phase
            // and use those instead
            //  
            if ( xcorr.has(refEv.id, event.id, refPhase.stationId, refPhase.procInfo.type) )
            {

                const auto& xcdata = xcorr.get(refEv.id, event.id, refPhase.stationId, refPhase.procInfo.type);
                diffTime -= xcdata.lag;
                weight *= xcorrObsWeight;
                isXcorr = true;
            }
            else
            {
                weight *= absTTDiffObsWeight;
            } 

            solver.addObservation(refEv.id, event.id, refPhase.stationId,  phaseTypeAsChar,
                                  diffTime, weight, true, !keepNeighboursFixed, isXcorr);
        }
    }
}

void
HypoDD::ObservationParams::add(TravelTimeTableInterfacePtr ttt, const Event& event,
                               const Station& station, char phaseType )
{
    const std::string key = std::to_string(event.id) + "@" + station.id + ":" + phaseType;
    if ( _entries.find(key) == _entries.end() )
    {
        TravelTime tt = ttt->compute(string(1, phaseType).c_str(),
                                     event.latitude, event.longitude, event.depth, 
                                     station.latitude, station.longitude, station.elevation);
        _entries[key] = Entry( {event, station, phaseType, tt.time} );
    }
}

const HypoDD::ObservationParams::Entry&
HypoDD::ObservationParams::get(unsigned eventId, const std::string stationId, char phaseType ) const
{
    const std::string key = std::to_string(eventId) + "@" + stationId + ":" + phaseType;
    return _entries.at(key);
}


void
HypoDD::ObservationParams::addToSolver(Solver& solver) const
{
    for ( const auto& kv : _entries )
    {
        const ObservationParams::Entry& e = kv.second;
        solver.addObservationParams(e.event.id, e.station.id, e.phaseType,
                                    e.event.latitude, e.event.longitude, e.event.depth,
                                    e.station.latitude, e.station.longitude, e.station.elevation,
                                    e.travelTime);
    }
}


CatalogPtr
HypoDD::updateRelocatedEvents(const Solver& solver,
                              const CatalogCPtr& catalog,
                              const std::list<NeighboursPtr>& neighCluster,
                              ObservationParams& obsparams ) const
{
    unordered_map<string,Station> stations = catalog->getStations();
    map<unsigned,Event> events = catalog->getEvents();
    unordered_multimap<unsigned,Phase> phases = catalog->getPhases();
    unsigned relocatedEvs = 0;

    for (const NeighboursPtr& neighbours : neighCluster)
    {
        Event& event = events.at(neighbours->refEvId);
        event.relocInfo.isRelocated = false;

        double deltaLat, deltaLon, deltaDepth, deltaTT;
        if ( ! solver.getEventChanges(event.id, deltaLat, deltaLon, deltaDepth, deltaTT) )
        {
            continue;
        }

        if ( event.depth + deltaDepth < 0 )
        {
            SEISCOMP_DEBUG("Ignoring airquake event %s", string(event).c_str());
            continue;
        }

        relocatedEvs++;

        event.latitude  += deltaLat;
        event.longitude += deltaLon;
        event.depth     += deltaDepth;
        event.time      += Core::TimeSpan(deltaTT);
        event.relocInfo.isRelocated = true;
        event.relocInfo.numNeighbours = neighbours->numNeighbours;

        event.relocInfo.usedP = 0;
        event.relocInfo.usedS = 0;
        event.relocInfo.numCCp = 0;
        event.relocInfo.numCCs = 0;
        event.relocInfo.numCTp = 0;
        event.relocInfo.numCTs = 0; 
        event.relocInfo.meanObsWeight = 0;
        event.relocInfo.meanFinalObsWeight = 0;
        event.rms = 0.;

        unsigned rmsCount = 0;
        unsigned pCount = 0;
        auto eqlrng = phases.equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it) 
        {
            Phase& phase = it->second;
            const Station& station = stations.at(phase.stationId);
            char phaseTypeAsChar = static_cast<char>(phase.procInfo.type);

            phase.relocInfo.isRelocated = false;

            unsigned startingObservations, startingXcorrObservations, totalFinalObservations;
            double meanAPrioriWeight, meanFinalWeight;
            if ( ! solver.getObservationParamsChanges(event.id, station.id, phaseTypeAsChar,
                                                      startingObservations, 
                                                      startingXcorrObservations,
                                                      totalFinalObservations,
                                                      meanAPrioriWeight,
                                                      meanFinalWeight) )
            {
                continue;
            }

            event.relocInfo.meanObsWeight      += meanAPrioriWeight;
            event.relocInfo.meanFinalObsWeight += meanFinalWeight;
            ++pCount;

            if (  meanFinalWeight == 0 || totalFinalObservations == 0 )
                continue;

            phase.relocInfo.isRelocated = true;
            phase.relocInfo.finalWeight = phase.procInfo.weight * meanFinalWeight / meanAPrioriWeight;
            phase.relocInfo.numObservs           = startingObservations;
            phase.relocInfo.numXcorrObservs      = startingXcorrObservations;
            phase.relocInfo.meanObsWeight        = meanAPrioriWeight;
            phase.relocInfo.meanFinalObsWeight   = meanFinalWeight;

            try {
                obsparams.add(_ttt, event, station, phaseTypeAsChar);
                double travelTime = obsparams.get(event.id, station.id, phaseTypeAsChar).travelTime;
                phase.relocInfo.residual = travelTime - (phase.time - event.time).length();
                rmsCount++;
            } catch ( exception &e ) { phase.relocInfo.residual = 0; }

            event.rms += (phase.relocInfo.residual * phase.relocInfo.residual);
            if ( phase.procInfo.type == Phase::Type::P ) 
            {
                event.relocInfo.usedP++;
                event.relocInfo.numCCp        += phase.relocInfo.numXcorrObservs;
                event.relocInfo.numCTp        += phase.relocInfo.numObservs;
            }
            if ( phase.procInfo.type == Phase::Type::S )
            {
                event.relocInfo.usedS++;
                event.relocInfo.numCCs        += phase.relocInfo.numXcorrObservs;
                event.relocInfo.numCTs        += phase.relocInfo.numObservs;
            }
        }

        if ( rmsCount > 0 ) event.rms = std::sqrt(event.rms / rmsCount);
        if ( pCount > 0 )
        {
            event.relocInfo.meanObsWeight      /= pCount;
            event.relocInfo.meanFinalObsWeight /= pCount;
        }
    }

    SEISCOMP_INFO("Successfully relocated %u events", relocatedEvs);

    return new Catalog(stations, events, phases);
}


void
HypoDD::addMissingEventPhases(const Event& refEv, CatalogPtr& refEvCatalog,
                              const CatalogCPtr& searchCatalog,
                              const NeighboursPtr& neighbours)
{
    std::vector<Phase> newPhases = findMissingEventPhases(refEv, refEvCatalog,
                                                          searchCatalog, neighbours);

    for (Phase& ph : newPhases)
    {
        refEvCatalog->updatePhase(ph, true);
        const Station& station = searchCatalog->getStations().at(ph.stationId);
        refEvCatalog->addStation(station);
    }
}



std::vector<Phase>
HypoDD::findMissingEventPhases(const Event& refEv, CatalogPtr& refEvCatalog,
                               const CatalogCPtr& searchCatalog,
                               const NeighboursPtr& neighbours)
{
    //
    // find stations for which the refEv doesn't have phases
    //
    vector<MissingStationPhase> missingPhases = getMissingPhases(refEv, refEvCatalog, 
                                                                 searchCatalog);

    //
    // for each missed phase try to detect it
    //
    std::vector<Phase> newPhases;
    for ( const MissingStationPhase& pair : missingPhases )
    {
        const Station& station = searchCatalog->getStations().at(pair.first);
        const Phase::Type phaseType = pair.second;

        //
        // loop through each other event and select the ones who have a manually picked phase for the missing station
        //
        vector<HypoDD::PhasePeer> peers = findPhasePeers(station, phaseType,
                                                         searchCatalog, neighbours);
        if ( peers.size() <= 0 )
        {
            continue;
        }

        // compute velocity using existing background catalog phases
        double phaseVelocity = 0;
        for (const PhasePeer& peer  : peers)
        {
            const Event& event = peer.first;
            const Phase& phase = peer.second;
            double travelTime = (phase.time - event.time).length();
            double stationDistance = computeDistance(event, station);
            double vel = stationDistance / travelTime;
            phaseVelocity += vel;
        }
        phaseVelocity /= peers.size();

        Phase refEvNewPhase = createThoreticalPhase(station, phaseType, refEv, peers, phaseVelocity);

        newPhases.push_back(refEvNewPhase);
    }

    return newPhases;
}



vector<HypoDD::MissingStationPhase>
HypoDD::getMissingPhases(const Event& refEv, CatalogPtr& refEvCatalog,
                         const CatalogCPtr& searchCatalog) const
{
    const auto& refEvPhases = refEvCatalog->getPhases().equal_range(refEv.id);

    //
    // loop through stations and find those for which the refEv doesn't have phases
    //
    vector<MissingStationPhase> missingPhases;
    for (const auto& kv : searchCatalog->getStations() )
    {
        const Station& station = kv.second;

        bool foundP = false, foundS = false;
        for (auto it = refEvPhases.first; it != refEvPhases.second; ++it)
        {
            const Phase& phase = it->second;
            if ( station.networkCode  == phase.networkCode &&
                 station.stationCode  == phase.stationCode &&
                 station.locationCode == phase.locationCode)
            {
                if ( phase.procInfo.type == Phase::Type::P ) foundP = true;
                if ( phase.procInfo.type == Phase::Type::S ) foundS = true;
            }
            if ( foundP and foundS ) break;
        }
        if ( ! foundP || ! foundS )
        {
            if ( ! foundP )
                missingPhases.push_back( MissingStationPhase(station.id,Phase::Type::P) );
            if ( ! foundS )
                missingPhases.push_back( MissingStationPhase(station.id,Phase::Type::S) );
        }
    }

    return missingPhases;
}



vector<HypoDD::PhasePeer>
HypoDD::findPhasePeers(const Station& station, const Phase::Type& phaseType,
                       const CatalogCPtr& searchCatalog, 
                       const NeighboursPtr& neighbours) const
{
    //
    // loop through each other event and select the manual phases for the station we
    // are interested in
    //
    vector<PhasePeer> phasePeers;

    for ( unsigned neighEvId : neighbours->ids )
    {
        const Event& event = searchCatalog->getEvents().at(neighEvId);

        if ( neighbours->has(neighEvId, station.id, phaseType) )
        {
            const Phase& phase = searchCatalog->searchPhase(neighEvId, station.id, phaseType)->second;

            if ( station.networkCode  == phase.networkCode  &&
                 station.stationCode  == phase.stationCode  &&
                 station.locationCode == phase.locationCode &&
                 phaseType            == phase.procInfo.type )
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



Phase
HypoDD::createThoreticalPhase(const Station& station,
                              const Phase::Type& phaseType,
                              const Event& refEv,
                              const vector<HypoDD::PhasePeer>& peers,
                              double phaseVelocity)
{
    const auto xcorrCfg = _cfg.xcorr.at(phaseType);

    // store most recent channelCode used
    struct {
        string channelCode;
        Core::Time time;
    } streamInfo = {"", Core::Time()};

    for (const PhasePeer& peer : peers)
    {
        //const Event& event = peer.first;
        const Phase& phase = peer.second;
        // get the closest in time to refEv stream information
        if ( (refEv.time - phase.time).abs() < (refEv.time - streamInfo.time).abs() )
            streamInfo = {phase.channelCode, phase.time};
    }

    // initialize the new phase
    Phase refEvNewPhase;

    refEvNewPhase.eventId      = refEv.id;
    refEvNewPhase.stationId    = station.id;
    refEvNewPhase.networkCode  = station.networkCode;
    refEvNewPhase.stationCode  = station.stationCode;
    refEvNewPhase.locationCode = station.locationCode;
    refEvNewPhase.channelCode  = WfMngr::getBandAndInstrumentCodes(streamInfo.channelCode) + xcorrCfg.components[0];
    refEvNewPhase.isManual     = false;
    refEvNewPhase.procInfo.type = phaseType;

    // use phase velocity to compute phase time
    double stationDistance = computeDistance(refEv, station);
    refEvNewPhase.time = refEv.time + Core::TimeSpan(stationDistance / phaseVelocity);

    refEvNewPhase.lowerUncertainty = Catalog::DEFAULT_AUTOMATIC_PICK_UNCERTAINTY;
    refEvNewPhase.upperUncertainty = Catalog::DEFAULT_AUTOMATIC_PICK_UNCERTAINTY;
    refEvNewPhase.procInfo.weight = Catalog::computePickWeight(refEvNewPhase);
    refEvNewPhase.procInfo.source = Phase::Source::THEORETICAL;
    refEvNewPhase.type = stringify("%ct", static_cast<char>(phaseType));

    return refEvNewPhase;
}


XCorrCache
HypoDD::buildXCorrCache(CatalogPtr& catalog,
                        const std::list<NeighboursPtr>& neighCluster,
                        bool computeTheoreticalPhases)
{
    XCorrCache xcorr;
    _counters = {0};
    _wf->resetCounters();

    for (const NeighboursPtr& neighbours : neighCluster)
    {
        const Event& refEv = catalog->getEvents().at(neighbours->refEvId);

        // Compute theoretical phases for stations that have no picks. The cross correlation will
        // be used to detect and fix pick time
        if ( computeTheoreticalPhases )
        {
            addMissingEventPhases(refEv, catalog, catalog, neighbours);
        }

        buildXcorrDiffTTimePairs(catalog, neighbours, refEv, xcorr);

        // Update theoretical and automatic phase pick time and uncertainties based on
        // cross-correlation results
        // Also drop theoretical phases wihout any good cross correlation result
        fixPhases(catalog, refEv, xcorr);
    }

    printCounters();

    return xcorr;
}


/*
 * Compute and store to XCorrCache cross-correlated differential travel times
 * for pairs of earthquake
 */
void 
HypoDD::buildXcorrDiffTTimePairs(CatalogPtr& catalog,
                                 const NeighboursPtr& neighbours,
                                 const Event& refEv,
                                 XCorrCache& xcorr)
{
    SEISCOMP_INFO("Computing cross-correlation differential travel times for event %s",
                  string(refEv).c_str() );

    WfMngr::WfCache wfTmpCache;

    CacheType permCache = _useCatalogDiskCache ? CacheType::PERMANENT : CacheType::NONE;
    CacheType tempCache = (_useCatalogDiskCache && _waveformCacheAll) ? CacheType::TEMP : CacheType::NONE; 

    // xcorr settings depending on the phase type
    map<Phase::Source, PhaseXCorrCfg> phCfgs = {
        {Phase::Source::CATALOG,      {permCache, &_wfCache,}},
        {Phase::Source::RT_EVENT,     {tempCache, &wfTmpCache,}},
        {Phase::Source::THEORETICAL,  {tempCache, &wfTmpCache,}}
    };

    // keep track of refEv distant to stations
    multimap<double,string> stationByDistance; // <distance, stationid>
    unordered_set<string> computedStations;

    //
    // loop through reference event phases
    //
    auto eqlrngRef = catalog->getPhases().equal_range(refEv.id);
    for (auto itRef = eqlrngRef.first; itRef != eqlrngRef.second; ++itRef)
    {
        const Phase& refPhase = itRef->second;
        const Station& station = catalog->getStations().at(refPhase.stationId);

        //
        // skip stations too far away
        //
        double stationDistance = computeDistance(refEv, station);
        if ( stationDistance > _cfg.ddObservations2.xcorrMaxEvStaDist &&
             _cfg.ddObservations2.xcorrMaxEvStaDist >= 0 )
            continue;

        //
        // loop through neighbouring events and cross correlate phase pairs
        //
        for ( unsigned neighEvId : neighbours->ids )
        {
            const Event& event = catalog->getEvents().at(neighEvId);

            //
            // skip events too far away
            //
            double interEventDistance = computeDistance(refEv, event);
            if ( interEventDistance > _cfg.ddObservations2.xcorrMaxInterEvDist &&
                 _cfg.ddObservations2.xcorrMaxInterEvDist >= 0 )
                continue;

            if ( neighbours->has(neighEvId, refPhase.stationId, refPhase.procInfo.type) )
            {
                const Phase& phase = catalog->searchPhase(event.id, refPhase.stationId,
                                                          refPhase.procInfo.type)->second;

                // for non-manual phases the SNR is checked later, once the pick time
                // has been fixed using xcorr results (we also trust pick times for all
                // catalog phases)
                PhaseXCorrCfg refPhCfg = phCfgs.at(refPhase.procInfo.source);
                refPhCfg.allowSnrCheck = refPhase.isManual || 
                                         (refPhase.procInfo.source == Phase::Source::CATALOG);

                // Catalog phases always allow SNR check, since will not fix those
                PhaseXCorrCfg phaseCfg = phCfgs.at(phase.procInfo.source);
                phaseCfg.allowSnrCheck = true;

                double coeff, lag;
                if ( xcorrPhases(refEv, refPhase, refPhCfg, event, phase, phaseCfg, coeff, lag) )
                {
                    //
                    // Store good xcorr results
                    //
                    auto& entry = xcorr.getForUpdate(refEv.id, refPhase.stationId, refPhase.procInfo.type);
                    entry.update(event, phase, coeff, lag);
                }

                // keep trace of  events/station distance for every xcorr performed
                if ( computedStations.find(refPhase.stationId) == computedStations.end() )
                {
                    stationByDistance.emplace(stationDistance, refPhase.stationId);
                    computedStations.insert(refPhase.stationId);
                }
            }
        }

        if ( xcorr.has(refEv.id, refPhase.stationId, refPhase.procInfo.type) )
        {
            // finalize statistics
            auto& entry = xcorr.getForUpdate(refEv.id, refPhase.stationId, refPhase.procInfo.type);
            entry.computeStats();

            // discard phases with low SNR (if not already done at the previous step)
            if ( _cfg.snr.minSnr > 0 && ! refPhase.isManual && 
                (refPhase.procInfo.source != Phase::Source::CATALOG) )
            {
                const auto xcorrCfg = _cfg.xcorr.at(refPhase.procInfo.type);
                Core::TimeWindow tw = xcorrTimeWindowShort(refPhase);

                GenericRecordCPtr trace;
                for (const string& component : xcorrCfg.components )
                {
                    Phase tmpPh = refPhase;
                    tmpPh.time  -= Core::TimeSpan(entry.mean_lag);
                    tmpPh.channelCode = WfMngr::getBandAndInstrumentCodes(tmpPh.channelCode) + component;
                    trace = _wf->getWaveform(tw, refEv, tmpPh, nullptr, tempCache, true);
                    if ( trace ) break;
                }

                if ( ! trace )
                {
                    xcorr.remove(refEv.id, refPhase.stationId, refPhase.procInfo.type);
                }
            }
        }
    }

    // Print some useful information
    for (const auto& kv : stationByDistance)
    {
        const double stationDistance = kv.first;
        const Station& station = catalog->getStations().at(kv.second);

        bool goodPXcorr = xcorr.has(refEv.id, station.id, Phase::Type::P);
        bool goodSXcorr = xcorr.has(refEv.id, station.id, Phase::Type::S);

        if ( ! goodPXcorr && ! goodSXcorr )
        {
            SEISCOMP_INFO("xcorr: event %5s sta %4s %5s dist %7.2f [km] - low corr coeff pairs",
                          string(refEv).c_str(), station.networkCode.c_str(),
                          station.stationCode.c_str(), stationDistance);
        }
        else
        {
            if ( goodPXcorr )
            {
                const auto& pdata = xcorr.get(refEv.id, station.id, Phase::Type::P);
                SEISCOMP_INFO("xcorr: event %5s sta %4s %5s dist %7.2f [km] - "
                          "%d P phases, mean coeff %.2f lag %.2f (events: %s)",
                          string(refEv).c_str(), station.networkCode.c_str(),
                          station.stationCode.c_str(), stationDistance,
                          pdata.ccCount, pdata.mean_coeff, pdata.mean_lag,
                          pdata.peersStr.c_str());
            }
            if ( goodSXcorr )
            {
                const auto& sdata = xcorr.get(refEv.id, station.id, Phase::Type::S);
                SEISCOMP_INFO("xcorr: event %5s sta %4s %5s dist %7.2f [km] - "
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
 * cross-correlation results.
 * Also drop theoretical phases wihout any good cross correlation result
 */
void HypoDD::fixPhases(CatalogPtr& catalog,
                       const Event& refEv,
                       XCorrCache& xcorr)
{
    unsigned newP = 0, newS = 0;

    std::vector<Phase> phasesToBeRemoved;
    std::vector<Phase> newPhases;

    auto eqlrng = catalog->getPhases().equal_range(refEv.id);
    for (auto it = eqlrng.first; it != eqlrng.second; it++)
    {
        const Phase& phase = it->second;
        bool goodXcorr = xcorr.has(refEv.id, phase.stationId, phase.procInfo.type);

        // nothing to do if we dont't have good xcorr results of if the phase is manual or from catalog
        if ( ! goodXcorr || phase.isManual || (phase.procInfo.source == Phase::Source::CATALOG) )
        {
            // remove thoretical phases wihtout good xcorr results
            if ( phase.procInfo.source == Phase::Source::THEORETICAL )
                phasesToBeRemoved.push_back(phase);
            continue;
        }

        const auto& pdata = xcorr.get(refEv.id, phase.stationId, phase.procInfo.type);

        //
        // Set new phase time and uncertainty
        //
        Phase newPhase(phase);
        newPhase.time  -= Core::TimeSpan(pdata.mean_lag);
        newPhase.lowerUncertainty = pdata.mean_lag - pdata.min_lag;
        newPhase.upperUncertainty = pdata.max_lag - pdata.mean_lag;
        newPhase.procInfo.weight = Catalog::computePickWeight(newPhase);
        newPhase.procInfo.source = Phase::Source::XCORR;
        newPhase.type = stringify("%cx", static_cast<char>(newPhase.procInfo.type));

        if ( phase.procInfo.source == Phase::Source::THEORETICAL )
        {
            newP += newPhase.procInfo.type == Phase::Type::P ? 1 : 0;
            newS += newPhase.procInfo.type == Phase::Type::S ? 1 : 0; 
        }

        newPhases.push_back(newPhase);
    }

    //
    // Replace automatic/theoretical phases with xcorr detected ones
    //
    for (const Phase& ph : phasesToBeRemoved)
    {
        catalog->removePhase(ph.eventId, ph.stationId, ph.procInfo.type);
    }

    for (Phase& ph : newPhases)
    {
        catalog->updatePhase(ph, true);
    }

    const auto& refEvPhases = catalog->getPhases().equal_range(refEv.id);
    SEISCOMP_INFO("Event %s total phases %lu: created %u (%u P and %u S) from theoretical picks",
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
HypoDD::xcorrTimeWindowLong(const Phase& phase) const
{
    const auto xcorrCfg = _cfg.xcorr.at(phase.procInfo.type);
    Core::TimeWindow tw = xcorrTimeWindowShort(phase);
    tw.setStartTime( tw.startTime() - Core::TimeSpan(xcorrCfg.maxDelay) );
    tw.setEndTime(   tw.endTime()   + Core::TimeSpan(xcorrCfg.maxDelay) );
    return tw;
}


Core::TimeWindow
HypoDD::xcorrTimeWindowShort(const Phase& phase) const
{
    const auto xcorrCfg = _cfg.xcorr.at(phase.procInfo.type);
    double shortDuration = xcorrCfg.endOffset - xcorrCfg.startOffset;
    Core::TimeSpan shortTimeCorrection = Core::TimeSpan(xcorrCfg.startOffset);
    return Core::TimeWindow(phase.time + shortTimeCorrection, shortDuration);
}


bool
HypoDD::xcorrPhases(const Event& event1, const Phase& phase1, PhaseXCorrCfg& phCfg1,
                    const Event& event2, const Phase& phase2, PhaseXCorrCfg& phCfg2,
                    double& coeffOut, double& lagOut)
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
    // Try to use the same channels in cross correlation, in case the two phases differ
    // but do not change the catalog phase channels
    //
    const string channelCodeRoot1 = WfMngr::getBandAndInstrumentCodes(phase1.channelCode);
    const string channelCodeRoot2 = WfMngr::getBandAndInstrumentCodes(phase2.channelCode);

    string commonChRoot;

    if ( channelCodeRoot1 == channelCodeRoot2 )
    {
        commonChRoot = channelCodeRoot1;
    }
    else if (phase1.procInfo.source == Phase::Source::CATALOG &&
             phase2.procInfo.source != Phase::Source::CATALOG )
    {
        DataModel::ThreeComponents dummy;
        DataModel::SensorLocation *loc2 = Catalog::findSensorLocation(phase2.networkCode, phase2.stationCode, phase2.locationCode, phase2.time);
        if ( loc2 && getThreeComponents(dummy, loc2, channelCodeRoot1.c_str(), phase2.time) )
        {
            // phase 2 has the same channels of phase 1
            commonChRoot = channelCodeRoot1;
        } 
    }
    else if (phase1.procInfo.source != Phase::Source::CATALOG &&
             phase2.procInfo.source == Phase::Source::CATALOG )
    { 
        DataModel::ThreeComponents dummy;
        DataModel::SensorLocation *loc1 = Catalog::findSensorLocation(phase1.networkCode, phase1.stationCode, phase1.locationCode, phase1.time);
        if ( loc1 && getThreeComponents(dummy, loc1, channelCodeRoot2.c_str(), phase1.time) )
        {
            // phase 1 has the same channels of phase 2
            commonChRoot = channelCodeRoot2;
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
        Phase tmpPh1 = phase1;
        Phase tmpPh2 = phase2;

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

        performed = _xcorrPhases(event1, tmpPh1, phCfg1, event2, tmpPh2, phCfg2, coeffOut, lagOut);

        coeffOut = std::abs(coeffOut);

        goodCoeff = ( performed && coeffOut >= xcorrCfg.minCoef );

        // if the xcorr was successfull and the coeffiecnt is good then stopi here
        if ( goodCoeff )
        {
            break;
        }
    }

    //
    // Deal with counters
    //
    bool isS = ( phase1.procInfo.type == Phase::Type::S );
    bool isTheoretical = ( phase1.procInfo.source == Phase::Source::XCORR       ||
                           phase2.procInfo.source == Phase::Source::XCORR       ||
                           phase1.procInfo.source == Phase::Source::THEORETICAL ||
                           phase2.procInfo.source == Phase::Source::THEORETICAL );

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
HypoDD::_xcorrPhases(const Event& event1, const Phase& phase1, PhaseXCorrCfg& phCfg1,
                     const Event& event2, const Phase& phase2, PhaseXCorrCfg& phCfg2,
                     double& coeffOut, double& lagOut)
{
    coeffOut = lagOut = 0;

    auto xcorrCfg = _cfg.xcorr.at(phase1.procInfo.type);

    Core::TimeWindow tw1 = xcorrTimeWindowLong(phase1);
    Core::TimeWindow tw2 = xcorrTimeWindowLong(phase2);

    // load the long trace 1, because we want to cache the long version. Then we'll trim it.
    GenericRecordCPtr tr1 = _wf->getWaveform(tw1, event1, phase1, phCfg1.cache, phCfg1.type, phCfg1.allowSnrCheck);
    if ( !tr1 )
    {
        return false;
    }

    // load the long trace 2, because we want to cache the long version. Then we'll trim it
    GenericRecordCPtr tr2 = _wf->getWaveform(tw2, event2, phase2, phCfg2.cache, phCfg2.type, phCfg2.allowSnrCheck);
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

    coeffOut  = xcorr_coeff;
    lagOut    = xcorr_lag;

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
            if ( coeff < prevCoeff && notDecreasing )
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
            denomS += smpsS[idxS] * smpsS[idxS];

            int idxL = idxS + (smpsLsize-smpsSsize)/2 + delay;
            if (idxL < 0 || idxL >= smpsLsize)
                continue;

            numer  += smpsS[idxS] * smpsL[idxL];
            denomL += smpsL[idxL] * smpsL[idxL];
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


namespace {

    struct XCorrEvalStats {
        unsigned total = 0;
        unsigned goodCC = 0;
        double deviation = 0;
        double absDeviation = 0;
        double meanCoeff = 0;
        unsigned meanCount = 0;

        XCorrEvalStats& operator+=(XCorrEvalStats const& rhs)&
        {
          total += rhs.total;
          goodCC += rhs.goodCC;
          deviation += rhs.deviation;
          absDeviation += rhs.absDeviation;
          meanCoeff += rhs.meanCoeff;
          meanCount += rhs.meanCount;
          return *this;
        }

        friend XCorrEvalStats operator+(XCorrEvalStats lhs, XCorrEvalStats const& rhs)
        {
          lhs+=rhs;
          return lhs;
        }

        void normalize()
        {
            if ( goodCC != 0 )
            {
              deviation    /= goodCC;
              absDeviation /= goodCC;
              meanCoeff    /= goodCC;
              meanCount    /= goodCC;
            }
        }

        string describe(bool theoretical) const
        {
            XCorrEvalStats tmp = *this;
            tmp.normalize();
            string log = stringify("#pha %6d pha good CC %3.f%% avg coeff %.2f avg #matches/ph %2d",
                              tmp.total, (tmp.goodCC * 100. / tmp.total),
                              tmp.meanCoeff, tmp.meanCount); 
            if ( theoretical )
            {
                log += stringify(" time-diff %3.f [msec] abs time-diff %3.f [msec]",
                                 tmp.deviation*1000, tmp.absDeviation*1000);
            }
            return log;
        }
    };
}


void
HypoDD::evalXCorr()
{
    bool theoretical = false; // this is useful for testing the ability of detecting phases

    XCorrEvalStats totalStats;
    map<string,XCorrEvalStats> statsByStation; // key station id
    map<int,XCorrEvalStats> statsByInterEvDistance; // key distance
    map<int,XCorrEvalStats> statsByStaDistance; // key distance
    const double EV_DIST_STEP = 0.1; // km
    const double STA_DIST_STEP = 3; // km

    auto printStats = [&](string title, bool theoretical)
    {
        string log = title + "\n";
        log += stringify("Cumulative stats: %s\n", totalStats.describe(theoretical).c_str());

        log += stringify("Stats by inter-event distance in %.2f km step\n", EV_DIST_STEP);
        for ( const auto& kv : statsByInterEvDistance)
        {
            XCorrEvalStats tmp = kv.second;
            tmp.normalize();
            log += stringify("Inter-event dist %.2f-%-.2f [km]: ",
                              kv.first*EV_DIST_STEP, (kv.first+1)*EV_DIST_STEP);
            log += stringify("#CC %6d good CC %3.f%% avg coeff %.2f avg #matches/ev %2d\n",
                              tmp.total, (tmp.goodCC * 100. / tmp.total),
                              tmp.meanCoeff, tmp.meanCount);
            if ( theoretical )
                log += stringify(" time-diff %3.f [msec] abs time-diff %3.f [msec]",
                                 tmp.deviation*1000, tmp.absDeviation*1000);
        }

        log += stringify("Stats by event to station distance in %.2f km step\n", STA_DIST_STEP);
        for ( const auto& kv : statsByStaDistance)
        {
            log += stringify("Station dist %3d-%-3d [km]: %s\n", int(kv.first*STA_DIST_STEP),
                            int((kv.first+1)*STA_DIST_STEP),
                            kv.second.describe(theoretical).c_str());
        }

        log += stringify("Stats by station\n");
        for ( const auto& kv : statsByStation)
        {
            log += stringify("%-12s: %s\n", kv.first.c_str(),
                             kv.second.describe(theoretical).c_str());
        }
        SEISCOMP_WARNING("%s", log.c_str() );
    };

    _counters = {0};
    _wf->resetCounters(); 
    int loop = 0;

    for (const auto& kv : _ddbgc->getEvents() )
    {
        const Event& event = kv.second;

        // find the neighbouring events
        NeighboursPtr neighbours;
        try {
            neighbours = selectNeighbouringEvents(
                _ddbgc, event, _ddbgc, _cfg.ddObservations2.minWeight,
                _cfg.ddObservations2.minESdist, _cfg.ddObservations2.maxESdist,
                _cfg.ddObservations2.minEStoIEratio, _cfg.ddObservations2.minDTperEvt,
                _cfg.ddObservations2.maxDTperEvt, _cfg.ddObservations2.minNumNeigh,
                _cfg.ddObservations2.maxNumNeigh, _cfg.ddObservations2.numEllipsoids,
                _cfg.ddObservations2.maxEllipsoidSize, false);
        } catch ( ... ) { continue; }

        CatalogPtr catalog;

        if ( theoretical )
        {
            catalog = neighbours->fromNeighbours(_ddbgc, false);
 
            // create theoretical phases for this event
            // beware: no event phases are present in catalog so the event
            // will end up with only theoretical phases
            addMissingEventPhases(event, catalog, _ddbgc, neighbours);
        }
        else
        {
            catalog = neighbours->fromNeighbours(_ddbgc, true);
        }

        // cross correlate every neighbour phase with corresponding event theoretical phase
        XCorrCache xcorr;
        buildXcorrDiffTTimePairs(catalog, neighbours, event, xcorr);

        // Update theoretical and automatic phase pick time and uncertainties based on
        // cross-correlation results
        // Also drop theoretical phases wihout any good cross correlation result
        if ( theoretical )
            fixPhases(catalog, event, xcorr);

        //
        // Compare the detected phases with the actual event phases (manual or automatic)
        //
        XCorrEvalStats evStats;
        map<unsigned,XCorrEvalStats> statsByNeighbour; // key neighbour id

        for ( const auto& kv : neighbours->allPhases() )
            for ( Phase::Type phaseType : kv.second )
        {
            const string stationId = kv.first;
            const Phase& catalogPhase = _ddbgc->searchPhase(event.id, stationId, phaseType)->second;

            XCorrEvalStats phStaStats;
            phStaStats.total = 1;

            if ( xcorr.has(event.id, stationId, phaseType) )
            {
                const Phase& detectedPhase = catalog->searchPhase(event.id, stationId,
                                                                  phaseType)->second;
                phStaStats.goodCC = 1;
                double deviation = (catalogPhase.time - detectedPhase.time).length();
                phStaStats.deviation = deviation;
                phStaStats.absDeviation = std::abs(deviation);
                auto& pdata = xcorr.get(event.id, stationId, phaseType);
                phStaStats.meanCoeff = pdata.mean_coeff;
                phStaStats.meanCount = pdata.ccCount;
            }

            evStats += phStaStats;
            statsByStation[catalogPhase.stationId] += phStaStats;

            const Station& station = _ddbgc->getStations().at(catalogPhase.stationId);
            double stationDistance = computeDistance(event, station);
            statsByStaDistance[ int(stationDistance/STA_DIST_STEP) ] += phStaStats;

            //
            //  collect stats by neighbour
            //
            for ( unsigned neighEvId : neighbours->ids )
            {
                if ( neighbours->has(neighEvId, stationId, phaseType) )
                {
                    XCorrEvalStats& neighStats = statsByNeighbour[neighEvId];
                    neighStats.total++;
                    if ( xcorr.has(event.id, neighEvId, stationId, phaseType) )
                    {
                        const auto& pdata = xcorr.get(event.id, neighEvId, stationId, phaseType);
                        neighStats.goodCC++;
                        neighStats.meanCount++;
                        neighStats.meanCoeff += pdata.coeff;
                        double deviation = phStaStats.deviation;
                        deviation -= xcorr.get(event.id, stationId, phaseType).mean_lag - pdata.lag;
                        neighStats.deviation += deviation;
                        neighStats.absDeviation += std::abs(deviation);
                    }
                }
            }
        }

        //
        //  collect stats by inter event distance
        //
        for ( auto& kv : statsByNeighbour )
        {
            const Event& neighbEv = catalog->getEvents().at(kv.first);
            XCorrEvalStats& neighStats = kv.second;
            neighStats.meanCount *= neighStats.meanCount;
            double interEvDistance = computeDistance(event, neighbEv); 
            statsByInterEvDistance[ int(interEvDistance/EV_DIST_STEP) ] += neighStats;
        }

        // total stats
        totalStats += evStats;
        SEISCOMP_WARNING("Event %-5s mag %3.1f %s", string(event).c_str(), event.magnitude,
                         evStats.describe(theoretical).c_str());

        if ( ++loop % 50 == 0 )
        {
            printStats("<<<Progressive stats>>>", theoretical);
        }
    }

    printStats("<<<Final stats>>>", theoretical);
    printCounters();

}


} // HDD
} // Seiscomp


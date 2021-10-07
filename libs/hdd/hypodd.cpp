/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as          *
 * published by the Free Software Foundation, either version 3 of the      *
 * License, or (at your option) any later version.                         *
 *                                                                         *
 * This software is distributed in the hope that it will be useful,        *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#include "hypodd.h"
#include "sccatalog.h"
#include "utils.h"

#include <boost/filesystem.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <seiscomp3/client/inventory.h>
#include <seiscomp3/core/datetime.h>
#include <seiscomp3/core/strings.h>
#include <seiscomp3/core/typedarray.h>
#include <seiscomp3/io/recordinput.h>
#include <seiscomp3/utils/files.h>
#include <stdexcept>

#define SEISCOMP_COMPONENT HDD
#include <seiscomp3/logging/file.h>
#include <seiscomp3/logging/log.h>

using namespace std;
using namespace Seiscomp;
using Seiscomp::Core::stringify;
using Event   = HDD::Catalog::Event;
using Phase   = HDD::Catalog::Phase;
using Station = HDD::Catalog::Station;
using HDD::Waveform::getBandAndInstrumentCodes;

namespace Seiscomp {
namespace HDD {

HypoDD::HypoDD(const CatalogCPtr &catalog,
               const Config &cfg,
               const string &workingDir)
    : _workingDir(workingDir), _cfg(cfg)
{
  setCatalog(catalog);

  if (!Util::pathExists(_workingDir))
  {
    if (!Util::createPath(_workingDir))
    {
      string msg = "Unable to create working directory: " + _workingDir;
      throw runtime_error(msg);
    }
  }

  _cacheDir = (boost::filesystem::path(_workingDir) / "wfcache").string();
  if (!Util::pathExists(_cacheDir))
  {
    if (!Util::createPath(_cacheDir))
    {
      string msg = "Unable to create cache directory: " + _cacheDir;
      throw runtime_error(msg);
    }
  }

  _tmpCacheDir = (boost::filesystem::path(_workingDir) / "tmpcache").string();
  if (!Util::pathExists(_tmpCacheDir))
  {
    if (!Util::createPath(_tmpCacheDir))
    {
      string msg = "Unable to create cache directory: " + _tmpCacheDir;
      throw runtime_error(msg);
    }
  }

  _wfDebugDir = (boost::filesystem::path(_workingDir) / "wfdebug").string();

  setUseCatalogWaveformDiskCache(true);
  setWaveformCacheAll(false);
  setWaveformDebug(false);
}

void HypoDD::setCatalog(const CatalogCPtr &catalog)
{
  _srcCat = catalog;
  _bgCat  = Catalog::filterPhasesAndSetWeights(
      *_srcCat, Phase::Source::CATALOG, _cfg.validPphases, _cfg.validSphases);
}

void HypoDD::setUseCatalogWaveformDiskCache(bool cache)
{
  _useCatalogWaveformDiskCache = cache;
  createWaveformCache();
}

void HypoDD::createWaveformCache()
{
  _wfAccess.unloadableWfs.clear();
  _wfAccess.loader    = nullptr;
  _wfAccess.diskCache = nullptr;
  _wfAccess.extraLen  = nullptr;
  _wfAccess.snrFilter = nullptr;
  _wfAccess.memCache  = nullptr;

  _wfAccess.loader = new Waveform::Loader(_cfg.recordStreamURL);

  if (_useCatalogWaveformDiskCache)
  {
    _wfAccess.diskCache =
        new Waveform::DiskCachedLoader(_wfAccess.loader, _cacheDir);
    _wfAccess.extraLen =
        new Waveform::ExtraLenLoader(_wfAccess.diskCache, DISK_TRACE_MIN_LEN);

    if (_cfg.snr.minSnr > 0)
    {
      _wfAccess.snrFilter = new Waveform::SnrFilteredLoader(
          _wfAccess.extraLen, _cfg.snr.minSnr, _cfg.snr.noiseStart,
          _cfg.snr.noiseEnd, _cfg.snr.signalStart, _cfg.snr.signalEnd);
      _wfAccess.memCache = new Waveform::MemCachedLoader(_wfAccess.snrFilter);
    }
    else
    {
      _wfAccess.memCache = new Waveform::MemCachedLoader(_wfAccess.extraLen);
    }
  }
  else
  {
    if (_cfg.snr.minSnr > 0)
    {
      _wfAccess.snrFilter = new Waveform::SnrFilteredLoader(
          _wfAccess.loader, _cfg.snr.minSnr, _cfg.snr.noiseStart,
          _cfg.snr.noiseEnd, _cfg.snr.signalStart, _cfg.snr.signalEnd);
      _wfAccess.memCache = new Waveform::MemCachedLoader(_wfAccess.snrFilter);
    }
    else
    {
      _wfAccess.memCache = new Waveform::MemCachedLoader(_wfAccess.loader);
    }
  }
}

void HypoDD::setWaveformDebug(bool debug)
{
  _waveformDebug = debug;
  if (_waveformDebug)
  {
    if (!Util::pathExists(_wfDebugDir))
    {
      if (!Util::createPath(_wfDebugDir))
      {
        string msg =
            "Unable to create waveform debug directory: " + _wfDebugDir;
        throw runtime_error(msg);
      }
    }
  }

  if (_wfAccess.snrFilter)
    _wfAccess.snrFilter->setDebugDirectory(_waveformDebug ? _wfDebugDir : "");
  _wfAccess.memCache->setDebugDirectory(_waveformDebug ? _wfDebugDir : "");
}

string HypoDD::generateWorkingSubDir(const Event &ev) const
{
  static UniformRandomer ran(0, 9999);
  string id = stringify(
      "%s_%05d_%06d_%s_%04zu",
      ev.time.toString("%Y%m%d%H%M%S").c_str(),           // origin time
      int(ev.latitude * 1000),                            // Latitude
      int(ev.longitude * 1000),                           // Longitude
      Core::Time::GMT().toString("%Y%m%d%H%M%S").c_str(), // creation time
      ran.next()                                          // random number
  );
  return id;
}

void HypoDD::preloadWaveforms()
{
  SEISCOMP_INFO("Preloading catalog waveform data (%lu events to load)",
                _bgCat->getEvents().size());

  resetCounters();

  unsigned numPhases = 0, numSPhases = 0, numEvents = 0;

  //
  // preload waveforms on disk and cache them in memory (pre-processed)
  //
  for (const auto &kv : _bgCat->getEvents())
  {
    const Event &event = kv.second;
    auto eqlrng        = _bgCat->getPhases().equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      const Phase &phase  = it->second;
      Core::TimeWindow tw = xcorrTimeWindowLong(phase);
      const auto xcorrCfg = _cfg.xcorr.at(phase.procInfo.type);

      for (string component : xcorrCfg.components)
      {
        Phase tmpPh = phase;
        tmpPh.channelCode =
            getBandAndInstrumentCodes(tmpPh.channelCode) + component;
        getWaveform(tw, event, tmpPh, _wfAccess.memCache);
      }

      numPhases++;
      if (phase.procInfo.type == Phase::Type::S) numSPhases++;
    }

    if (++numEvents % (_bgCat->getEvents().size() / 100) == 0)
    {
      SEISCOMP_INFO("Loaded %lu%% of waveforms",
                    (numEvents * 100 / _bgCat->getEvents().size()));
    }
  }

  updateCounters();
  SEISCOMP_INFO(
      "Finished preloading catalog waveform data: total phases %u (P %.f%%, S "
      "%.f%%) phases with Signal to Noise ratio too low %u (%.f%%), "
      "phases data not available %u (%.f%%), "
      "(waveforms downloaded %u, waveforms loaded from disk cache %u)",
      numPhases, ((numPhases - numSPhases) * 100. / numPhases),
      (numSPhases * 100. / numPhases), _counters.wf_snr_low,
      (_counters.wf_snr_low * 100. / numPhases), _counters.wf_no_avail,
      (_counters.wf_no_avail * 100. / numPhases), _counters.wf_downloaded,
      _counters.wf_disk_cached);
}

CatalogPtr HypoDD::relocateMultiEvents(const ClusteringOptions &clustOpt,
                                       const SolverOptions &solverOpt)
{
  SEISCOMP_INFO("Starting HypoDD relocator in multiple events mode");

  if (!_ttt) _ttt = TravelTimeTable::create(_cfg.ttt.type, _cfg.ttt.model);

  CatalogPtr catToReloc(new Catalog(*_bgCat));

  // prepare a folder for debug files
  string catalogWorkingDir =
      (boost::filesystem::path(_workingDir) / "catalog").string();
  if (_saveProcessing)
  {
    if (!Util::pathExists(catalogWorkingDir))
    {
      if (!Util::createPath(catalogWorkingDir))
      {
        string msg = "Unable to create working directory: " + catalogWorkingDir;
        throw runtime_error(msg);
      }
    }
    SEISCOMP_INFO("Working dir %s", catalogWorkingDir.c_str());
  }

  // prepare file logger
  std::shared_ptr<Logging::Output> processingInfoOutput;
  if (_saveProcessing)
  {
    processingInfoOutput =
        std::shared_ptr<Logging::Output>(new Logging::FileOutput(
            (boost::filesystem::path(catalogWorkingDir) / "info.log")
                .string()
                .c_str()));
    processingInfoOutput->subscribe(Seiscomp::Logging::_SCInfoChannel);
    processingInfoOutput->subscribe(Seiscomp::Logging::_SCWarningChannel);
    processingInfoOutput->subscribe(Seiscomp::Logging::_SCErrorChannel);
  }

  // find Neighbours for each event in the catalog
  list<NeighboursPtr> allNeighbours = selectNeighbouringEventsCatalog(
      catToReloc, clustOpt.minWeight, clustOpt.minESdist, clustOpt.maxESdist,
      clustOpt.minEStoIEratio, clustOpt.minDTperEvt, clustOpt.maxDTperEvt,
      clustOpt.minNumNeigh, clustOpt.maxNumNeigh, clustOpt.numEllipsoids,
      clustOpt.maxEllipsoidSize, true);

  // Organize the neighbours by not connected clusters. In addition,
  // don't report the same pair multiple times (e.g. ev1-ev2 and ev2-ev1)
  // since we only need one observation per pair in the DD solver.
  deque<list<NeighboursPtr>> clusters =
      clusterizeNeighbouringEvents(allNeighbours);

  SEISCOMP_INFO("Found %lu event clusters", clusters.size());

  if (_saveProcessing)
  {
    catToReloc->writeToFile(
        (boost::filesystem::path(catalogWorkingDir) / "starting-event.csv")
            .string(),
        (boost::filesystem::path(catalogWorkingDir) / "starting-phase.csv")
            .string(),
        (boost::filesystem::path(catalogWorkingDir) / "starting-station.csv")
            .string());
  }

  //
  // relocate one cluster a time
  //
  CatalogPtr relocatedCatalog(new Catalog());

  for (unsigned clusterId = 0; clusterId < clusters.size(); clusterId++)
  {
    const list<NeighboursPtr> &neighCluster = clusters[clusterId];

    SEISCOMP_INFO("Relocating cluster %u (%lu events)", clusterId + 1,
                  neighCluster.size());

    if (_saveProcessing)
    {
      CatalogPtr catToDump(new Catalog());
      for (const NeighboursPtr &n : neighCluster)
        catToDump->add(n->refEvId, *catToReloc, true);
      string prefix = "cluster-" + to_string(clusterId + 1);
      catToDump->writeToFile(
          (boost::filesystem::path(catalogWorkingDir) / (prefix + "-event.csv"))
              .string(),
          (boost::filesystem::path(catalogWorkingDir) / (prefix + "-phase.csv"))
              .string(),
          (boost::filesystem::path(catalogWorkingDir) /
           (prefix + "-station.csv"))
              .string());
    }

    // Perform cross-correlation which also detects picks around theoretical
    // arrival times. The catalog will be updated with the corresponding
    // theoretical phases.
    const XCorrCache xcorr = buildXCorrCache(
        catToReloc, neighCluster, _useArtificialPhases,
        clustOpt.xcorrMaxEvStaDist, clustOpt.xcorrMaxInterEvDist);

    // the actual relocation
    CatalogPtr relocatedCluster =
        relocate(catToReloc, neighCluster, solverOpt, false, xcorr);

    relocatedCatalog->add(*relocatedCluster, true);

    if (_saveProcessing)
    {
      string prefix = "relocated-cluster-" + to_string(clusterId + 1);
      relocatedCluster->writeToFile(
          (boost::filesystem::path(catalogWorkingDir) / (prefix + "-event.csv"))
              .string(),
          (boost::filesystem::path(catalogWorkingDir) / (prefix + "-phase.csv"))
              .string(),
          (boost::filesystem::path(catalogWorkingDir) /
           (prefix + "-station.csv"))
              .string());
    }
  }

  if (_saveProcessing)
  {
    relocatedCatalog->writeToFile(
        (boost::filesystem::path(catalogWorkingDir) / "relocated-event.csv")
            .string(),
        (boost::filesystem::path(catalogWorkingDir) / "relocated-phase.csv")
            .string(),
        (boost::filesystem::path(catalogWorkingDir) / "relocated-station.csv")
            .string());
  }

  return relocatedCatalog;
}

CatalogPtr HypoDD::relocateSingleEvent(const CatalogCPtr &singleEvent,
                                       const ClusteringOptions &clustOpt1,
                                       const ClusteringOptions &clustOpt2,
                                       const SolverOptions &solverOpt)
{
  if (!_ttt) _ttt = TravelTimeTable::create(_cfg.ttt.type, _cfg.ttt.model);

  const CatalogCPtr bgCat = _bgCat;

  // there must be only one event in the catalog, the origin to relocate
  const Event &evToRelocate = singleEvent->getEvents().begin()->second;
  const auto &evToRelocatePhases =
      singleEvent->getPhases().equal_range(evToRelocate.id);

  SEISCOMP_INFO(
      "Starting HypoDD relocator in single event mode: event %s lat %.6f lon "
      "%.6f depth %.4f mag %.2f time %s #phases %ld",
      string(evToRelocate).c_str(), evToRelocate.latitude,
      evToRelocate.longitude, evToRelocate.depth, evToRelocate.magnitude,
      evToRelocate.time.iso().c_str(),
      std::distance(evToRelocatePhases.first, evToRelocatePhases.second));

  // prepare a folder for debug files
  string baseWorkingDir;
  if (_saveProcessing)
  {
    do
    {
      baseWorkingDir = generateWorkingSubDir(evToRelocate);
      baseWorkingDir =
          (boost::filesystem::path(_workingDir) / baseWorkingDir).string();
    } while (Util::pathExists(baseWorkingDir));

    if (!Util::createPath(baseWorkingDir))
    {
      string msg = "Unable to create working directory: " + baseWorkingDir;
      throw runtime_error(msg);
    }
    SEISCOMP_INFO("Working dir %s", baseWorkingDir.c_str());
  }

  // prepare file logger
  std::shared_ptr<Logging::Output> processingInfoOutput;
  if (_saveProcessing)
  {
    processingInfoOutput =
        std::shared_ptr<Logging::Output>(new Logging::FileOutput(
            (boost::filesystem::path(baseWorkingDir) / "info.log")
                .string()
                .c_str()));
    processingInfoOutput->subscribe(Seiscomp::Logging::_SCInfoChannel);
    processingInfoOutput->subscribe(Seiscomp::Logging::_SCWarningChannel);
    processingInfoOutput->subscribe(Seiscomp::Logging::_SCErrorChannel);
  }

  SEISCOMP_INFO(
      "Performing step 1: initial location refinement (no cross-correlation)");

  string eventWorkingDir =
      (boost::filesystem::path(baseWorkingDir) / "step1").string();

  CatalogPtr evToRelocateCat =
      Catalog::filterPhasesAndSetWeights(*singleEvent, Phase::Source::RT_EVENT,
                                         _cfg.validPphases, _cfg.validSphases);

  CatalogPtr relocatedEvCat =
      relocateEventSingleStep(bgCat, evToRelocateCat, eventWorkingDir,
                              clustOpt1, solverOpt, false, false);

  if (relocatedEvCat)
  {
    const Event &ev = relocatedEvCat->getEvents().begin()->second;
    SEISCOMP_INFO("Step 1 relocation successful, new location: "
                  "lat %.6f lon %.6f depth %.4f time %s",
                  ev.latitude, ev.longitude, ev.depth, ev.time.iso().c_str());
    SEISCOMP_INFO("Relocation report: %s",
                  relocationReport(relocatedEvCat).c_str());

    evToRelocateCat = relocatedEvCat;
  }
  else
  {
    SEISCOMP_ERROR("Failed to perform step 1 origin relocation");
  }

  SEISCOMP_INFO("Performing step 2: relocation with cross-correlation");

  eventWorkingDir =
      (boost::filesystem::path(baseWorkingDir) / "step2").string();

  CatalogPtr relocatedEvWithXcorr =
      relocateEventSingleStep(bgCat, evToRelocateCat, eventWorkingDir,
                              clustOpt2, solverOpt, true, _useArtificialPhases);

  if (relocatedEvWithXcorr)
  {
    Event ev = relocatedEvWithXcorr->getEvents().begin()->second;
    SEISCOMP_INFO("Step 2 relocation successful, new location: "
                  "lat %.6f lon %.6f depth %.4f time %s",
                  ev.latitude, ev.longitude, ev.depth, ev.time.iso().c_str());
    SEISCOMP_INFO("Relocation report: %s",
                  relocationReport(relocatedEvWithXcorr).c_str());

    // update the "origin change information" taking into consideration
    // the first relocation step, too
    if (relocatedEvCat)
    {
      const Event &prevRelocEv = relocatedEvCat->getEvents().begin()->second;
      if (prevRelocEv.relocInfo.isRelocated)
      {
        ev.relocInfo.locChange += prevRelocEv.relocInfo.locChange;
        ev.relocInfo.depthChange += prevRelocEv.relocInfo.depthChange;
        ev.relocInfo.timeChange += prevRelocEv.relocInfo.timeChange;
        ev.relocInfo.startRms = prevRelocEv.relocInfo.startRms;
        relocatedEvWithXcorr->updateEvent(ev);
      }
    }

    SEISCOMP_INFO("Total Changes: location=%.2f[km] depth=%.2f[km] "
                  "time=%.3f[sec] Rms=%.3f[sec] (before/after %.3f/%.3f)",
                  ev.relocInfo.locChange, ev.relocInfo.depthChange,
                  ev.relocInfo.timeChange, (ev.rms - ev.relocInfo.startRms),
                  ev.relocInfo.startRms, ev.rms);
  }
  else
  {
    SEISCOMP_ERROR("Failed to perform step 2 origin relocation");
  }

  if (!relocatedEvWithXcorr) throw runtime_error("Failed origin relocation");

  return relocatedEvWithXcorr;
}

CatalogPtr HypoDD::relocateEventSingleStep(const CatalogCPtr bgCat,
                                           const CatalogCPtr &evToRelocateCat,
                                           const string &workingDir,
                                           const ClusteringOptions &clustOpt,
                                           const SolverOptions &solverOpt,
                                           bool doXcorr,
                                           bool computeTheoreticalPhases)
{
  if (_saveProcessing)
  {
    if (!Util::createPath(workingDir))
    {
      string msg = "Unable to create working directory: " + workingDir;
      throw runtime_error(msg);
    }
    SEISCOMP_INFO("Working dir %s", workingDir.c_str());

    evToRelocateCat->writeToFile(
        (boost::filesystem::path(workingDir) / "single-event.csv").string(),
        (boost::filesystem::path(workingDir) / "single-event-phase.csv")
            .string(),
        (boost::filesystem::path(workingDir) / "single-event-station.csv")
            .string());
  }

  CatalogPtr relocatedEvCat;

  try
  {
    // extract event to relocate
    const Event &evToRelocate = evToRelocateCat->getEvents().begin()->second;

    //
    // select neighbouring events
    //
    bool keepUnmatchedPhases = doXcorr; // useful for detecting missed picks

    NeighboursPtr neighbours = selectNeighbouringEvents(
        bgCat, evToRelocate, evToRelocateCat, clustOpt.minWeight,
        clustOpt.minESdist, clustOpt.maxESdist, clustOpt.minEStoIEratio,
        clustOpt.minDTperEvt, clustOpt.maxDTperEvt, clustOpt.minNumNeigh,
        clustOpt.maxNumNeigh, clustOpt.numEllipsoids, clustOpt.maxEllipsoidSize,
        keepUnmatchedPhases);

    //
    // prepare catalog to relocate
    //
    CatalogPtr catalog = neighbours->toCatalog(bgCat);
    unsigned evToRelocateNewId =
        catalog->add(evToRelocate.id, *evToRelocateCat, false);
    neighbours->refEvId = evToRelocateNewId;

    if (_saveProcessing)
    {
      catalog->writeToFile(
          (boost::filesystem::path(workingDir) / "starting-event.csv").string(),
          (boost::filesystem::path(workingDir) / "starting-phase.csv").string(),
          (boost::filesystem::path(workingDir) / "starting-station.csv")
              .string());
    }

    XCorrCache xcorr;
    if (doXcorr)
    {
      // Perform cross-correlation, which also detects picks around theoretical
      // arrival times. The catalog will be updated with the corresponding
      // phases.
      xcorr = buildXCorrCache(catalog, {neighbours}, computeTheoreticalPhases,
                              clustOpt.xcorrMaxEvStaDist,
                              clustOpt.xcorrMaxInterEvDist);
    }

    // the actual relocation
    relocatedEvCat = relocate(catalog, {neighbours}, solverOpt, true, xcorr);

    if (_saveProcessing)
    {
      relocatedEvCat->writeToFile(
          (boost::filesystem::path(workingDir) / "relocated-event.csv")
              .string(),
          (boost::filesystem::path(workingDir) / "relocated-phase.csv")
              .string(),
          (boost::filesystem::path(workingDir) / "relocated-station.csv")
              .string());
    }
  }
  catch (exception &e)
  {
    SEISCOMP_ERROR("%s", e.what());
  }

  return relocatedEvCat;
}

CatalogPtr HypoDD::relocate(const CatalogCPtr &catalog,
                            const std::list<NeighboursPtr> &neighCluster,
                            const SolverOptions &solverOpt,
                            bool keepNeighboursFixed,
                            const XCorrCache &xcorr) const
{
  SEISCOMP_INFO("Building and solving double-difference system...");

  //
  // iterate the solver computation multiple times
  //
  CatalogCPtr finalCatalog = catalog;
  unordered_map<unsigned, NeighboursPtr> finalNeighCluster;
  ObservationParams obsparams;
  for (unsigned iteration = 0; iteration < solverOpt.algoIterations;
       iteration++)
  {
    //
    // compute parameters for this loop iteration
    //
    auto interpolate = [&](double start, double end) -> double {
      if (solverOpt.algoIterations < 2) return (start + end) / 2;
      return start + (end - start) * iteration / (solverOpt.algoIterations - 1);
    };

    double dampingFactor =
        interpolate(solverOpt.dampingFactorStart, solverOpt.dampingFactorEnd);
    double downWeightingByResidual =
        interpolate(solverOpt.downWeightingByResidualStart,
                    solverOpt.downWeightingByResidualEnd);
    double absLocConstraint   = interpolate(solverOpt.absLocConstraintStart,
                                          solverOpt.absLocConstraintEnd);
    double absTTDiffObsWeight = interpolate(1.0, solverOpt.absTTDiffObsWeight);
    double xcorrObsWeight     = interpolate(1.0, solverOpt.xcorrObsWeight);

    SEISCOMP_INFO("Solving iteration %u num events %lu. Parameters: "
                  "observWeight TT/CC=%.2f/%.2f dampingFactor=%.2f "
                  "downWeightingByResidual=%.2f absLocConstraint=%.2f",
                  iteration, neighCluster.size(), absTTDiffObsWeight,
                  xcorrObsWeight, dampingFactor, downWeightingByResidual,
                  absLocConstraint);

    // create a solver and then add observations
    Solver solver(solverOpt.type);

    //
    // Add absolute travel time/cross-correlation differences to the solver
    // (the observations)
    //
    for (const NeighboursPtr &neighbours : neighCluster)
    {
      addObservations(solver, absTTDiffObsWeight, xcorrObsWeight, finalCatalog,
                      neighbours, keepNeighboursFixed,
                      solverOpt.usePickUncertainty, xcorr, obsparams);
    }
    obsparams.addToSolver(solver);

    //
    // solve the system
    //
    try
    {
      solver.solve(solverOpt.solverIterations, absLocConstraint, dampingFactor,
                   downWeightingByResidual, solverOpt.L2normalization);
    }
    catch (exception &e)
    {
      SEISCOMP_INFO("Cannot solve the double-difference system, stop here (%s)",
                    e.what());
      break;
    }

    // prepare for next iteration
    obsparams = ObservationParams();

    // update event parameters
    finalCatalog = updateRelocatedEvents(
        solver, finalCatalog, neighCluster, obsparams,
        std::max(absTTDiffObsWeight, xcorrObsWeight), finalNeighCluster);
  }

  // compute last bit of statistics for the relocated events
  return updateRelocatedEventsFinalStats(catalog, finalCatalog,
                                         finalNeighCluster);
}

string HypoDD::relocationReport(const CatalogCPtr &relocatedEv)
{
  const Event &event = relocatedEv->getEvents().begin()->second;
  if (!event.relocInfo.isRelocated) return "Event not relocated";

  return stringify(
      "Origin changes: location=%.2f[km] depth=%.2f[km] time=%.3f[sec] "
      "Rms change [sec]: %.3f (before/after %.3f/%.3f) "
      "Neighbours=%u Used Phases: P=%u S=%u "
      "Stations distance [km]: min=%.1f median=%.1f max=%.1f "
      "Neighbours mean distace to centroid [km]: location=%.2f depth=%.2f "
      "Origin distace to neighbours centroid [km]: location=%.2f depth=%.2f "
      "DD observations: %u (CC P/S %u/%u TT P/S %u/%u) "
      "DD observations residuals [msec]: before=%.f+/-%.1f after=%.f+/-%.1f",
      event.relocInfo.locChange, event.relocInfo.depthChange,
      event.relocInfo.timeChange, (event.rms - event.relocInfo.startRms),
      event.relocInfo.startRms, event.rms, event.relocInfo.neighbours.amount,
      event.relocInfo.phases.usedP, event.relocInfo.phases.usedS,
      event.relocInfo.phases.stationDistMin,
      event.relocInfo.phases.stationDistMedian,
      event.relocInfo.phases.stationDistMax,
      event.relocInfo.neighbours.meanDistToCentroid,
      event.relocInfo.neighbours.meanDepthDistToCentroid,
      event.relocInfo.neighbours.eventDistToCentroid,
      event.relocInfo.neighbours.eventDepthDistToCentroid,
      (event.relocInfo.ddObs.numCCp + event.relocInfo.ddObs.numCCs +
       event.relocInfo.ddObs.numTTp + event.relocInfo.ddObs.numTTs),
      event.relocInfo.ddObs.numCCp, event.relocInfo.ddObs.numCCs,
      event.relocInfo.ddObs.numTTp, event.relocInfo.ddObs.numTTs,
      event.relocInfo.ddObs.startResidualMedian * 1000,
      event.relocInfo.ddObs.startResidualMAD * 1000,
      event.relocInfo.ddObs.finalResidualMedian * 1000,
      event.relocInfo.ddObs.finalResidualMAD * 1000);
}

/*
 * Add both the absolute travel time differences and the differential
 * travel times from the cross-correlation for pairs of earthquakes to the
 * solver.
 */
void HypoDD::addObservations(Solver &solver,
                             double absTTDiffObsWeight,
                             double xcorrObsWeight,
                             const CatalogCPtr &catalog,
                             const NeighboursPtr &neighbours,
                             bool keepNeighboursFixed,
                             bool usePickUncertainty,
                             const XCorrCache &xcorr,
                             ObservationParams &obsparams) const
{
  // copy event because we'll update it
  const Event &refEv = catalog->getEvents().at(neighbours->refEvId);

  //
  // loop through reference event phases
  //
  auto eqlrng = catalog->getPhases().equal_range(refEv.id);
  for (auto it = eqlrng.first; it != eqlrng.second; ++it)
  {
    const Phase &refPhase  = it->second;
    const Station &station = catalog->getStations().at(refPhase.stationId);
    char phaseTypeAsChar   = static_cast<char>(refPhase.procInfo.type);

    //
    // loop through neighbouring events and look for the matching phase
    //
    for (unsigned neighEvId : neighbours->ids)
    {
      const Event &event = catalog->getEvents().at(neighEvId);

      if (!neighbours->has(neighEvId, refPhase.stationId,
                           refPhase.procInfo.type))
        continue;

      const Phase &phase = catalog
                               ->searchPhase(event.id, refPhase.stationId,
                                             refPhase.procInfo.type)
                               ->second;
      //
      // compute travel times for both event and `refEv`
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

      if (!obsparams.add(_ttt, refEv, station, refPhase, true) ||
          !obsparams.add(_ttt, event, station, phase, !keepNeighboursFixed))
      {
        SEISCOMP_DEBUG("Skipping observation (ev %u-%u sta %s phase %c)",
                       refEv.id, event.id, station.id.c_str(), phaseTypeAsChar);
        continue;
      }

      //
      // compute absolute travel time differences to the solver
      //
      double diffTime = ref_travel_time - travel_time;
      double weight =
          usePickUncertainty
              ? (refPhase.procInfo.weight + phase.procInfo.weight) / 2.0
              : 1.0;
      bool isXcorr = false;

      //
      // Check if we have cross-correlation results for the current
      // event/`refEv` pair at station/phase and use those instead.
      //
      if (xcorr.has(refEv.id, event.id, refPhase.stationId,
                    refPhase.procInfo.type))
      {

        const auto &xcdata = xcorr.get(refEv.id, event.id, refPhase.stationId,
                                       refPhase.procInfo.type);
        diffTime -= xcdata.lag;
        weight *= xcorrObsWeight;
        isXcorr = true;
      }
      else
      {
        weight *= absTTDiffObsWeight;
      }

      solver.addObservation(refEv.id, event.id, refPhase.stationId,
                            phaseTypeAsChar, diffTime, weight, isXcorr);
    }
  }
}

bool HypoDD::ObservationParams::add(HDD::TravelTimeTablePtr ttt,
                                    const Event &event,
                                    const Station &station,
                                    const Phase &phase,
                                    bool computeEvChanges)
{
  char phaseType = static_cast<char>(phase.procInfo.type);
  const std::string key =
      std::to_string(event.id) + "@" + station.id + ":" + phaseType;
  if (_entries.find(key) == _entries.end())
  {
    try
    {
      double travelTime, takeOfAngleAzim, takeOfAngleDip, velocityAtSrc;
      ttt->compute(event, station, string(1, phaseType), travelTime,
                   takeOfAngleAzim, takeOfAngleDip, velocityAtSrc);
      double ttResidual = travelTime - (phase.time - event.time).length();
      _entries[key]     = Entry{event,          station,       phaseType,
                            travelTime,     ttResidual,    takeOfAngleAzim,
                            takeOfAngleDip, velocityAtSrc, computeEvChanges};
    }
    catch (exception &e)
    {
      SEISCOMP_WARNING(
          "Travel Time Table error: %s (Event lat %.6f lon %.6f depth %.6f "
          "Station lat %.6f lon %.6f elevation %.f )",
          e.what(), event.latitude, event.longitude, event.depth,
          station.latitude, station.longitude, station.elevation);
      return false;
    }
  }
  return true;
}

const HypoDD::ObservationParams::Entry &HypoDD::ObservationParams::get(
    unsigned eventId, const std::string stationId, char phaseType) const
{
  const std::string key =
      std::to_string(eventId) + "@" + stationId + ":" + phaseType;
  return _entries.at(key);
}

void HypoDD::ObservationParams::addToSolver(Solver &solver) const
{
  for (const auto &kv : _entries)
  {
    const ObservationParams::Entry &e = kv.second;
    solver.addObservationParams(
        e.event.id, e.station.id, e.phaseType, e.event.latitude,
        e.event.longitude, e.event.depth, e.station.latitude,
        e.station.longitude, e.station.elevation, e.computeEvChanges,
        e.travelTime, e.travelTimeResidual, e.takeOfAngleAzim, e.takeOfAngleDip,
        e.velocityAtSrc);
  }
}

CatalogPtr HypoDD::updateRelocatedEvents(
    const Solver &solver,
    const CatalogCPtr &catalog,
    const std::list<NeighboursPtr> &neighCluster,
    ObservationParams &obsparams,
    double pickWeightScaler,
    std::unordered_map<unsigned, NeighboursPtr> &finalNeighCluster // output
) const
{
  unordered_map<string, Station> stations    = catalog->getStations();
  map<unsigned, Event> events                = catalog->getEvents();
  unordered_multimap<unsigned, Phase> phases = catalog->getPhases();
  unsigned relocatedEvs                      = 0;
  vector<double> allRms;

  //
  // loop through each event with its corresponding neighbours cluster
  //
  for (const NeighboursPtr &neighbours : neighCluster)
  {
    Event &event = events.at(neighbours->refEvId);

    // get relocation changes computed by the solver for the current event
    double deltaLat, deltaLon, deltaDepth, deltaTT;
    if (!solver.getEventChanges(event.id, deltaLat, deltaLon, deltaDepth,
                                deltaTT))
    {
      event.relocInfo.isRelocated = false;
      continue;
    }

    // check for airquakes (can we trust ttt/vel.model above surface?)
    if (event.depth + deltaDepth < 0)
    {
      // allow 100 meters change
      if (deltaDepth > 0.100)
      {
        SEISCOMP_DEBUG("Ignoring airquake event %s", string(event).c_str());
        continue;
      }
      // do not move the depth in this case
      deltaDepth = 0;
    }

    NeighboursPtr finalNeighbours(new Neighbours());
    finalNeighbours->refEvId = event.id;
    bool isFirstIteration    = !event.relocInfo.isRelocated;

    //
    // update event location/time and compute statistics
    //
    relocatedEvs++;
    event.latitude += deltaLat;
    event.longitude += deltaLon;
    event.depth += deltaDepth;
    event.time += Core::TimeSpan(deltaTT);
    event.rms                    = 0;
    event.relocInfo.isRelocated  = true;
    event.relocInfo.phases       = {0};
    event.relocInfo.neighbours   = {0};
    event.relocInfo.ddObs.numTTp = 0;
    event.relocInfo.ddObs.numTTs = 0;
    event.relocInfo.ddObs.numCCp = 0;
    event.relocInfo.ddObs.numCCs = 0;

    set<unsigned> neighbourIds;
    vector<double> obsResiduals;
    unsigned rmsCount = 0;
    auto eqlrng       = phases.equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      Phase &phase           = it->second;
      const Station &station = stations.at(phase.stationId);
      char phaseTypeAsChar   = static_cast<char>(phase.procInfo.type);

      phase.relocInfo.isRelocated = false;

      unsigned startTTObs, startCCObs, finalTotalObs;
      double meanObsResidual, meanAPrioriWeight, meanFinalWeight;

      if (!solver.getObservationParamsChanges(
              event.id, station.id, phaseTypeAsChar, startTTObs, startCCObs,
              finalTotalObs, meanAPrioriWeight, meanFinalWeight,
              meanObsResidual, neighbourIds))
      {
        continue;
      }

      if (finalTotalObs == 0) continue;

      phase.relocInfo.isRelocated = true;
      phase.relocInfo.finalWeight =
          meanFinalWeight / pickWeightScaler; // range 0-1
      phase.relocInfo.numTTObs = startTTObs;
      phase.relocInfo.numCCObs = startCCObs;
      if (isFirstIteration)
        phase.relocInfo.startMeanObsResidual = meanObsResidual;
      phase.relocInfo.finalMeanObsResidual = meanObsResidual;
      obsResiduals.push_back(meanObsResidual);

      if (obsparams.add(_ttt, event, station, phase, true))
      {
        double travelTime =
            obsparams.get(event.id, station.id, phaseTypeAsChar).travelTime;
        phase.relocInfo.finalResidual =
            travelTime - (phase.time - event.time).length();
        rmsCount++;
      }
      else
      {
        phase.relocInfo.finalResidual = 0;
      }

      event.rms +=
          (phase.relocInfo.finalResidual * phase.relocInfo.finalResidual);
      if (phase.procInfo.type == Phase::Type::P)
      {
        event.relocInfo.phases.usedP++;
        event.relocInfo.ddObs.numCCp += phase.relocInfo.numCCObs;
        event.relocInfo.ddObs.numTTp += phase.relocInfo.numTTObs;
      }
      if (phase.procInfo.type == Phase::Type::S)
      {
        event.relocInfo.phases.usedS++;
        event.relocInfo.ddObs.numCCs += phase.relocInfo.numCCObs;
        event.relocInfo.ddObs.numTTs += phase.relocInfo.numTTObs;
      }

      for (unsigned nId : neighbourIds)
      {
        finalNeighbours->add(nId, station.id, phase.procInfo.type);
      }
    }

    if (rmsCount > 0)
    {
      event.rms = std::sqrt(event.rms / rmsCount);
      allRms.push_back(event.rms);
    }

    double residualMedian = computeMedian(obsResiduals);
    double residualMAD =
        computeMedianAbsoluteDeviation(obsResiduals, residualMedian);
    if (isFirstIteration)
    {
      event.relocInfo.ddObs.startResidualMedian = residualMedian;
      event.relocInfo.ddObs.startResidualMAD    = residualMAD;
    }
    event.relocInfo.ddObs.finalResidualMedian = residualMedian;
    event.relocInfo.ddObs.finalResidualMAD    = residualMAD;

    event.relocInfo.neighbours.amount = finalNeighbours->numNeighbours();
    finalNeighCluster[finalNeighbours->refEvId] = finalNeighbours;
  }

  const double allRmsMedian = computeMedian(allRms);
  const double allRmsMAD = computeMedianAbsoluteDeviation(allRms, allRmsMedian);

  SEISCOMP_INFO(
      "Successfully relocated %u events. RMS median %.4f [sec] median "
      "absolute deviation %.4f [sec]",
      relocatedEvs, allRmsMedian, allRmsMAD);

  return new Catalog(stations, events, phases);
}

CatalogPtr HypoDD::updateRelocatedEventsFinalStats(
    const CatalogCPtr &startCatalog,
    const CatalogCPtr &finalCatalog,
    const std::unordered_map<unsigned, NeighboursPtr> &neighCluster) const
{
  CatalogPtr catalogToReturn(new Catalog());
  vector<double> allRms;
  vector<double> stationDist;

  for (const auto &kv : neighCluster)
  {
    const NeighboursPtr &neighbours = kv.second;
    auto it = finalCatalog->getEvents().find(neighbours->refEvId);

    // If the event hasn't been relocated, remove it from the final catalog.
    if (it == finalCatalog->getEvents().end() ||
        !it->second.relocInfo.isRelocated)
    {
      continue;
    }

    CatalogPtr tmpCat = finalCatalog->extractEvent(neighbours->refEvId, true);

    const Event &startEvent = startCatalog->getEvents().at(neighbours->refEvId);
    Event finalEvent        = tmpCat->getEvents().at(neighbours->refEvId);

    //
    // Compute location/time change of start/final origin.
    //
    finalEvent.relocInfo.locChange =
        computeDistance(finalEvent.latitude, finalEvent.longitude,
                        startEvent.latitude, startEvent.longitude);
    finalEvent.relocInfo.depthChange = finalEvent.depth - startEvent.depth;
    finalEvent.relocInfo.timeChange =
        (finalEvent.time - startEvent.time).length();

    //
    // Compute starting event rms considering only the phases in the final
    // catalog.
    //
    unsigned rmsCount             = 0;
    finalEvent.relocInfo.startRms = 0;
    auto eqlrng = tmpCat->getPhases().equal_range(finalEvent.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      Phase finalPhase       = it->second;
      const Station &station = tmpCat->getStations().at(finalPhase.stationId);
      try
      {
        double travelTime;
        _ttt->compute(startEvent, station,
                      string(1, static_cast<char>(finalPhase.procInfo.type)),
                      travelTime);
        double residual =
            travelTime - (finalPhase.time - startEvent.time).length();
        finalPhase.relocInfo.startResidual = residual;
        tmpCat->updatePhase(finalPhase, false);
        finalEvent.relocInfo.startRms += residual * residual;
        rmsCount++;
      }
      catch (exception &e)
      {
        SEISCOMP_WARNING(
            "Travel Time Table error: %s (Event lat %.6f lon %.6f depth %.6f "
            "Station lat %.6f lon %.6f elevation %.f )",
            e.what(), startEvent.latitude, startEvent.longitude,
            startEvent.depth, station.latitude, station.longitude,
            station.elevation);
      }
    }

    if (rmsCount > 0)
    {
      finalEvent.relocInfo.startRms =
          std::sqrt(finalEvent.relocInfo.startRms / rmsCount);
      allRms.push_back(finalEvent.relocInfo.startRms);
    }

    //
    // compute station distances to final event
    //
    stationDist.clear();
    for (const auto &kv : neighbours->allPhases())
    {
      const string &stationId = kv.first;
      for (Phase::Type phaseType : kv.second)
      {
        // Check if this station is actually part of the event phases (remember
        // keepUnmatchedPhases option during clustering). Make sure that the
        // phase was used for relocation.
        auto it = tmpCat->searchPhase(finalEvent.id, stationId, phaseType);
        if (it != tmpCat->getPhases().end() && it->second.relocInfo.isRelocated)
        {
          const Station &station =
              tmpCat->getStations().at(it->second.stationId);
          stationDist.push_back(computeDistance(finalEvent, station));
          break;
        }
      }
    }

    finalEvent.relocInfo.phases.stationDistMedian = computeMedian(stationDist);
    auto min_max = std::minmax_element(stationDist.begin(), stationDist.end());
    finalEvent.relocInfo.phases.stationDistMin = *min_max.first;
    finalEvent.relocInfo.phases.stationDistMax = *min_max.second;

    //
    // compute starting neighbours to starting event distance
    //
    vector<double> latDiff, lonDiff, depthDiff;
    for (unsigned neighEvId : neighbours->ids)
    {
      const Event &neighEv = startCatalog->getEvents().at(neighEvId);
      latDiff.push_back(startEvent.latitude - neighEv.latitude);
      lonDiff.push_back(startEvent.longitude - neighEv.longitude);
      depthDiff.push_back(startEvent.depth - neighEv.depth);
    }
    double meanLatDiff   = computeMean(latDiff);
    double meanLonDiff   = computeMean(lonDiff);
    double meanDepthDiff = computeMean(depthDiff);
    double centroidLat   = startEvent.latitude + meanLatDiff;
    double centroidLon   = startEvent.longitude + meanLonDiff;
    // double centroidDepth = startEvent.depth + meanDepthDiff; not used

    finalEvent.relocInfo.neighbours.eventDistToCentroid = computeDistance(
        startEvent.latitude, startEvent.longitude, centroidLat, centroidLon);
    finalEvent.relocInfo.neighbours.eventDepthDistToCentroid = meanDepthDiff;
    finalEvent.relocInfo.neighbours.meanDepthDistToCentroid =
        computeMeanAbsoluteDeviation(depthDiff, meanDepthDiff);

    finalEvent.relocInfo.neighbours.meanDistToCentroid = 0;
    for (unsigned neighEvId : neighbours->ids)
    {
      const Event &neighEv = startCatalog->getEvents().at(neighEvId);
      finalEvent.relocInfo.neighbours.meanDistToCentroid += computeDistance(
          neighEv.latitude, neighEv.longitude, centroidLat, centroidLon);
    }
    finalEvent.relocInfo.neighbours.meanDistToCentroid /=
        neighbours->numNeighbours();

    tmpCat->updateEvent(finalEvent, false);
    catalogToReturn->add(finalEvent.id, *tmpCat, true);
  }

  const double allRmsMedian = computeMedian(allRms);
  const double allRmsMAD = computeMedianAbsoluteDeviation(allRms, allRmsMedian);

  SEISCOMP_INFO("Events RMS before relocation: median %.4f median absolute "
                "deviation %.4f",
                allRmsMedian, allRmsMAD);

  return catalogToReturn;
}

void HypoDD::addMissingEventPhases(const Event &refEv,
                                   CatalogPtr &refEvCatalog,
                                   const CatalogCPtr &searchCatalog,
                                   const NeighboursPtr &neighbours)
{
  std::vector<Phase> newPhases =
      findMissingEventPhases(refEv, refEvCatalog, searchCatalog, neighbours);

  for (Phase &ph : newPhases)
  {
    refEvCatalog->updatePhase(ph, true);
    const Station &station = searchCatalog->getStations().at(ph.stationId);
    refEvCatalog->addStation(station);
  }
}

std::vector<Phase>
HypoDD::findMissingEventPhases(const Event &refEv,
                               CatalogPtr &refEvCatalog,
                               const CatalogCPtr &searchCatalog,
                               const NeighboursPtr &neighbours)
{
  //
  // find stations for which the `refEv` doesn't have phases
  //
  vector<MissingStationPhase> missingPhases =
      getMissingPhases(refEv, refEvCatalog, searchCatalog);

  //
  // for each missed phase try to detect it
  //
  std::vector<Phase> newPhases;
  for (const MissingStationPhase &pair : missingPhases)
  {
    const Station &station      = searchCatalog->getStations().at(pair.first);
    const Phase::Type phaseType = pair.second;

    //
    // Loop through every other event and select the ones who have a manually
    // picked phase for the missing station.
    //
    vector<HypoDD::PhasePeer> peers =
        findPhasePeers(station, phaseType, searchCatalog, neighbours);
    if (peers.size() <= 0)
    {
      continue;
    }

    // compute velocity using existing background catalog phases
    double phaseVelocity = 0;
    for (const PhasePeer &peer : peers)
    {
      const Event &event     = peer.first;
      const Phase &phase     = peer.second;
      double travelTime      = (phase.time - event.time).length();
      double stationDistance = computeDistance(event, station);
      double vel             = stationDistance / travelTime;
      phaseVelocity += vel;
    }
    phaseVelocity /= peers.size();

    Phase refEvNewPhase =
        createThoreticalPhase(station, phaseType, refEv, peers, phaseVelocity);

    newPhases.push_back(refEvNewPhase);
  }

  return newPhases;
}

vector<HypoDD::MissingStationPhase>
HypoDD::getMissingPhases(const Event &refEv,
                         CatalogPtr &refEvCatalog,
                         const CatalogCPtr &searchCatalog) const
{
  const auto &refEvPhases = refEvCatalog->getPhases().equal_range(refEv.id);

  //
  // Loop through stations and find those for which the `refEv` doesn't have
  // phases.
  //
  vector<MissingStationPhase> missingPhases;
  for (const auto &kv : searchCatalog->getStations())
  {
    const Station &station = kv.second;

    bool foundP = false, foundS = false;
    for (auto it = refEvPhases.first; it != refEvPhases.second; ++it)
    {
      const Phase &phase = it->second;
      if (station.networkCode == phase.networkCode &&
          station.stationCode == phase.stationCode &&
          station.locationCode == phase.locationCode)
      {
        if (phase.procInfo.type == Phase::Type::P) foundP = true;
        if (phase.procInfo.type == Phase::Type::S) foundS = true;
      }
      if (foundP and foundS) break;
    }
    if (!foundP || !foundS)
    {
      if (!foundP)
        missingPhases.push_back(
            MissingStationPhase(station.id, Phase::Type::P));
      if (!foundS)
        missingPhases.push_back(
            MissingStationPhase(station.id, Phase::Type::S));
    }
  }

  return missingPhases;
}

vector<HypoDD::PhasePeer>
HypoDD::findPhasePeers(const Station &station,
                       const Phase::Type &phaseType,
                       const CatalogCPtr &searchCatalog,
                       const NeighboursPtr &neighbours) const
{
  //
  // Loop through every other event and select those manual phases of the
  // `station` we are interested in.
  //
  vector<PhasePeer> phasePeers;

  for (unsigned neighEvId : neighbours->ids)
  {
    const Event &event = searchCatalog->getEvents().at(neighEvId);

    if (neighbours->has(neighEvId, station.id, phaseType))
    {
      const Phase &phase =
          searchCatalog->searchPhase(neighEvId, station.id, phaseType)->second;

      if (station.networkCode == phase.networkCode &&
          station.stationCode == phase.stationCode &&
          station.locationCode == phase.locationCode &&
          phaseType == phase.procInfo.type)
      {
        if (phase.isManual)
        {
          phasePeers.push_back(PhasePeer(event, phase));
        }
        break;
      }
    }
  }

  return phasePeers;
}

Phase HypoDD::createThoreticalPhase(const Station &station,
                                    const Phase::Type &phaseType,
                                    const Event &refEv,
                                    const vector<HypoDD::PhasePeer> &peers,
                                    double phaseVelocity)
{
  const auto xcorrCfg = _cfg.xcorr.at(phaseType);

  // store most recent `channelCode` used
  struct
  {
    string channelCode;
    Core::Time time;
  } streamInfo = {"", Core::Time()};

  for (const PhasePeer &peer : peers)
  {
    // const Event& event = peer.first;
    const Phase &phase = peer.second;
    // get the closest stream to the `refEv` (w.r.t. time)
    if ((refEv.time - phase.time).abs() < (refEv.time - streamInfo.time).abs())
      streamInfo = {phase.channelCode, phase.time};
  }

  // initialize the new phase
  Phase refEvNewPhase;

  refEvNewPhase.eventId      = refEv.id;
  refEvNewPhase.stationId    = station.id;
  refEvNewPhase.networkCode  = station.networkCode;
  refEvNewPhase.stationCode  = station.stationCode;
  refEvNewPhase.locationCode = station.locationCode;
  refEvNewPhase.channelCode =
      getBandAndInstrumentCodes(streamInfo.channelCode) +
      xcorrCfg.components[0];
  refEvNewPhase.isManual      = false;
  refEvNewPhase.procInfo.type = phaseType;

  // use phase velocity to compute phase time
  double stationDistance = computeDistance(refEv, station);
  refEvNewPhase.time =
      refEv.time + Core::TimeSpan(stationDistance / phaseVelocity);

  refEvNewPhase.lowerUncertainty = Catalog::DEFAULT_AUTOMATIC_PICK_UNCERTAINTY;
  refEvNewPhase.upperUncertainty = Catalog::DEFAULT_AUTOMATIC_PICK_UNCERTAINTY;
  refEvNewPhase.procInfo.weight  = Catalog::computePickWeight(refEvNewPhase);
  refEvNewPhase.procInfo.source  = Phase::Source::THEORETICAL;
  refEvNewPhase.type = stringify("%ct", static_cast<char>(phaseType));

  return refEvNewPhase;
}

XCorrCache HypoDD::buildXCorrCache(CatalogPtr &catalog,
                                   const std::list<NeighboursPtr> &neighCluster,
                                   bool computeTheoreticalPhases,
                                   double xcorrMaxEvStaDist,
                                   double xcorrMaxInterEvDist)
{
  XCorrCache xcorr;
  resetCounters();

  unsigned long performed = 0;

  for (const NeighboursPtr &neighbours : neighCluster)
  {
    const Event &refEv = catalog->getEvents().at(neighbours->refEvId);

    // Compute theoretical phases for stations that have no picks. The
    // cross-correlation will be used to detect and fix the pick time.
    if (computeTheoreticalPhases)
    {
      addMissingEventPhases(refEv, catalog, catalog, neighbours);
    }

    buildXcorrDiffTTimePairs(catalog, neighbours, refEv, xcorrMaxEvStaDist,
                             xcorrMaxInterEvDist, xcorr);

    // Update theoretical and automatic phase pick time and uncertainties based
    // on cross-correlation results. Also, drop theoretical phases wihout any
    // good cross-correlation result.
    fixPhases(catalog, refEv, xcorr);

    SEISCOMP_INFO("Cross-correlation completion %.1f%%",
                  (++performed / (double)neighCluster.size()) * 100);
  }

  printCounters();

  return xcorr;
}

/*
 * Compute and store to `XCorrCache` cross-correlated differential travel times
 * for pairs of the earthquake.
 */
void HypoDD::buildXcorrDiffTTimePairs(CatalogPtr &catalog,
                                      const NeighboursPtr &neighbours,
                                      const Event &refEv,
                                      double xcorrMaxEvStaDist,
                                      double xcorrMaxInterEvDist,
                                      XCorrCache &xcorr)
{
  SEISCOMP_INFO(
      "Computing cross-correlation differential travel times for event %s",
      string(refEv).c_str());

  //
  // Prepare the waveform loaders for single-event. Since they are
  // temporary and discarded after the relocation we don't want to
  // make use of the background catalog caching
  //
  Waveform::LoaderPtr seLdr = preloadNonCatalogWaveforms(
      catalog, neighbours, refEv, xcorrMaxEvStaDist, xcorrMaxInterEvDist);

  Waveform::LoaderPtr seWfLdr = seLdr;
  if (_cfg.snr.minSnr > 0)
  {
    seWfLdr = new Waveform::SnrFilteredLoader(
        seWfLdr, _cfg.snr.minSnr, _cfg.snr.noiseStart, _cfg.snr.noiseEnd,
        _cfg.snr.signalStart, _cfg.snr.signalEnd);
    if (_waveformDebug) seWfLdr->setDebugDirectory(_wfDebugDir);
  }
  seWfLdr = new Waveform::MemCachedLoader(seWfLdr);
  if (_waveformDebug) seWfLdr->setDebugDirectory(_wfDebugDir);

  Waveform::LoaderPtr seWfLdrNoSnr = new Waveform::MemCachedLoader(seLdr);
  if (_waveformDebug) seWfLdrNoSnr->setDebugDirectory(_wfDebugDir);

  // keep track of the `refEv` distance to stations
  multimap<double, string> stationByDistance; // <distance, stationid>
  unordered_set<string> computedStations;

  //
  // loop through reference event phases
  //
  auto eqlrngRef = catalog->getPhases().equal_range(refEv.id);
  for (auto itRef = eqlrngRef.first; itRef != eqlrngRef.second; ++itRef)
  {
    const Phase &refPhase  = itRef->second;
    const Station &station = catalog->getStations().at(refPhase.stationId);

    //
    // skip stations too far away
    //
    double stationDistance = computeDistance(refEv, station);
    if (stationDistance > xcorrMaxEvStaDist && xcorrMaxEvStaDist >= 0) continue;

    //
    // select appropriate waveform loader for refPhase
    //
    Waveform::LoaderPtr refLdr;
    if (refPhase.procInfo.source == Phase::Source::CATALOG)
      refLdr = _wfAccess.memCache;
    else if (refPhase.isManual)
      refLdr = seWfLdr;
    else // not from catalog, not manual
    {
      // For "untrusted" phases the SNR is checked AFTER the cross-correlation
      // once the pick time has been adjusted using the lag computed by the
      // cross-corr against a "trusted" phase (manual or catalog).
      refLdr = seWfLdrNoSnr;
    }

    //
    // loop through neighbouring events and cross-correlate with `refPhase`
    //
    for (unsigned neighEvId : neighbours->ids)
    {
      const Event &event = catalog->getEvents().at(neighEvId);

      //
      // skip events too far away
      //
      double interEventDistance = computeDistance(refEv, event);
      if (interEventDistance > xcorrMaxInterEvDist && xcorrMaxInterEvDist >= 0)
        continue;

      if (neighbours->has(neighEvId, refPhase.stationId,
                          refPhase.procInfo.type))
      {
        const Phase &phase = catalog
                                 ->searchPhase(event.id, refPhase.stationId,
                                               refPhase.procInfo.type)
                                 ->second;

        // In single-event mode `refPhase` is real-time and `phase` is from the
        // catalog. In multi-event mode both are from the catalog.
        if (phase.procInfo.source != Phase::Source::CATALOG)
          throw runtime_error(
              "Internal logic error: phase is not from catalog");

        double coeff, lag;
        if (xcorrPhases(refEv, refPhase, refLdr, event, phase,
                        _wfAccess.memCache, coeff, lag))
        {
          bool goodSNR = true;

          // Check the SNR (using the pick time adjusted with the
          // cross-correlation lag) of those phases, where the check hadn't
          // been performed, yet.
          if (_cfg.snr.minSnr > 0 && !refPhase.isManual &&
              refPhase.procInfo.source != Phase::Source::CATALOG)
          {
            const auto xcorrCfg = _cfg.xcorr.at(refPhase.procInfo.type);

            // Make sure that at least one of the components allowed for this
            // phase type has a good SNR.
            goodSNR = false;
            for (const string &component : xcorrCfg.components)
            {
              Phase tmpPh = refPhase;
              tmpPh.channelCode =
                  getBandAndInstrumentCodes(tmpPh.channelCode) + component;
              Core::Time adjustedPickTime = refPhase.time - Core::TimeSpan(lag);
              Core::TimeWindow snrWin =
                  _wfAccess.snrFilter->snrTimeWindow(adjustedPickTime);
              GenericRecordCPtr trace = seWfLdrNoSnr->get(
                  snrWin, tmpPh, refEv, true, _cfg.wfFilter.filterStr,
                  _cfg.wfFilter.resampleFreq);
              if (trace &&
                  _wfAccess.snrFilter->goodSnr(trace, adjustedPickTime))
              {
                goodSNR = true;
                break;
              }
            }
          }

          // store good cross-correlation results
          if (goodSNR)
          {
            auto &entry1 = xcorr.getForUpdate(refEv.id, refPhase.stationId,
                                              refPhase.procInfo.type);
            entry1.update(event, phase, coeff, lag);
            auto &entry2 = xcorr.getForUpdate(event.id, phase.stationId,
                                              phase.procInfo.type);
            entry2.update(refEv, refPhase, coeff, lag);
          }
        }

        // keep track of events/station distance for every cross-correlation
        // performed
        if (computedStations.find(refPhase.stationId) == computedStations.end())
        {
          stationByDistance.emplace(stationDistance, refPhase.stationId);
          computedStations.insert(refPhase.stationId);
        }
      }
    }

    // finalize statistics
    if (xcorr.has(refEv.id, refPhase.stationId, refPhase.procInfo.type))
    {
      auto &entry = xcorr.getForUpdate(refEv.id, refPhase.stationId,
                                       refPhase.procInfo.type);
      entry.computeStats();
    }
  }

  // keep track of couters
  updateCounters(nullptr, nullptr,
                 dynamic_cast<Waveform::SnrFilteredLoader *>(seWfLdr.get()));

  // print some useful information
  for (const auto &kv : stationByDistance)
  {
    const double stationDistance = kv.first;
    const Station &station       = catalog->getStations().at(kv.second);

    bool goodPXcorr = xcorr.has(refEv.id, station.id, Phase::Type::P);
    bool goodSXcorr = xcorr.has(refEv.id, station.id, Phase::Type::S);

    if (!goodPXcorr && !goodSXcorr)
    {
      SEISCOMP_INFO(
          "xcorr: event %5s sta %4s %5s dist %7.2f [km] - low corr coeff pairs",
          string(refEv).c_str(), station.networkCode.c_str(),
          station.stationCode.c_str(), stationDistance);
    }
    else
    {
      if (goodPXcorr)
      {
        const auto &pdata = xcorr.get(refEv.id, station.id, Phase::Type::P);
        SEISCOMP_INFO("xcorr: event %5s sta %4s %5s dist %7.2f [km] - "
                      "%d P phases, mean coeff %.2f lag %.2f (events: %s)",
                      string(refEv).c_str(), station.networkCode.c_str(),
                      station.stationCode.c_str(), stationDistance,
                      pdata.ccCount, pdata.mean_coeff, pdata.mean_lag,
                      pdata.peersStr.c_str());
      }
      if (goodSXcorr)
      {
        const auto &sdata = xcorr.get(refEv.id, station.id, Phase::Type::S);
        SEISCOMP_INFO("xcorr: event %5s sta %4s %5s dist %7.2f [km] - "
                      "%d S phases, mean coeff %.2f lag %.2f (events: %s)",
                      string(refEv).c_str(), station.networkCode.c_str(),
                      station.stationCode.c_str(), stationDistance,
                      sdata.ccCount, sdata.mean_coeff, sdata.mean_lag,
                      sdata.peersStr.c_str());
      }
    }
  }
}

Waveform::LoaderPtr
HypoDD::preloadNonCatalogWaveforms(CatalogPtr &catalog,
                                   const NeighboursPtr &neighbours,
                                   const Event &refEv,
                                   double xcorrMaxEvStaDist,
                                   double xcorrMaxInterEvDist) const
{
  //
  // For single-event relocation in real-time we want to load the waveforms
  // in batch otherwise the seedlink server gets stuck and becomes unresponsive
  // due to the multiple connections requests, one for each event phase
  //
  Waveform::BatchLoaderPtr batchLoader =
      new Waveform::BatchLoader(_cfg.recordStreamURL);

  Waveform::LoaderPtr returnedLdr = batchLoader;

  // make sure to use disk cache if enabled for non catalog phases
  Waveform::DiskCachedLoaderPtr diskLoader;
  if (_useCatalogWaveformDiskCache && _waveformCacheAll)
  {
    diskLoader  = new Waveform::DiskCachedLoader(batchLoader, _tmpCacheDir);
    returnedLdr = new Waveform::ExtraLenLoader(diskLoader, DISK_TRACE_MIN_LEN);
  }

  // Load a large enough waveform to be able to check the SNR after
  // the updating of the pick time
  if (_cfg.snr.minSnr > 0)
  {
    double maxDelay = std::max({_cfg.xcorr.at(Phase::Type::P).maxDelay,
                                _cfg.xcorr.at(Phase::Type::S).maxDelay});
    double extraBefore =
        std::min({_cfg.snr.noiseStart, _cfg.snr.signalStart}) - maxDelay;
    extraBefore = extraBefore > 0 ? 0 : std::abs(extraBefore);
    double extraAfter =
        std::max({_cfg.snr.noiseEnd, _cfg.snr.signalEnd}) + maxDelay;
    extraAfter = extraAfter < 0 ? 0 : extraAfter;
    returnedLdr =
        new Waveform::ExtraLenLoader(returnedLdr, extraBefore, extraAfter);
  }

  //
  // loop through reference event phases
  //
  auto eqlrngRef = catalog->getPhases().equal_range(refEv.id);
  for (auto itRef = eqlrngRef.first; itRef != eqlrngRef.second; ++itRef)
  {
    const Phase &refPhase  = itRef->second;
    const Station &station = catalog->getStations().at(refPhase.stationId);

    // We deal only with real-time event data
    if (refPhase.procInfo.source == Phase::Source::CATALOG) continue;

    //
    // skip stations too far away
    //
    double stationDistance = computeDistance(refEv, station);
    if (stationDistance > xcorrMaxEvStaDist && xcorrMaxEvStaDist >= 0) continue;

    //
    // loop through neighbouring events and cross-correlate with `refPhase`
    //
    for (unsigned neighEvId : neighbours->ids)
    {
      const Event &event = catalog->getEvents().at(neighEvId);

      //
      // skip events too far away
      //
      double interEventDistance = computeDistance(refEv, event);
      if (interEventDistance > xcorrMaxInterEvDist && xcorrMaxInterEvDist >= 0)
        continue;

      if (neighbours->has(neighEvId, refPhase.stationId,
                          refPhase.procInfo.type))
      {
        //
        // For each match load the reference event phase waveforms
        //
        Core::TimeWindow tw = xcorrTimeWindowLong(refPhase);
        const auto xcorrCfg = _cfg.xcorr.at(refPhase.procInfo.type);

        for (string component : xcorrCfg.components)
        {
          Phase tmpPh = refPhase;
          tmpPh.channelCode =
              getBandAndInstrumentCodes(tmpPh.channelCode) + component;

          // This doesn't really load the trace but force the request to reach
          // the loader, which will load all the traces later
          returnedLdr->get(tw, tmpPh, refEv);
        }
      }
    }
  }

  // This will actually dowanload the waveworms all at once
  batchLoader->load();

  // keep track of couters
  updateCounters(batchLoader, diskLoader, nullptr);

  return returnedLdr;
}

/*
 * Update theoretical and automatic phase pick times and uncertainties based on
 * cross-correlation results. Drop theoretical phases not passing the
 * cross-correlation verification.
 */
void HypoDD::fixPhases(CatalogPtr &catalog,
                       const Event &refEv,
                       XCorrCache &xcorr)
{
  unsigned totP = 0, totS = 0;
  unsigned newP = 0, newS = 0;

  std::vector<Phase> phasesToBeRemoved;
  std::vector<Phase> newPhases;

  auto eqlrng = catalog->getPhases().equal_range(refEv.id);
  for (auto it = eqlrng.first; it != eqlrng.second; it++)
  {
    const Phase &phase = it->second;
    bool goodXcorr = xcorr.has(refEv.id, phase.stationId, phase.procInfo.type);

    if (phase.procInfo.type == Phase::Type::P) totP++;
    if (phase.procInfo.type == Phase::Type::S) totS++;

    // nothing to do if we dont't have good xcorr results of if the phase is
    // manual or from catalog
    if (!goodXcorr || phase.isManual ||
        (phase.procInfo.source == Phase::Source::CATALOG))
    {
      // remove thoretical phases without good cross-correlation results
      if (phase.procInfo.source == Phase::Source::THEORETICAL)
        phasesToBeRemoved.push_back(phase);
      continue;
    }

    const auto &pdata =
        xcorr.get(refEv.id, phase.stationId, phase.procInfo.type);

    //
    // set new phase time and uncertainty
    //
    Phase newPhase(phase);
    newPhase.time -= Core::TimeSpan(pdata.mean_lag);
    newPhase.lowerUncertainty = pdata.mean_lag - pdata.min_lag;
    newPhase.upperUncertainty = pdata.max_lag - pdata.mean_lag;
    newPhase.procInfo.weight  = Catalog::computePickWeight(newPhase);
    newPhase.procInfo.source  = Phase::Source::XCORR;
    newPhase.type = stringify("%cx", static_cast<char>(newPhase.procInfo.type));

    if (phase.procInfo.source == Phase::Source::THEORETICAL)
    {
      if (newPhase.procInfo.type == Phase::Type::P) newP++;
      if (newPhase.procInfo.type == Phase::Type::S) newS++;
    }

    newPhases.push_back(newPhase);
  }

  //
  // replace automatic/theoretical phases with those detected by means of
  // cross-correlation
  //
  for (const Phase &ph : phasesToBeRemoved)
  {
    if (ph.procInfo.type == Phase::Type::P) totP--;
    if (ph.procInfo.type == Phase::Type::S) totS--;
    catalog->removePhase(ph.eventId, ph.stationId, ph.procInfo.type);
  }

  for (Phase &ph : newPhases)
  {
    catalog->updatePhase(ph, true);
  }

  SEISCOMP_INFO("Event %s total phases %u (%u P and %u S): created %u (%u P "
                "and %u S) from theoretical picks",
                string(refEv).c_str(), (totP + totS), totP, totS, (newP + newS),
                newP, newS);
}

void HypoDD::resetCounters()
{
  _counters = {0};
  if (_wfAccess.loader)
  {
    _wfAccess.loader->_counters_wf_no_avail   = 0;
    _wfAccess.loader->_counters_wf_downloaded = 0;
  }
  if (_wfAccess.diskCache)
  {
    _wfAccess.diskCache->_counters_wf_cached = 0;
  }
  if (_wfAccess.snrFilter)
  {
    _wfAccess.snrFilter->_counters_wf_snr_low = 0;
  }
}

void HypoDD::updateCounters() const
{
  updateCounters(_wfAccess.loader, _wfAccess.diskCache, _wfAccess.snrFilter);
}

void HypoDD::updateCounters(Waveform::LoaderPtr loader,
                            Waveform::DiskCachedLoaderPtr diskCache,
                            Waveform::SnrFilteredLoaderPtr snrFilter) const
{
  if (loader)
  {
    _counters.wf_downloaded += loader->_counters_wf_downloaded;
    _counters.wf_no_avail += loader->_counters_wf_no_avail;
  }
  if (diskCache) _counters.wf_disk_cached += diskCache->_counters_wf_cached;
  if (snrFilter) _counters.wf_snr_low += snrFilter->_counters_wf_snr_low;
}

void HypoDD::printCounters() const
{
  updateCounters();

  unsigned performed      = _counters.xcorr_performed;
  unsigned performed_s    = _counters.xcorr_performed_s;
  unsigned performed_p    = performed - performed_s;
  unsigned perf_theo      = _counters.xcorr_performed_theo;
  unsigned perf_theo_s    = _counters.xcorr_performed_s_theo;
  unsigned perf_theo_p    = perf_theo - perf_theo_s;
  unsigned good_cc        = _counters.xcorr_good_cc;
  unsigned good_cc_s      = _counters.xcorr_good_cc_s;
  unsigned good_cc_p      = good_cc - good_cc_s;
  unsigned good_cc_theo   = _counters.xcorr_good_cc_theo;
  unsigned good_cc_s_theo = _counters.xcorr_good_cc_s_theo;
  unsigned good_cc_p_theo = good_cc_theo - good_cc_s_theo;

  unsigned wf_snr_low     = _counters.wf_snr_low;
  unsigned wf_no_avail    = _counters.wf_no_avail;
  unsigned wf_disk_cached = _counters.wf_disk_cached;
  unsigned wf_downloaded  = _counters.wf_downloaded;

  SEISCOMP_INFO("Cross-correlation performed %u, "
                "phases with SNR ratio too low %u, "
                "phases not available %u (waveforms downloaded %u, "
                "waveforms loaded from disk cache %u)",
                performed, wf_snr_low, wf_no_avail, wf_downloaded,
                wf_disk_cached);

  SEISCOMP_INFO("Total xcorr %u (P %.f%%, S %.f%%) success %.f%% (%u/%u). "
                "Successful P %.f%% (%u/%u). Successful S %.f%% (%u/%u)",
                performed, (performed_p * 100. / performed),
                (performed_s * 100. / performed), (good_cc * 100. / performed),
                good_cc, performed, (good_cc_p * 100. / performed_p), good_cc_p,
                performed_p, (good_cc_s * 100. / performed_s), good_cc_s,
                performed_s);

  if (perf_theo > 0)
  {
    unsigned perf_real      = performed - perf_theo,
             perf_real_s    = performed_s - perf_theo_s,
             perf_real_p    = perf_real - perf_real_s,
             good_cc_real   = good_cc - good_cc_theo,
             good_cc_s_real = good_cc_s - good_cc_s_theo,
             good_cc_p_real = good_cc_real - good_cc_s_real;

    SEISCOMP_INFO(
        "xcorr on actual picks %u/%u (P %.f%%, S %.f%%) success %.f%% (%u/%u). "
        "Successful P %.f%% (%u/%u). Successful S %.f%% (%u/%u)",
        perf_real, performed, (perf_real_p * 100. / perf_real),
        (perf_real_s * 100. / perf_real), (good_cc_real * 100. / perf_real),
        good_cc_real, perf_real, (good_cc_p_real * 100. / perf_real_p),
        good_cc_p_real, perf_real_p, (good_cc_s_real * 100. / perf_real_s),
        good_cc_s_real, perf_real_s);

    SEISCOMP_INFO(
        "xcorr on theoretical picks %u/%u (P %.f%%, S %.f%%) success %.f%% "
        "(%u/%u). Successful P %.f%% (%u/%u). Successful S %.f%% (%u/%u)",
        perf_theo, performed, (perf_theo_p * 100. / perf_theo),
        (perf_theo_s * 100. / perf_theo), (good_cc_theo * 100. / perf_theo),
        good_cc_theo, perf_theo, (good_cc_p_theo * 100. / perf_theo_p),
        good_cc_p_theo, perf_theo_p, (good_cc_s_theo * 100. / perf_theo_s),
        good_cc_s_theo, perf_theo_s);
  }
}

Core::TimeWindow HypoDD::xcorrTimeWindowLong(const Phase &phase) const
{
  const auto xcorrCfg = _cfg.xcorr.at(phase.procInfo.type);
  Core::TimeWindow tw = xcorrTimeWindowShort(phase);
  tw.setStartTime(tw.startTime() - Core::TimeSpan(xcorrCfg.maxDelay));
  tw.setEndTime(tw.endTime() + Core::TimeSpan(xcorrCfg.maxDelay));
  return tw;
}

Core::TimeWindow HypoDD::xcorrTimeWindowShort(const Phase &phase) const
{
  const auto xcorrCfg  = _cfg.xcorr.at(phase.procInfo.type);
  double shortDuration = xcorrCfg.endOffset - xcorrCfg.startOffset;
  Core::TimeSpan shortTimeCorrection = Core::TimeSpan(xcorrCfg.startOffset);
  return Core::TimeWindow(phase.time + shortTimeCorrection, shortDuration);
}

bool HypoDD::xcorrPhases(const Event &event1,
                         const Phase &phase1,
                         Waveform::LoaderPtr ph1Cache,
                         const Event &event2,
                         const Phase &phase2,
                         Waveform::LoaderPtr ph2Cache,
                         double &coeffOut,
                         double &lagOut)
{
  if (phase1.procInfo.type != phase2.procInfo.type)
  {
    SEISCOMP_ERROR(
        "Internal logic error: trying to cross-correlate mismatching "
        "phases (%s and %s)",
        string(phase1).c_str(), string(phase2).c_str());
    return false;
  }

  auto xcorrCfg  = _cfg.xcorr.at(phase1.procInfo.type);
  bool performed = false;
  bool goodCoeff = false;

  //
  // Try to use the same channels for the cross-correlation. In case the two
  // phases differ, do not change the catalog phase channels.
  //
  const string channelCodeRoot1 = getBandAndInstrumentCodes(phase1.channelCode);
  const string channelCodeRoot2 = getBandAndInstrumentCodes(phase2.channelCode);

  string commonChRoot;

  if (channelCodeRoot1 == channelCodeRoot2)
  {
    commonChRoot = channelCodeRoot1;
  }
  else if (phase1.procInfo.source == Phase::Source::CATALOG &&
           phase2.procInfo.source != Phase::Source::CATALOG)
  {
    DataModel::ThreeComponents dummy;
    DataModel::SensorLocation *loc2 =
        ScCatalog::findSensorLocation(phase2.networkCode, phase2.stationCode,
                                      phase2.locationCode, phase2.time);
    if (loc2 &&
        getThreeComponents(dummy, loc2, channelCodeRoot1.c_str(), phase2.time))
    {
      // phase 2 has the same channels as phase 1
      commonChRoot = channelCodeRoot1;
    }
  }
  else if (phase1.procInfo.source != Phase::Source::CATALOG &&
           phase2.procInfo.source == Phase::Source::CATALOG)
  {
    DataModel::ThreeComponents dummy;
    DataModel::SensorLocation *loc1 =
        ScCatalog::findSensorLocation(phase1.networkCode, phase1.stationCode,
                                      phase1.locationCode, phase1.time);
    if (loc1 &&
        getThreeComponents(dummy, loc1, channelCodeRoot2.c_str(), phase1.time))
    {
      // phase 1 has the same channels as phase 2
      commonChRoot = channelCodeRoot2;
    }
  }

  if (commonChRoot.empty())
  {
    SEISCOMP_DEBUG(
        "Cannot find common channels to cross-correlate %s and %s (%s and %s)",
        channelCodeRoot1.c_str(), channelCodeRoot2.c_str(),
        string(phase1).c_str(), string(phase2).c_str());
  }

  //
  // Perform the cross-correlation on all registered components until we get a
  // good correlation coefficient.
  //
  for (const string &component : xcorrCfg.components)
  {
    Phase tmpPh1 = phase1;
    Phase tmpPh2 = phase2;

    // overwrite phases' component for the cross-correlation
    if (commonChRoot.empty())
    {
      tmpPh1.channelCode = channelCodeRoot1 + component;
      tmpPh2.channelCode = channelCodeRoot2 + component;
    }
    else
    {
      tmpPh1.channelCode = commonChRoot + component;
      tmpPh2.channelCode = commonChRoot + component;
    }

    performed = _xcorrPhases(event1, tmpPh1, ph1Cache, event2, tmpPh2, ph2Cache,
                             coeffOut, lagOut);

    coeffOut = std::abs(coeffOut);

    goodCoeff = (performed && coeffOut >= xcorrCfg.minCoef);

    // If the cross-correlation was successful and the coefficient is good,
    // stop here.
    if (goodCoeff)
    {
      break;
    }
  }

  //
  // deal with counters
  //
  bool isS           = (phase1.procInfo.type == Phase::Type::S);
  bool isTheoretical = (phase1.procInfo.source == Phase::Source::XCORR ||
                        phase2.procInfo.source == Phase::Source::XCORR ||
                        phase1.procInfo.source == Phase::Source::THEORETICAL ||
                        phase2.procInfo.source == Phase::Source::THEORETICAL);

  if (performed)
  {
    _counters.xcorr_performed++;
    if (isTheoretical) _counters.xcorr_performed_theo++;
    if (isS)
    {
      _counters.xcorr_performed_s++;
      if (isTheoretical) _counters.xcorr_performed_s_theo++;
    }

    if (goodCoeff)
    {
      _counters.xcorr_good_cc++;
      if (isTheoretical) _counters.xcorr_good_cc_theo++;
      if (isS)
      {
        _counters.xcorr_good_cc_s++;
        if (isTheoretical) _counters.xcorr_good_cc_s_theo++;
      }
    }
  }

  return goodCoeff;
}

bool HypoDD::_xcorrPhases(const Event &event1,
                          const Phase &phase1,
                          Waveform::LoaderPtr ph1Cache,
                          const Event &event2,
                          const Phase &phase2,
                          Waveform::LoaderPtr ph2Cache,
                          double &coeffOut,
                          double &lagOut)
{
  coeffOut = lagOut = 0;

  auto xcorrCfg = _cfg.xcorr.at(phase1.procInfo.type);

  Core::TimeWindow tw1 = xcorrTimeWindowLong(phase1);
  Core::TimeWindow tw2 = xcorrTimeWindowLong(phase2);

  // Load the long `tr1`, because we want to cache the long version. Then we'll
  // trim it.
  GenericRecordCPtr tr1 = getWaveform(tw1, event1, phase1, ph1Cache);
  if (!tr1)
  {
    return false;
  }

  // Load the long `tr2`, because we want to cache the long version. Then we'll
  // trim it.
  GenericRecordCPtr tr2 = getWaveform(tw2, event2, phase2, ph2Cache);
  if (!tr2)
  {
    return false;
  }

  // Trust the manual pick on `phase2`: keep `tr2` short and cross-correlate it
  // with the larger `tr1` window.
  double xcorr_coeff = 0, xcorr_lag = 0;

  if (phase2.isManual || (!phase1.isManual && !phase2.isManual))
  {
    // Trim `tr2` to shorter length; we want to cross-correlate the short one
    // with the long one.
    GenericRecordPtr tr2Short = new GenericRecord(*tr2);
    Core::TimeWindow tw2Short = xcorrTimeWindowShort(phase2);
    if (!Waveform::trim(*tr2Short, tw2Short))
    {
      SEISCOMP_DEBUG("Cannot trim phase2 waveform, skipping cross-correlation "
                     "for phase pair phase1='%s', phase2='%s'",
                     string(phase1).c_str(), string(phase2).c_str());
      return false;
    }

    if (!Waveform::xcorr(tr1, tr2Short, xcorrCfg.maxDelay, true, xcorr_lag,
                         xcorr_coeff))
    {
      return false;
    }
  }

  // Trust the manual pick on `phase1`: keep `tr1` short and cross-correlate it
  // with a larger `tr2` window.
  double xcorr_coeff2 = 0, xcorr_lag2 = 0;

  if (phase1.isManual || (!phase1.isManual && !phase2.isManual))
  {
    // Trim `tr1` to shorter length; we want to cross-correlate the short with
    // the long one.
    GenericRecordPtr tr1Short = new GenericRecord(*tr1);
    Core::TimeWindow tw1Short = xcorrTimeWindowShort(phase1);
    if (!Waveform::trim(*tr1Short, tw1Short))
    {
      SEISCOMP_DEBUG("Cannot trim phase1 waveform, skipping cross-correlation "
                     "for phase pair phase1='%s', phase2='%s'",
                     string(phase1).c_str(), string(phase2).c_str());
      return false;
    }

    if (!Waveform::xcorr(tr1Short, tr2, xcorrCfg.maxDelay, true, xcorr_lag2,
                         xcorr_coeff2))
    {
      return false;
    }
  }

  if (std::abs(xcorr_coeff2) > std::abs(xcorr_coeff))
  {
    // swap
    xcorr_coeff = xcorr_coeff2;
    xcorr_lag   = xcorr_lag2;
  }

  coeffOut = xcorr_coeff;
  lagOut   = xcorr_lag;

  return true;
}

GenericRecordCPtr HypoDD::getWaveform(const Core::TimeWindow &tw,
                                      const Catalog::Event &ev,
                                      const Catalog::Phase &ph,
                                      Waveform::LoaderPtr wfLoader)
{
  string wfDesc = stringify(
      "Waveform for Phase '%s' and Time slice from %s length %.2f sec",
      string(ph).c_str(), tw.startTime().iso().c_str(), tw.length());

  const string wfId = Waveform::waveformId(ph, tw);

  // Check if we have already excluded the trace because we couldn't load it
  // (-> save time).
  if (_wfAccess.unloadableWfs.find(wfId) != _wfAccess.unloadableWfs.end())
  {
    return nullptr;
  }

  // try to load the waveform
  GenericRecordCPtr trace = wfLoader->get(
      tw, ph, ev, true, _cfg.wfFilter.filterStr, _cfg.wfFilter.resampleFreq);

  if (!trace)
  {
    _wfAccess.unloadableWfs.insert(wfId);
    return nullptr;
  }

  return trace;
}

namespace {

struct XCorrEvalStats
{
  unsigned total  = 0;
  unsigned goodCC = 0;
  vector<double> ccCoeff;
  vector<double> ccCount;
  vector<double> timeDiff;

  void addBadCC() { total++; }

  void addGoodCC(double ccCoeff, unsigned ccCount, double timeDiff)
  {
    this->total++;
    this->goodCC++;
    this->ccCoeff.push_back(ccCoeff);
    this->ccCount.push_back(ccCount);
    this->timeDiff.push_back(timeDiff);
  }

  XCorrEvalStats &operator+=(XCorrEvalStats const &rhs) &
  {
    total += rhs.total;
    goodCC += rhs.goodCC;
    ccCoeff.insert(ccCoeff.end(), rhs.ccCoeff.begin(), rhs.ccCoeff.end());
    ccCount.insert(ccCount.end(), rhs.ccCount.begin(), rhs.ccCount.end());
    timeDiff.insert(timeDiff.end(), rhs.timeDiff.begin(), rhs.timeDiff.end());
    return *this;
  }

  friend XCorrEvalStats operator+(XCorrEvalStats lhs, XCorrEvalStats const &rhs)
  {
    lhs += rhs;
    return lhs;
  }

  string describeShort() const
  {
    double meanCoeff    = computeMean(ccCoeff);
    double meanCoeffMAD = computeMeanAbsoluteDeviation(ccCoeff, meanCoeff);
    double meanCount    = computeMean(ccCount);
    double meanCountMAD = computeMeanAbsoluteDeviation(ccCount, meanCount);
    string log = stringify("#pha %6d pha good CC %3.f%% coeff %.2f (+/-%.2f) "
                           "goodCC/ph %4.1f (+/-%.1f)",
                           total, (goodCC * 100. / total), meanCoeff,
                           meanCoeffMAD, meanCount, meanCountMAD);
    double meanDev    = computeMean(timeDiff);
    double meanDevMAD = computeMeanAbsoluteDeviation(timeDiff, meanDev);
    log += stringify(" time-diff [msec] %3.f (+/-%.f)", meanDev * 1000,
                     meanDevMAD * 1000);
    return log;
  }

  string describe() const
  {
    double meanCoeff    = computeMean(ccCoeff);
    double meanCoeffMAD = computeMeanAbsoluteDeviation(ccCoeff, meanCoeff);
    double meanCount    = computeMean(ccCount);
    double meanCountMAD = computeMeanAbsoluteDeviation(ccCount, meanCount);
    string log = stringify("%9d %5.f%% % 4.2f (%4.2f)   %4.1f (%4.1f)", total,
                           (goodCC * 100. / total), meanCoeff, meanCoeffMAD,
                           meanCount, meanCountMAD);
    double meanDev    = computeMean(timeDiff);
    double meanDevMAD = computeMeanAbsoluteDeviation(timeDiff, meanDev);
    log += stringify("%8.f (%3.f)", meanDev * 1000, meanDevMAD * 1000);
    return log;
  }
};
} // namespace

void HypoDD::evalXCorr(const ClusteringOptions &clustOpt, bool theoretical)
{
  XCorrEvalStats totalStats, pPhaseStats, sPhaseStats;
  map<string, XCorrEvalStats> statsByStation;      // key station id
  map<int, XCorrEvalStats> statsByInterEvDistance; // key distance
  map<int, XCorrEvalStats> statsByStaDistance;     // key distance
  const double EV_DIST_STEP  = 0.1;                // km
  const double STA_DIST_STEP = 3;                  // km

  auto printStats = [&](string title) {
    string log = title + "\n";
    log +=
        stringify("Cumulative stats: %s\n", totalStats.describeShort().c_str());
    log += stringify("Cumulative stats P ph: %s\n",
                     pPhaseStats.describeShort().c_str());
    log += stringify("Cumulative stats S ph: %s\n",
                     sPhaseStats.describeShort().c_str());

    log += stringify(
        "Cross-correlated phases by inter-event distance in %.2f km step\n",
        EV_DIST_STEP);
    log += stringify(" EvDist [km]  #Phases GoodCC AvgCoeff(+/-) "
                     "GoodCC/Ph(+/-) time-diff[msec] (+/-)\n");
    for (const auto &kv : statsByInterEvDistance)
    {
      log += stringify("%5.2f-%-5.2f %s\n", kv.first * EV_DIST_STEP,
                       (kv.first + 1) * EV_DIST_STEP,
                       kv.second.describe().c_str());
    }

    log += stringify("Cross-correlated phases by event to station distance in "
                     "%.2f km step\n",
                     STA_DIST_STEP);
    log += stringify("StaDist [km]  #Phases GoodCC AvgCoeff(+/-) "
                     "GoodCC/Ph(+/-) time-diff[msec] (+/-)\n");
    for (const auto &kv : statsByStaDistance)
    {
      log += stringify("%3d-%-3d     %s\n", int(kv.first * STA_DIST_STEP),
                       int((kv.first + 1) * STA_DIST_STEP),
                       kv.second.describe().c_str());
    }

    log += stringify("Cross-correlations by station\n");
    log += stringify("Station       #Phases GoodCC AvgCoeff(+/-) "
                     "GoodCC/Ph(+/-) time-diff[msec] (+/-)\n");
    for (const auto &kv : statsByStation)
    {
      log += stringify("%-12s %s\n", kv.first.c_str(),
                       kv.second.describe().c_str());
    }
    SEISCOMP_WARNING("%s", log.c_str());
  };

  resetCounters();
  int loop = 0;

  for (const auto &kv : _bgCat->getEvents())
  {
    const Event &event = kv.second;

    // find the neighbouring events
    NeighboursPtr neighbours;
    try
    {
      neighbours = selectNeighbouringEvents(
          _bgCat, event, _bgCat, clustOpt.minWeight, clustOpt.minESdist,
          clustOpt.maxESdist, clustOpt.minEStoIEratio, clustOpt.minDTperEvt,
          clustOpt.maxDTperEvt, clustOpt.minNumNeigh, clustOpt.maxNumNeigh,
          clustOpt.numEllipsoids, clustOpt.maxEllipsoidSize, false);
    }
    catch (...)
    {
      continue;
    }

    CatalogPtr catalog;

    if (theoretical)
    {
      // create theoretical phases for this event instead of fetching its
      // phases from the catalog
      catalog = neighbours->toCatalog(_bgCat, false);
      addMissingEventPhases(event, catalog, _bgCat, neighbours);
    }
    else
    {
      catalog = neighbours->toCatalog(_bgCat, true);
    }

    // Cross-correlate every neighbour phase with its corresponding event
    // theoretical phase.
    XCorrCache xcorr;
    buildXcorrDiffTTimePairs(catalog, neighbours, event,
                             clustOpt.xcorrMaxEvStaDist,
                             clustOpt.xcorrMaxInterEvDist, xcorr);

    // Update theoretical and automatic phase pick time and uncertainties based
    // on cross-correlation results. Drop theoretical phases wihout any good
    // cross-correlation result.
    if (theoretical)
    {
      fixPhases(catalog, event, xcorr);
    }

    //
    // Compare the detected phases with the actual event phases (manual or
    // automatic).
    //
    XCorrEvalStats evStats;

    for (const auto &kv : neighbours->allPhases())
      for (Phase::Type phaseType : kv.second)
      {
        //
        //  collect stats by event, station, station distance
        //
        const string stationId = kv.first;
        const Phase &catalogPhase =
            _bgCat->searchPhase(event.id, stationId, phaseType)->second;

        XCorrEvalStats phStaStats;
        double phaseTimeDiff = 0;

        if (xcorr.has(event.id, stationId, phaseType))
        {
          const auto &xentry = xcorr.get(event.id, stationId, phaseType);
          if (theoretical)
          {
            const Phase &detectedPhase =
                catalog->searchPhase(event.id, stationId, phaseType)->second;
            phaseTimeDiff = (catalogPhase.time - detectedPhase.time).length();
          }
          else
          {
            phaseTimeDiff = xentry.mean_lag;
          }
          phStaStats.addGoodCC(xentry.mean_coeff, xentry.ccCount,
                               phaseTimeDiff);
        }
        else
        {
          phStaStats.addBadCC();
        }

        evStats += phStaStats;
        totalStats += phStaStats;
        if (phaseType == Phase::Type::P) pPhaseStats += phStaStats;
        if (phaseType == Phase::Type::S) sPhaseStats += phStaStats;
        statsByStation[catalogPhase.stationId] += phStaStats;

        const Station &station =
            _bgCat->getStations().at(catalogPhase.stationId);
        double stationDistance = computeDistance(event, station);
        statsByStaDistance[int(stationDistance / STA_DIST_STEP)] += phStaStats;

        //
        //  collect stats by inter-event distance
        //
        map<unsigned, XCorrEvalStats> tmpStatsByInterEvDistance;

        for (unsigned neighEvId : neighbours->ids)
        {
          if (neighbours->has(neighEvId, stationId, phaseType))
          {
            const Event &neighbEv  = catalog->getEvents().at(neighEvId);
            double interEvDistance = computeDistance(event, neighbEv);
            XCorrEvalStats &interEvDistStats =
                tmpStatsByInterEvDistance[int(interEvDistance / EV_DIST_STEP)];
            interEvDistStats.total = 1;
            if (xcorr.has(event.id, neighEvId, stationId, phaseType))
            {
              const auto &xpi =
                  xcorr.get(event.id, neighEvId, stationId, phaseType);
              interEvDistStats.goodCC = 1;
              interEvDistStats.ccCount.push_back(1);
              interEvDistStats.ccCoeff.push_back(xpi.coeff);
              if (theoretical)
              {
                const auto &xentry = xcorr.get(event.id, stationId, phaseType);
                double timeDiff = phaseTimeDiff - (xentry.mean_lag - xpi.lag);
                interEvDistStats.timeDiff.push_back(timeDiff);
              }
              else
              {
                interEvDistStats.timeDiff.push_back(xpi.lag);
              }
            }
          }
        }

        for (auto &kv : tmpStatsByInterEvDistance)
        {
          const double interEvDistanceBucket = kv.first;
          XCorrEvalStats &newStats           = kv.second;
          if (newStats.goodCC > 0)
          {
            newStats.ccCount  = {std::accumulate(newStats.ccCount.begin(),
                                                newStats.ccCount.end(), 0.)};
            newStats.ccCoeff  = {computeMean(newStats.ccCoeff)};
            newStats.timeDiff = {computeMean(newStats.timeDiff)};
          }
          statsByInterEvDistance[interEvDistanceBucket] += newStats;
        }
      }

    SEISCOMP_WARNING("Event %-5s mag %3.1f %s", string(event).c_str(),
                     event.magnitude, evStats.describeShort().c_str());

    if (++loop % 100 == 0)
    {
      printStats("<PROGRESSIVE STATS>");
    }
  }

  printStats("<FINAL STATS>");
  printCounters();
}

} // namespace HDD
} // namespace Seiscomp

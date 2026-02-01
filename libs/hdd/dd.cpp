/***************************************************************************
 * MIT License                                                             *
 *                                                                         *
 * Copyright (C) by ETHZ/SED                                               *
 *                                                                         *
 * Permission is hereby granted, free of charge, to any person obtaining a *
 * copy of this software and associated documentation files (the           *
 * “Software”), to deal in the Software without restriction, including     *
 * without limitation the rights to use, copy, modify, merge, publish,     *
 * distribute, sublicense, and/or sell copies of the Software, and to      *
 * permit persons to whom the Software is furnished to do so, subject to   *
 * the following conditions:                                               *
 *                                                                         *
 * The above copyright notice and this permission notice shall be          *
 * included in all copies or substantial portions of the Software.         *
 *                                                                         *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,         *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      *
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  *
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY    *
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,    *
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE       *
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                  *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#include "dd.h"
#include "log.h"
#include "random.h"
#include "utils.h"
#include "xcorr.h"
#include <limits>

using namespace std;
using namespace HDD::Logger;
using std::chrono::duration;
using Event     = HDD::Catalog::Event;
using Phase     = HDD::Catalog::Phase;
using Station   = HDD::Catalog::Station;
using Transform = HDD::Waveform::Processor::Transform;
using HDD::Waveform::getBandAndInstrumentCodes;
using HDD::Waveform::getOrientationCode;

namespace {

struct WfCounters
{
  unsigned downloaded;
  unsigned no_avail;
  unsigned disk_cached;

  void update(HDD::Waveform::Loader *loader)
  {
    HDD::Waveform::BasicLoader *basic =
        dynamic_cast<HDD::Waveform::BasicLoader *>(loader);
    if (basic)
    {
      downloaded += basic->_counters_wf_downloaded;
      basic->_counters_wf_downloaded = 0;
      no_avail += basic->_counters_wf_no_avail;
      basic->_counters_wf_no_avail = 0;
      return;
    }

    HDD::Waveform::BatchLoader *batch =
        dynamic_cast<HDD::Waveform::BatchLoader *>(loader);
    if (batch)
    {
      downloaded += batch->_counters_wf_downloaded;
      batch->_counters_wf_downloaded = 0;
      no_avail += batch->_counters_wf_no_avail;
      batch->_counters_wf_no_avail = 0;
      return;
    }
  }

  void update(HDD::Waveform::DiskCachedLoader *diskCache)
  {
    if (diskCache)
    {
      disk_cached += diskCache->_counters_wf_cached;
      diskCache->_counters_wf_cached = 0;
    }
  }
};

string generateWorkingSubDirPath(const string &prefix)
{
  HDD::UTCTime now = HDD::UTCClock::now();
  int year, month, day, hour, min, sec, usec;
  HDD::UTCClock::toDate(now, year, month, day, hour, min, sec, usec);

  HDD::UniformRandomer ran(0, 9999);
  ran.setSeed(HDD::durToSec(now.time_since_epoch()));

  // preifx_creationTime_randomNumber
  string id = HDD::strf("%s_%04d%02d%02d%02d%02d%02d_%04zu", prefix.c_str(),
                        year, month, day, hour, min, sec, ran.next());
  return id;
}

void writeDoubleDifferenceToFile(
    const std::vector<HDD::Solver::DoubleDifference> &dds,
    const HDD::Catalog &cat,
    const std::string &file)
{
  ofstream os(file);
  os << "eventId1,eventId2,networkCode,stationCode,locationCode,phaseType,"
        "weight,xcorrUsed,xcorrCoefficient,doubleDifferenceResidual,"
        "interEventDistance"
     << endl;

  for (const auto &e : dds)
  {
    const HDD::Catalog::Station &sta = cat.getStations().at(e.staId);
    os << HDD::strf("%u,%u,%s,%s,%s,%s,%g,%s,%g,%g,%g", e.evId1, e.evId2,
                    sta.networkCode.c_str(), sta.stationCode.c_str(),
                    sta.locationCode.c_str(), e.phase.c_str(), e.weight,
                    e.xcorrUsed ? "true" : "false", e.xcorrCoeff,
                    e.doubleDifference, e.interEventDistance)
       << endl;
  }
}

} // namespace

namespace HDD {

DD::DD(const Catalog &catalog,
       const Config &cfg,
       std::unique_ptr<HDD::TravelTimeTable> ttt,
       std::unique_ptr<HDD::Waveform::Proxy> wf)
    : _cfg(cfg),
      _bgCat(Catalog::fillProcessingInfo(catalog,
                                         Phase::Source::CATALOG,
                                         _cfg.validPphases,
                                         _cfg.validSphases,
                                         _cfg.pickUncertaintyClasses)),
      _evTree(createEventTree(_bgCat)), _ttt(std::move(ttt)),
      _proxy(std::move(wf))
{
  disableCatalogWaveformDiskCache();
}

void DD::disableCatalogWaveformDiskCache()
{
  _wfAccess.useCache = false;
  initWaveformAccess();
}

void DD::enableCatalogWaveformDiskCache(const std::string &cacheDir,
                                        double diskTraceMinLen)
{
  _wfAccess.useCache        = true;
  _wfAccess.cacheDir        = cacheDir;
  _wfAccess.diskTraceMinLen = diskTraceMinLen;
  if (!pathExists(_wfAccess.cacheDir))
  {
    if (!createDirectories(_wfAccess.cacheDir))
    {
      throw Exception("Unable to create cache directory: " +
                      _wfAccess.cacheDir);
    }
  }
  initWaveformAccess();
}

void DD::initWaveformAccess()
{
  _wfAccess.loader = make_shared<Waveform::BasicLoader>(_proxy);
  _wfAccess.diskCache.reset();
  _wfAccess.extraLen.reset();
  _wfAccess.basicProc.reset();
  _wfAccess.memCache.reset();

  shared_ptr<Waveform::Loader> currLdr = _wfAccess.loader;

  if (_wfAccess.useCache)
  {
    _wfAccess.diskCache = make_shared<Waveform::DiskCachedLoader>(
        _proxy, currLdr, _wfAccess.cacheDir);
    _wfAccess.extraLen = make_shared<Waveform::ExtraLenLoader>(
        _wfAccess.diskCache, _wfAccess.diskTraceMinLen);
    currLdr = _wfAccess.extraLen;
  }

  _wfAccess.basicProc = make_shared<Waveform::BasicProcessor>(
      _proxy, currLdr, _cfg.wfFilter.extraTraceLen);
  _wfAccess.memCache =
      make_shared<Waveform::MemCachedProc>(_wfAccess.basicProc);
}

void DD::replaceWaveformLoader(const shared_ptr<Waveform::Loader> &baseLdr)
{
  if (_wfAccess.useCache)
  {
    _wfAccess.diskCache->setAuxLoader(baseLdr);
  }
  else
  {
    _wfAccess.basicProc->setAuxLoader(baseLdr);
  }
}

void DD::unloadWaveforms() { initWaveformAccess(); }

const std::vector<std::string> DD::xcorrComponents(const XcorrOptions &xcorrOpt,
                                                   const Phase &phase) const
{
  const auto &xcorrCfg = xcorrOpt.phase.at(phase.procInfo.type);
  if (xcorrCfg.components.empty() ||
      (xcorrCfg.components.size() == 1 && xcorrCfg.components.at(0).empty()))
  {
    return {getOrientationCode(phase.channelCode)};
  }
  else
  {
    return xcorrCfg.components;
  }
}

void DD::loadCatalogWaveformDiskCache(const XcorrOptions &xcorrOpt)
{
  //
  // Load waveforms and store them on disk.
  // For better performance we want to load waveforms in batch (batchloader)
  //
  if (!_wfAccess.useCache)
  {
    logInfoF("Waveform disk cache is disabled: nothing to do.");
    return;
  }

  logInfoF("Loading catalog waveform data (%zu events to load)",
           _bgCat.getEvents().size());

  unsigned numEvents = 0;
  WfCounters wfcount{};

  for (const auto &kv : _bgCat.getEvents())
  {
    const Event &event = kv.second;

    logDebugF("Loading event %u waveforms...", event.id);

    // Use a BatchLoader for the current event
    shared_ptr<Waveform::BatchLoader> batchLoader =
        make_shared<Waveform::BatchLoader>(_proxy);
    replaceWaveformLoader(batchLoader);

    // The following loop won't try to load the waveforms. Instead it will
    // records in BatchLoader the required waveforms that will be loaded
    // in batch later
    auto eqlrng = _bgCat.getPhases().equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      const Phase &phase = it->second;
      if (phase.procInfo.type == Phase::Type::NO)
      {
        continue;
      }
      TimeWindow tw = xcorrTimeWindowLong(xcorrOpt, event, phase);
      for (const string &component : xcorrComponents(xcorrOpt, phase))
      {
        getWaveform(*_wfAccess.basicProc, tw, event, phase, component);
      }
    }

    // This will actually download the waveworms
    batchLoader->load();

    // now  store the waveforms to disk
    eqlrng = _bgCat.getPhases().equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      const Phase &phase = it->second;
      if (phase.procInfo.type == Phase::Type::NO)
      {
        continue;
      }
      TimeWindow tw = xcorrTimeWindowLong(xcorrOpt, event, phase);
      for (const string &component : xcorrComponents(xcorrOpt, phase))
      {
        // _wfAccess.basicProc because we don't want to cache the
        // whole catalog waveforms into memory with _wfAccess.memCache
        getWaveform(*_wfAccess.basicProc, tw, event, phase, component);
      }
    }

    wfcount.update(batchLoader.get());

    const unsigned onePercent =
        _bgCat.getEvents().size() < 100 ? 1 : (_bgCat.getEvents().size() / 100);
    if (++numEvents % onePercent == 0)
    {
      logInfoF("Loaded %.1f%% of catalog phase waveforms",
               (numEvents * 100.0 / _bgCat.getEvents().size()));
    }
  }

  // Restore original loader
  replaceWaveformLoader(_wfAccess.loader);

  wfcount.update(_wfAccess.diskCache.get());

  logInfoF(
      "Finished loading catalog waveform data: total events %zu."
      "Waveforms downloaded %u, not available %u, loaded from disk cache %u",
      _bgCat.getEvents().size(), wfcount.downloaded, wfcount.no_avail,
      wfcount.disk_cached);
}

void DD::dumpWaveforms(const XcorrOptions &xcorrOpt, const string &basePath)
{
  if (!basePath.empty() && !pathExists(basePath))
  {
    if (!createDirectories(basePath))
    {
      throw Exception("Unable to create waveform dump directory: " + basePath);
    }
  }

  auto waveformPath = [](const string &wfDebugDir, const Catalog::Event &ev,
                         const Catalog::Phase &ph,
                         const string &channelCode) -> string {
    string debugFile =
        HDD::strf("ev%u.%s.%s.%s.%s.%s.mseed", ev.id, ph.networkCode.c_str(),
                  ph.stationCode.c_str(), ph.locationCode.c_str(),
                  channelCode.c_str(), ph.type.c_str());
    return HDD::joinPath(wfDebugDir, debugFile);
  };

  for (const auto &kv : _bgCat.getEvents())
  {
    const Event &event = kv.second;
    auto eqlrng        = _bgCat.getPhases().equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      const Phase &phase = it->second;

      if (phase.procInfo.type == Phase::Type::NO)
      {
        continue;
      }

      TimeWindow tw = xcorrTimeWindowShort(xcorrOpt, event, phase);

      for (const string &component : xcorrComponents(xcorrOpt, phase))
      {
        // _wfAccess.basicProc because we don't want to cache the
        // whole catalog waveforms into memory with _wfAccess.memCache
        shared_ptr<const Trace> tr =
            getWaveform(*_wfAccess.basicProc, tw, event, phase, component);

        if (!tr) continue;

        string channelCode =
            getBandAndInstrumentCodes(phase.channelCode) + component;

        const string wfPath = waveformPath(basePath, event, phase, channelCode);

        logInfoF("Writing %s", wfPath.c_str());
        try
        {
          _proxy->writeTrace(*tr, wfPath);
        }
        catch (exception &e)
        {
          logWarningF("Couldn't write waveform to disk %s: %s", wfPath.c_str(),
                      e.what());
        }
      }
    }
  }
}

std::list<unordered_map<unsigned, Neighbours>>
DD::findClusters(const ClusteringOptions &clustOpt)
{
  // find Neighbours for each event in the catalog
  unordered_map<unsigned, Neighbours> neighboursByEvent =
      selectNeighbouringEventsCatalog(
          _evTree, _bgCat, clustOpt.minEvStaDist, clustOpt.maxEvStaDist,
          clustOpt.minEvStaToInterEvRatio, clustOpt.minNumPhases,
          clustOpt.maxNumPhases, clustOpt.minNumNeigh, clustOpt.maxNumNeigh,
          clustOpt.numEllipsoids, clustOpt.maxNeighbourDist);

  // Organize the neighbours by not connected clusters
  return clusterizeNeighbouringEvents(std::move(neighboursByEvent));
}

Catalog DD::relocateMultiEvents(
    std::list<unordered_map<unsigned, Neighbours>> &clusters,
    XCorrCache &xcorrData,
    const ClusteringOptions &clustOpt,
    const XcorrOptions &xcorrOpt,
    const SolverOptions &solverOpt,
    bool saveProcessing,
    string processingDataDir)
{
  logInfoF(
      "Starting DD relocator in multiple events mode. The catalog contains "
      "%zu events",
      _bgCat.getEvents().size());

  string logFile;
  if (saveProcessing)
  {
    // prepare a folder for processing files
    if (processingDataDir.empty())
    {
      do
      {
        processingDataDir = generateWorkingSubDirPath("multievent");
      } while (pathExists(processingDataDir));
    }

    if (!pathExists(processingDataDir))
    {
      if (!createDirectories(processingDataDir))
      {
        throw Exception("Unable to create directory: " + processingDataDir);
      }
    }
    logInfoF("Saving processing data in %s", processingDataDir.c_str());

    // prepare file logger
    logFile = joinPath(processingDataDir, "relocation.log");
    addFileLogger(logFile, Level::info);
  }

  if (clusters.empty())
  {
    // find Neighbours for each event in the catalog
    unordered_map<unsigned, Neighbours> neighboursByEvent =
        selectNeighbouringEventsCatalog(
            _evTree, _bgCat, clustOpt.minEvStaDist, clustOpt.maxEvStaDist,
            clustOpt.minEvStaToInterEvRatio, clustOpt.minNumPhases,
            clustOpt.maxNumPhases, clustOpt.minNumNeigh, clustOpt.maxNumNeigh,
            clustOpt.numEllipsoids, clustOpt.maxNeighbourDist);

    // Organize the neighbours by non-connected clusters
    clusters = clusterizeNeighbouringEvents(std::move(neighboursByEvent));
  }

  logInfoF("Found %zu event clusters with the following number of events:",
           clusters.size());
  for (const auto &cluster : clusters)
  {
    logInfoF(" %zu events", cluster.size());
  }

  if (saveProcessing)
  {
    _bgCat.writeToFile(joinPath(processingDataDir, "input-event.csv"),
                       joinPath(processingDataDir, "input-phase.csv"),
                       joinPath(processingDataDir, "input-station.csv"));
  }

  //
  // relocate one cluster at the time
  //
  Catalog relocatedCatalog{};

  unsigned clusterId = 1;
  for (const auto &cluster : clusters)
  {
    logInfoF("Relocating cluster %u (%zu events)", clusterId, cluster.size());

    string clusterLogFile;
    if (saveProcessing)
    {
      string prefix  = strf("cluster-%u", clusterId);
      clusterLogFile = joinPath(processingDataDir, prefix + "-relocation.log");
      addFileLogger(clusterLogFile, Level::info);
      Neighbours::writeToFile(
          cluster, _bgCat, joinPath(processingDataDir, prefix + "-pair.csv"));
      Catalog catToDump;
      for (const auto &kv : cluster)
      {
        catToDump.add(kv.first, _bgCat, true);
      }
      catToDump.writeToFile(
          joinPath(processingDataDir, (prefix + "-event.csv")),
          joinPath(processingDataDir, (prefix + "-phase.csv")),
          joinPath(processingDataDir, (prefix + "-station.csv")));
    }

    // Perform cross-correlation
    XCorrCache xcorr =
        buildXCorrCache(_bgCat, cluster, xcorrOpt, true, xcorrData);

    // the actual relocation
    vector<Solver::DoubleDifference> startDDs, finalDDs;
    Catalog relocatedCluster = relocate(_bgCat, cluster, xcorrOpt, solverOpt,
                                        false, xcorr, startDDs, finalDDs);

    relocatedCatalog.add(relocatedCluster, true);

    if (saveProcessing)
    {
      string prefix = strf("cluster-%u", clusterId);
      writeDoubleDifferenceToFile(
          startDDs, _bgCat,
          joinPath(processingDataDir,
                   (prefix + "-initial-double-difference.csv")));
      writeDoubleDifferenceToFile(
          finalDDs, _bgCat,
          joinPath(processingDataDir,
                   (prefix + "-final-double-difference.csv")));
      relocatedCluster.writeToFile(
          joinPath(processingDataDir, (prefix + "-reloc-event.csv")),
          joinPath(processingDataDir, (prefix + "-reloc-phase.csv")));
      removeFileLogger(clusterLogFile);
    }
    clusterId++;

    // make sure xcorrData stores all the cross-correlation results
    xcorrData.add(xcorr);
  }

  if (saveProcessing)
  {
    relocatedCatalog.writeToFile(
        joinPath(processingDataDir, "reloc-event.csv"),
        joinPath(processingDataDir, "reloc-phase.csv"));
    xcorrData.writeToFile(_bgCat, joinPath(processingDataDir, "xcorr.csv"));
    removeFileLogger(logFile);
  }

  return relocatedCatalog;
}

Catalog DD::relocateSingleEvent(const Catalog &singleEvent,
                                bool isManual,
                                const ClusteringOptions &clustOpt,
                                const XcorrOptions &xcorrOpt,
                                const SolverOptions &solverOpt,
                                bool saveProcessing,
                                string processingDataDir)
{
  if (saveProcessing)
  {
    // prepare a folder for processing files
    if (processingDataDir.empty())
    {
      do
      {
        processingDataDir = generateWorkingSubDirPath("singleevent");
      } while (pathExists(processingDataDir));
    }

    if (!pathExists(processingDataDir))
    {
      if (!createDirectories(processingDataDir))
      {
        throw Exception("Unable to create directory: " + processingDataDir);
      }
    }
    logInfoF("Saving processing data in %s", processingDataDir.c_str());

    singleEvent.writeToFile(
        joinPath(processingDataDir, "single-event.csv"),
        joinPath(processingDataDir, "single-event-phase.csv"),
        joinPath(processingDataDir, "single-event-station.csv"));
  }

  // prepare file logger
  string logFile = joinPath(processingDataDir, "relocation.log");
  if (saveProcessing)
  {
    addFileLogger(logFile, Level::info);
  }

  Catalog evToRelocateCat = Catalog::fillProcessingInfo(
      singleEvent,
      (isManual ? Phase::Source::RT_EVENT_MANUAL
                : Phase::Source::RT_EVENT_AUTOMATIC),
      _cfg.validPphases, _cfg.validSphases, _cfg.pickUncertaintyClasses);

  // extract event to relocate
  const Event &evToRelocate = evToRelocateCat.getEvents().begin()->second;

  //
  // select neighbouring events
  //
  Neighbours neighbours = selectNeighbouringEvents(
      _evTree, _bgCat, evToRelocate, evToRelocateCat, clustOpt.minEvStaDist,
      clustOpt.maxEvStaDist, clustOpt.minEvStaToInterEvRatio,
      clustOpt.minNumPhases, clustOpt.maxNumPhases, clustOpt.maxNumNeigh,
      clustOpt.numEllipsoids, clustOpt.maxNeighbourDist);

  logInfoF("Found %zu neighbouring events", neighbours.amount());

  // check if enough neighbors were found
  if (neighbours.amount() < clustOpt.minNumNeigh)
  {
    logInfo("Insufficient number of neighbors " + neighbours.amount());
    return Catalog{};
  }

  //
  // prepare catalog to relocate
  //
  Catalog catalog = neighbours.toCatalog(_bgCat);
  unsigned evToRelocateNewId =
      catalog.add(evToRelocate.id, evToRelocateCat, false);
  neighbours.setReferenceId(evToRelocateNewId);

  if (saveProcessing)
  {
    catalog.writeToFile(joinPath(processingDataDir, "input-event.csv"),
                        joinPath(processingDataDir, "input-phase.csv"),
                        joinPath(processingDataDir, "input-station.csv"));
  }

  unordered_map<unsigned, Neighbours> cluster;
  cluster.emplace(neighbours.referenceId(), std::move(neighbours));

  if (saveProcessing)
  {
    Neighbours::writeToFile(cluster, catalog,
                            joinPath(processingDataDir, "pair.csv"));
  }

  // Perform cross-correlation
  XCorrCache xcorr = buildXCorrCache(catalog, cluster, xcorrOpt, false);

  // the actual relocation
  vector<Solver::DoubleDifference> startDDs, finalDDs;
  Catalog relocatedEvCat = relocate(catalog, cluster, xcorrOpt, solverOpt, true,
                                    xcorr, startDDs, finalDDs);

  if (saveProcessing)
  {
    writeDoubleDifferenceToFile(
        startDDs, catalog,
        joinPath(processingDataDir, "initial-double-difference.csv"));
    writeDoubleDifferenceToFile(
        finalDDs, catalog,
        joinPath(processingDataDir, "final-double-difference.csv"));
    relocatedEvCat.writeToFile(joinPath(processingDataDir, "reloc-event.csv"),
                               joinPath(processingDataDir, "reloc-phase.csv"));
  }

  removeFileLogger(logFile);

  return relocatedEvCat;
}

Catalog DD::relocate(const Catalog &catalog,
                     const unordered_map<unsigned, Neighbours> &cluster,
                     const XcorrOptions &xcorrOpt,
                     const SolverOptions &solverOpt,
                     bool keepNeighboursFixed,
                     const XCorrCache &xcorr,
                     std::vector<Solver::DoubleDifference> &startDDs,
                     std::vector<Solver::DoubleDifference> &finalDDs) const
{
  logInfo("Building and solving double-difference system...");

  //
  // iterate the solver computation multiple times
  //
  Catalog currentCatalog(catalog), prevCatalog{};

  for (unsigned iteration = 0; iteration <= solverOpt.algoIterations;
       iteration++)
  {
    bool isFirstIteration = (iteration == 0);
    bool isLastIteration  = (iteration == solverOpt.algoIterations);

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
    double absLocConstraint = interpolate(solverOpt.absLocConstraintStart,
                                          solverOpt.absLocConstraintEnd);

    logInfoF("Iteration %u num events %zu. Parameters: dampingFactor=%.2f "
             "downWeightingByResidual=%.2f absLocConstraint=%.2f",
             iteration, cluster.size(), dampingFactor, downWeightingByResidual,
             absLocConstraint);

    //
    // create a solver and then add observations
    //
    Solver solver(solverOpt.type);

    for (const auto &kv : cluster)
    {
      unsigned eventId = kv.first;
      bool tttOk =
          addObservations(solver, currentCatalog, kv.second,
                          keepNeighboursFixed, solverOpt.usePickUncertainties,
                          solverOpt.xcorrWeightScaler, xcorrOpt, xcorr);

      if (!tttOk && !prevCatalog.empty() &&
          currentCatalog.getEvents().at(eventId).relocInfo.isRelocated)
      {
        // No trave time information: the event was reloacated but the new
        // location is outside the travel time table boundaries (e.g. an
        // air-quake when the ttt do not support negative depths).
        // Restore the previous location of the event
        logInfoF("Event %u relocated, but it is now outside the travel time "
                 "table boundaries. Revert it to the previous location.",
                 eventId);
        currentCatalog.removeEvent(eventId);
        currentCatalog.add(eventId, prevCatalog, true);
        addObservations(solver, currentCatalog, kv.second, keepNeighboursFixed,
                        solverOpt.usePickUncertainties,
                        solverOpt.xcorrWeightScaler, xcorrOpt, xcorr);
      }
    }

    //
    // Compute and log event rms
    //
    currentCatalog = computeEventResiduals(solver, currentCatalog, solverOpt,
                                           cluster, isFirstIteration);

    //
    // Prepare the solver, which logs the DD residuals
    //
    solver.prepare(absLocConstraint, downWeightingByResidual);

    if (isFirstIteration)
    {
      startDDs = solver.getDoubleDifferences();
    }

    if (isLastIteration)
    {
      // Last iteration is just fo logging the residuals and collect latest
      // double differences. So exit without solving the system
      finalDDs = solver.getDoubleDifferences();
      break;
    }

    //
    // solve the system
    //
    logInfoF("Solving...");
    try
    {
      solver.solve(solverOpt.solverIterations, dampingFactor,
                   solverOpt.L2normalization);
    }
    catch (Exception &e)
    {
      logInfoF("Cannot solve the double-difference system, stop here (%s)",
               e.what());
      break;
    }

    //
    // update event location/time
    //
    prevCatalog    = std::move(currentCatalog);
    currentCatalog = updateRelocatedEvents(solver, prevCatalog, solverOpt,
                                           cluster, isFirstIteration);
  }

  //
  // Build and return the catalog with only relocated events/phases
  //
  Catalog catalogToReturn{};
  for (const auto &kv : currentCatalog.getEvents())
  {
    const Event &event = kv.second;
    if (event.relocInfo.isRelocated)
    {
      catalogToReturn.add(event.id, currentCatalog, true);
    }
  }
  return catalogToReturn;
}

/*
 * Add both the absolute travel time differences and the differential
 * travel times from the cross-correlation for pairs of earthquakes to the
 * solver.
 */
bool DD::addObservations(Solver &solver,
                         const Catalog &catalog,
                         const Neighbours &neighbours,
                         bool keepNeighboursFixed,
                         bool usePickUncertainties,
                         double xcorrWeightScaler,
                         const XcorrOptions &xcorrOpt,
                         const XCorrCache &xcorr) const
{
  const Event &refEv = catalog.getEvents().at(neighbours.referenceId());

  // Detect an event that moved to a location outside the ttt range
  unsigned tttAttempted = 0, tttAvailable = 0;
  string badStations;

  // stationId, phase, neighbourId
  vector<tuple<string, string, unsigned>> nPhases = neighbours.phases();

  // sort the nPhases so that we can search the reference event phase only once
  // and reuse it for multiple neighbours (see loop)
  std::sort(nPhases.begin(), nPhases.end());

  //
  // Loop through neighbours events that matched refEv phases
  //
  const Phase *refPhase  = nullptr;
  const Station *station = nullptr;
  bool tttError;
  for (const auto &t : nPhases)
  {
    const string stationId   = std::get<0>(t);
    const string phaseType   = std::get<1>(t);
    const unsigned neighEvId = std::get<2>(t);

    if (!station || station->id != stationId)
    {
      station = &catalog.getStations().at(stationId);
    }

    if (!refPhase || refPhase->stationId != stationId ||
        refPhase->type != phaseType)
    {
      refPhase = &catalog.searchPhase(refEv.id, stationId, phaseType)->second;
      tttError = false;

      // Query Travel Time Table Reference Phase
      tttAttempted++;
      if (!addObservationParams(solver, *_ttt, refEv, *station, *refPhase,
                                true))
      {
        tttError = true;
        badStations += " " + phaseType + "@" + stationId;
      }
      else
      {
        tttAvailable++;
      }
    }

    if (tttError)
    {
      logDebugF("Skipping observation (ev %u-%u sta %s phase %s)", refEv.id,
                neighEvId, station->id.c_str(), phaseType.c_str());
      continue;
    }

    const auto &xcorrCfg = xcorrOpt.phase.at(refPhase->procInfo.type);

    const Event &event = catalog.getEvents().at(neighEvId);

    const Phase &phase =
        catalog.searchPhase(event.id, stationId, phaseType)->second;

    if (!addObservationParams(solver, *_ttt, event, *station, phase,
                              !keepNeighboursFixed))
    {
      logDebugF("Skipping observation (ev %u-%u sta %s phase %s)", refEv.id,
                event.id, station->id.c_str(), phaseType.c_str());
      continue;
    }

    //
    // compute travel times for both event and `refEv`
    //
    UTCTime::duration ref_travel_time = refPhase->time - refEv.time;
    if (ref_travel_time.count() < 0)
    {
      logDebugF("Ignoring phase %s with negative travel time ",
                string(*refPhase).c_str());
      continue;
    }

    UTCTime::duration travel_time = phase.time - event.time;
    if (travel_time.count() < 0)
    {
      logDebugF("Ignoring phase %s with negative travel time ",
                string(phase).c_str());
      continue;
    }

    //
    // compute absolute travel time differences to the solver
    //
    duration<double> diffTime = ref_travel_time - travel_time;
    double weight =
        usePickUncertainties
            ? ((refPhase->procInfo.classWeight + phase.procInfo.classWeight) /
               2.0)
            : 1.0;

    //
    // Check if we have cross-correlation results for the current
    // event/refEv pair at station/phase and refine the differential
    // time and update the weight
    //
    bool xcorrUsed    = false;
    double xcorrCoeff = 0;

    if (!xcorr.empty()) // xcorr is enabled
    {

      if (xcorr.has(refEv.id, event.id, stationId, phaseType))
      {
        const auto &e = xcorr.get(refEv.id, event.id, stationId, phaseType);
        if (e.valid)
        {
          xcorrCoeff = e.coeff;
          if (e.coeff >= xcorrCfg.minCoef)
          {
            xcorrUsed = true;
            diffTime -= duration<double>(e.lag);
            double weightGain = (weight * xcorrWeightScaler) - weight;
            weight += weightGain * e.coeff;
          }
        }
      }
    }

    solver.addObservation(refEv.id, event.id, stationId, phaseType,
                          diffTime.count(), weight, xcorrUsed, xcorrCoeff);
  }

  if (tttAttempted > 0 && tttAttempted != tttAvailable)
  {
    logWarningF("TTT: %u/%u queries failed for Event %u at lat %.6f "
                "lon %.6f depth %.6f: %s",
                tttAttempted - tttAvailable, tttAttempted, refEv.id,
                refEv.latitude, refEv.longitude, refEv.depth,
                badStations.c_str());
  }

  return tttAttempted == 0 || tttAvailable > 0;
}

bool DD::addObservationParams(Solver &solver,
                              TravelTimeTable &ttt,
                              const Catalog::Event &event,
                              const Catalog::Station &station,
                              const Catalog::Phase &phase,
                              bool computeEvChanges) const
{
  bool dummy1;
  double dummy2;
  //
  // check if this data has been already added to the solver
  //
  if (!solver.getObservationParams(event.id, station.id, phase.type, dummy1,
                                   dummy2, dummy2, dummy2, dummy2, dummy2))
  {
    //
    // Compute travel time information
    //
    double travelTime, takeOffAngleAzim, takeOffAngleDip, dtdd, dtdh;
    try
    {
      string phaseName = phase.type;
      if (_cfg.PSTableOnly)
      {
        if (phase.procInfo.type == Phase::Type::P)
        {
          phaseName = "P";
        }
        else if (phase.procInfo.type == Phase::Type::S)
        {
          phaseName = "S";
        }
      }
      ttt.compute(event, station, phaseName, travelTime, takeOffAngleAzim,
                  takeOffAngleDip, dtdd, dtdh);
    }
    catch (Exception &e)
    {
      logDebugF("TTT error: %s (Event lat %.6f lon %.6f depth %.6f "
                "Station %s)",
                e.what(), event.latitude, event.longitude, event.depth,
                station.id.c_str());
      return false;
    }
    double travelTimeResidual = travelTime - durToSec(phase.time - event.time);

    //
    // partial derivatives
    //
    dtdd               = dtdd / kmOfDegree(); // [sec/deg] -> [sec/km]
    const double angle = degToRad(90 - takeOffAngleAzim);

    const double dx = dtdd * -std::cos(angle);
    const double dy = dtdd * -std::sin(angle);
    const double dz = -dtdh;

    //
    // Populate the solver
    //
    solver.addEvent(event.id, event.latitude, event.longitude, event.depth);
    solver.addStation(station.id, station.latitude, station.longitude,
                      station.elevation);
    solver.addObservationParams(event.id, station.id, phase.type,
                                computeEvChanges, travelTime,
                                travelTimeResidual, dx, dy, dz);
  }
  return true;
}

Catalog
DD::computeEventResiduals(const Solver &solver,
                          const Catalog &catalog,
                          const SolverOptions &solverOpt,
                          const unordered_map<unsigned, Neighbours> &cluster,
                          bool isFirstIteration) const
{
  unordered_map<string, Station> stations    = catalog.getStations();
  map<unsigned, Event> events                = catalog.getEvents();
  unordered_multimap<unsigned, Phase> phases = catalog.getPhases();
  vector<double> allRms;

  //
  // loop through each event in the cluster
  //
  for (const auto &kv : cluster)
  {
    Event &event = events.at(kv.second.referenceId());

    if (isFirstIteration)
    {
      event.relocInfo.startRms = 0;
    }
    event.relocInfo.finalRms = 0;

    double sumSquaredResiduals = 0.0;
    double sumSquaredWeights   = 0.0;
    auto eqlrng                = phases.equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      Phase &phase           = it->second;
      const Station &station = stations.at(phase.stationId);

      if (isFirstIteration)
      {
        phase.relocInfo.startResidual = 0;
      }

      phase.relocInfo.finalResidual = 0;
      phase.relocInfo.weight =
          solverOpt.usePickUncertainties ? phase.procInfo.classWeight : 1;

      bool dummy1;
      double residual, dummy2;

      if (solver.getObservationParams(event.id, station.id, phase.type, dummy1,
                                      dummy2, residual, dummy2, dummy2, dummy2))
      {
        double weight = phase.relocInfo.weight;

        if (isFirstIteration)
        {
          phase.relocInfo.startResidual = residual;
        }
        else
        {
          phase.relocInfo.finalResidual = residual;
        }

        sumSquaredResiduals += (residual * weight) * (residual * weight);
        sumSquaredWeights += weight * weight;
      }
    }

    if (sumSquaredWeights > 0)
    {
      double rms = std::sqrt(sumSquaredResiduals / sumSquaredWeights);
      if (isFirstIteration)
      {
        event.relocInfo.startRms = rms;
      }
      else
      {
        event.relocInfo.finalRms = rms;
      }
      allRms.push_back(rms);
    }
  }

  double min, max, q1, q2, q3;
  compute5numberSummary(allRms, min, max, q1, q2, q3);
  logInfoF("Event RMS [msec]: min %.3f 1st quartile %.3f median %.3f 3rd "
           "quartile %.3f max %.3f ",
           min * 1000, q1 * 1000, q2 * 1000, q3 * 1000, max * 1000);

  return Catalog(stations, events, phases);
}

Catalog
DD::updateRelocatedEvents(const Solver &solver,
                          const Catalog &catalog,
                          const SolverOptions &solverOpt,
                          const unordered_map<unsigned, Neighbours> &cluster,
                          bool isFirstIteration) const
{
  unordered_map<string, Station> stations    = catalog.getStations();
  map<unsigned, Event> events                = catalog.getEvents();
  unordered_multimap<unsigned, Phase> phases = catalog.getPhases();
  unsigned relocatedEvs                      = 0;

  //
  // loop through each event in the cluster
  //
  for (const auto &kv : cluster)
  {
    Event &event = events.at(kv.second.referenceId());

    if (isFirstIteration)
    {
      event.relocInfo.isRelocated = false;
    }

    //
    // get relocation changes (km and sec) computed by the solver for
    // the current event
    //
    double deltaX, deltaY, deltaZ, deltaTime;
    if (!solver.getEventChanges(event.id, deltaX, deltaY, deltaZ, deltaTime))
    {
      continue;
    }

    double newDepth = event.depth - deltaZ; // depth opposite sign of Z
    UTCTime newTime = event.time + secToDur(deltaTime);

    //
    // convert delta lat lon from km to degree using computeCoordinates()
    //
    double distance = std::sqrt(square(deltaX) + square(deltaY));
    // Since azimuth=PI/2-angle we use the identity atan2(x,y)=PI/2-atan2(y,x)
    double azimuth = std::atan2(deltaX, deltaY);

    // Computes the coordinates (lat, lon) of the point which is at
    // azimuth [rad] and distance [km] as seen from the original event location
    double newLat, newLon;
    computeCoordinates(distance, azimuth, event.latitude, event.longitude,
                       newLat, newLon, event.depth);

    //
    // update event location/time and compute statistics
    //
    event.latitude              = newLat;
    event.longitude             = newLon;
    event.depth                 = newDepth;
    event.time                  = newTime;
    event.relocInfo.isRelocated = true;

    auto eqlrng = phases.equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      Phase &phase           = it->second;
      const Station &station = stations.at(phase.stationId);

      if (isFirstIteration)
      {
        phase.relocInfo.isRelocated = false;
      }

      if (!solver.isEventPhaseUsed(event.id, station.id, phase.type))
      {
        continue;
      }

      phase.relocInfo.isRelocated = true;
    }

    relocatedEvs++;
  }

  logInfoF("Successfully relocated %u events", relocatedEvs);

  return Catalog(stations, events, phases);
}

XCorrCache
DD::buildXCorrCache(const Catalog &catalog,
                    const unordered_map<unsigned, Neighbours> &cluster,
                    const XcorrOptions &xcorrOpt,
                    bool cacheRefEvWf,
                    const XCorrCache &precomputed)
{
  XCorrCache xcorr{};

  if (!xcorrOpt.enable)
  {
    logInfo("Cross-correlation is disabled");
    return xcorr;
  }

  logInfo("Computing differential times via cross-correlation...");

  unsigned long performed = 0;

  for (const auto &kv : cluster)
  {
    const Neighbours &neighbours = kv.second;
    const Event &refEv = catalog.getEvents().at(neighbours.referenceId());

    buildXcorrDiffTTimePairs(catalog, neighbours, refEv, cacheRefEvWf, xcorrOpt,
                             precomputed, xcorr);

    const unsigned onePerThousand =
        cluster.size() < 1000 ? 1 : (cluster.size() / 1000);
    if (++performed % onePerThousand == 0)
    {
      logInfoF("Cross-correlation completion %.1f%%",
               (performed * 100 / (double)cluster.size()));
    }
  }

  WfCounters wfcount{};
  wfcount.update(_wfAccess.loader.get());
  wfcount.update(_wfAccess.diskCache.get());

  logInfoF("Catalog waveform data: waveforms downloaded %u, not available "
           "%u, loaded from disk cache %u",
           wfcount.downloaded, wfcount.no_avail, wfcount.disk_cached);

  logXCorrSummary(cluster, xcorrOpt, xcorr);

  return xcorr;
}

/*
 * Compute and store to `XCorrCache` cross-correlated differential travel
 * times for pairs of the earthquake.
 */
void DD::buildXcorrDiffTTimePairs(const Catalog &catalog,
                                  const Neighbours &neighbours,
                                  const Event &refEv,
                                  bool cacheRefEvWf,
                                  const XcorrOptions &xcorrOpt,
                                  const XCorrCache &precomputed,
                                  XCorrCache &xcorr)
{
  logDebugF(
      "Computing cross-correlation differential travel times for event %s",
      string(refEv).c_str());
  //
  // Prepare the waveform loaders for single-event. Since its waveforms are
  // temporary and discarded after the relocation, we don't want to
  // make use of the background catalog cache (_wfAccess.memCache)
  //
  shared_ptr<Waveform::Processor> tmpMemCache;
  if (!cacheRefEvWf)
  {
    tmpMemCache =
        preloadNonCatalogWaveforms(catalog, neighbours, refEv, xcorrOpt);
  }

  //
  // loop through reference event phases
  //
  auto eqlrngRef = catalog.getPhases().equal_range(refEv.id);
  for (auto itRef = eqlrngRef.first; itRef != eqlrngRef.second; ++itRef)
  {
    const Phase &refPhase  = itRef->second;
    const Station &station = catalog.getStations().at(refPhase.stationId);

    //
    // skip stations too close or too far away
    //
    double stationDistance = computeDistance(refEv, station);
    if (stationDistance < xcorrOpt.minEvStaDist)
    {
      continue;
    }
    if (stationDistance > xcorrOpt.maxEvStaDist && xcorrOpt.maxEvStaDist >= 0)
    {
      continue;
    }

    //
    // select appropriate waveform loader for refPhase
    //
    shared_ptr<Waveform::Processor> refPrc =
        cacheRefEvWf ? _wfAccess.memCache : tmpMemCache;

    //
    // loop through neighbouring events and cross-correlate with `refPhase`
    //
    for (unsigned neighEvId : neighbours.ids())
    {
      const Event &event = catalog.getEvents().at(neighEvId);

      //
      // skip events too far away
      //
      double interEventDistance = computeDistance(refEv, event);
      if (interEventDistance > xcorrOpt.maxInterEvDist &&
          xcorrOpt.maxInterEvDist >= 0)
      {
        continue;
      }

      //
      // Check if this event pair has an entry for this station/phase
      //
      if (neighbours.has(neighEvId, refPhase.stationId, refPhase.type))
      {
        const Phase &phase =
            catalog.searchPhase(event.id, refPhase.stationId, refPhase.type)
                ->second;

        // In single-event mode `refPhase` is real-time and `phase` is from
        // the catalog. In multi-event mode both are from the catalog.
        if (phase.procInfo.source != Phase::Source::CATALOG)
        {
          throw Exception("Internal logic error: phase is not from catalog");
        }

        // skip cross-correlation if we already have the result in cache
        if (!xcorr.has(refEv.id, event.id, refPhase.stationId, refPhase.type))
        {
          try
          {
            // try to fetch the cross-correlation results from the precomputed
            // values
            const XCorrCache::Entry &e = precomputed.get(
                refEv.id, event.id, refPhase.stationId, refPhase.type);
            xcorr.add(refEv.id, event.id, refPhase.stationId, refPhase.type,
                      e.valid, e.coeff, e.lag);
          }
          catch (std::out_of_range &e)
          {
            // Do the actual cross-correlation if not found in precomputed
            double coeff, lag;
            string component;
            bool valid =
                xcorrPhases(xcorrOpt, refEv, refPhase, *refPrc, event, phase,
                            *_wfAccess.memCache, coeff, lag, component);

            // cache cross-correlation results
            xcorr.add(refEv.id, event.id, refPhase.stationId, refPhase.type,
                      valid, coeff, lag);
          }
        }
      }
    }
  }
}

shared_ptr<Waveform::Processor>
DD::preloadNonCatalogWaveforms(const Catalog &catalog,
                               const Neighbours &neighbours,
                               const Event &refEv,
                               const XcorrOptions &xcorrOpt)
{
  //
  // For single-event relocation we want to load the waveforms in batch
  // otherwise the seedlink server gets stuck and becomes unresponsive due to
  // the multiple connections requests, one for each event phase
  //
  shared_ptr<Waveform::BatchLoader> batchLoader(
      new Waveform::BatchLoader(_proxy));
  shared_ptr<Waveform::Processor> proc(new Waveform::BasicProcessor(
      _proxy, batchLoader, _cfg.wfFilter.extraTraceLen));

  //
  // loop through reference event phases
  //
  auto eqlrngRef = catalog.getPhases().equal_range(refEv.id);
  for (auto itRef = eqlrngRef.first; itRef != eqlrngRef.second; ++itRef)
  {
    const Phase &refPhase  = itRef->second;
    const Station &station = catalog.getStations().at(refPhase.stationId);
    const auto components  = xcorrComponents(xcorrOpt, refPhase);

    //
    // skip stations too close or too far away
    //
    double stationDistance = computeDistance(refEv, station);
    if (stationDistance < xcorrOpt.minEvStaDist)
    {
      continue;
    }
    if (stationDistance > xcorrOpt.maxEvStaDist && xcorrOpt.maxEvStaDist >= 0)
    {
      continue;
    }

    //
    // loop through neighbouring events and cross-correlate with `refPhase`
    //
    for (unsigned neighEvId : neighbours.ids())
    {
      const Event &event = catalog.getEvents().at(neighEvId);

      //
      // skip events too far away
      //
      double interEventDistance = computeDistance(refEv, event);
      if (interEventDistance > xcorrOpt.maxInterEvDist &&
          xcorrOpt.maxInterEvDist >= 0)
        continue;

      //
      // Check if this event pair has an entry for this station/phase
      //
      if (neighbours.has(neighEvId, refPhase.stationId, refPhase.type))
      {
        //
        // For each match load the reference event phase waveforms
        //
        TimeWindow tw = xcorrTimeWindowLong(xcorrOpt, refEv, refPhase);

        for (const string &component : components)
        {
          // This doesn't really load the trace but force the request to reach
          // the loader, which will load all the traces later
          getWaveform(*proc, tw, refEv, refPhase, component);
        }
      }
    }
  }

  // This will actually dowanload the waveforms all at once
  batchLoader->load();

  // print counters
  WfCounters wfcount{};
  wfcount.update(batchLoader.get());
  logInfoF("Event %s: waveforms downloaded %u, not available %u",
           string(refEv).c_str(), wfcount.downloaded, wfcount.no_avail);

  return make_shared<Waveform::MemCachedProc>(proc);
}

TimeWindow DD::xcorrTimeWindowLong(const XcorrOptions &xcorrOpt,
                                   const Event &event,
                                   const Phase &phase) const
{
  const auto &xcorrCfg = xcorrOpt.phase.at(phase.procInfo.type);
  TimeWindow tw        = xcorrTimeWindowShort(xcorrOpt, event, phase);
  tw.setStartTime(tw.startTime() - secToDur(xcorrCfg.maxDelay));
  tw.setEndTime(tw.endTime() + secToDur(xcorrCfg.maxDelay));
  return tw;
}

TimeWindow DD::xcorrTimeWindowShort(const XcorrOptions &xcorrOpt,
                                    const Event &event,
                                    const Phase &phase) const
{
  const auto &xcorrCfg = xcorrOpt.phase.at(phase.procInfo.type);
  const UTCTime::duration win_len =
      secToDur(xcorrCfg.endOffset - xcorrCfg.startOffset);
  const UTCTime::duration travel_time  = phase.time - event.time;
  const duration<double> win_extention = travel_time * xcorrCfg.winScaling;
  const double win_scaler              = (win_len + win_extention) / win_len;
  const duration<double> win_start =
      secToDur(xcorrCfg.startOffset) * win_scaler;
  const UTCTime final_start =
      phase.time + std::chrono::duration_cast<UTCTime::duration>(win_start);
  const UTCTime::duration final_win_len =
      std::chrono::duration_cast<UTCTime::duration>(win_len + win_extention);
  return TimeWindow(final_start, final_win_len);
}

bool DD::xcorrPhases(const XcorrOptions &xcorrOpt,
                     const Event &event1,
                     const Phase &phase1,
                     Waveform::Processor &ph1Cache,
                     const Event &event2,
                     const Phase &phase2,
                     Waveform::Processor &ph2Cache,
                     double &coeffOut,
                     double &lagOut,
                     string &componentOut)
{
  if (phase1.type != phase2.type)
  {
    logErrorF("Skipping cross-correlation: mismatching phases (%s and %s)",
              string(phase1).c_str(), string(phase2).c_str());
    return false;
  }

  const string channelCodeRoot1 = getBandAndInstrumentCodes(phase1.channelCode);
  const string channelCodeRoot2 = getBandAndInstrumentCodes(phase2.channelCode);

  if (channelCodeRoot1 != channelCodeRoot2)
  {
    logDebugF("Skipping cross-correlation: incompatible channels %s and %s "
              "(%s and %s)",
              channelCodeRoot1.c_str(), channelCodeRoot2.c_str(),
              string(phase1).c_str(), string(phase2).c_str());
    return false;
  }

  //
  // Try all registered components until we are able to perform
  // the cross-correlation
  //
  bool performed = false;
  for (const string &component : xcorrComponents(xcorrOpt, phase1))
  {
    performed =
        xcorrPhasesOneComponent(xcorrOpt, event1, phase1, ph1Cache, event2,
                                phase2, ph2Cache, component, coeffOut, lagOut);
    if (performed)
    {
      componentOut = component;
      break;
    }
  }
  return performed;
}

bool DD::xcorrPhasesOneComponent(const XcorrOptions &xcorrOpt,
                                 const Event &event1,
                                 const Phase &phase1,
                                 Waveform::Processor &ph1Cache,
                                 const Event &event2,
                                 const Phase &phase2,
                                 Waveform::Processor &ph2Cache,
                                 const string &component,
                                 double &coeffOut,
                                 double &lagOut)
{
  coeffOut = lagOut = 0;

  const auto &xcorrCfg = xcorrOpt.phase.at(phase1.procInfo.type);

  TimeWindow tw1 = xcorrTimeWindowLong(xcorrOpt, event1, phase1);
  TimeWindow tw2 = xcorrTimeWindowLong(xcorrOpt, event2, phase2);

  // Load the long `tr1`, because we want to cache the long version. Then
  // we'll trim it.
  shared_ptr<const Trace> tr1 =
      getWaveform(ph1Cache, tw1, event1, phase1, component);
  if (!tr1)
  {
    logDebugF("Skipping cross-correlation: no waveform data for phase %s",
              string(phase1).c_str());
    return false;
  }

  // Load the long `tr2`, because we want to cache the long version. Then
  // we'll trim it.
  shared_ptr<const Trace> tr2 =
      getWaveform(ph2Cache, tw2, event2, phase2, component);
  if (!tr2)
  {
    logDebugF("Skipping cross-correlation: no waveform data for phase %s",
              string(phase2).c_str());
    return false;
  }

  if (tr1->samplingFrequency() != tr2->samplingFrequency())
  {
    logWarningF(
        "Skipping cross-correlation: traces have different sampling freq "
        "(%f!=%f): phase1 %s and phase 2 %s",
        tr1->samplingFrequency(), tr2->samplingFrequency(),
        string(phase1).c_str(), string(phase2).c_str());
    return false;
  }

  // Trust the manual pick on `phase2`: keep `tr2` short and cross-correlate
  // it with the larger `tr1` window.
  double xcorr_coeff = numeric_limits<double>::quiet_NaN(), xcorr_lag = 0;

  if (phase2.trusted() || (!phase1.trusted() && !phase2.trusted()))
  {
    // Trim `tr2` to shorter length; we want to cross-correlate the short one
    // with the long one.
    Trace tr2Short(*tr2);
    TimeWindow tw2Short = xcorrTimeWindowShort(xcorrOpt, event2, phase2);
    if (!tr2Short.slice(tw2Short))
    {
      logDebugF("Skipping cross-correlation: cannot trim phase2 waveform "
                "for phase pair phase1='%s', phase2='%s'",
                string(phase1).c_str(), string(phase2).c_str());
      return false;
    }

    xcorr(*tr1, tr2Short, xcorrCfg.maxDelay, xcorr_lag, xcorr_coeff);
  }

  // Trust the manual pick on `phase1`: keep `tr1` short and cross-correlate
  // it with a larger `tr2` window.
  double xcorr_coeff2 = numeric_limits<double>::quiet_NaN(), xcorr_lag2 = 0;

  if (phase1.trusted() || (!phase1.trusted() && !phase2.trusted()))
  {
    // Trim `tr1` to shorter length; we want to cross-correlate the short with
    // the long one.
    Trace tr1Short(*tr1);
    TimeWindow tw1Short = xcorrTimeWindowShort(xcorrOpt, event1, phase1);
    if (!tr1Short.slice(tw1Short))
    {
      logDebugF("Skipping cross-correlation: cannot trim phase1 waveform "
                "for phase pair phase1='%s', phase2='%s'",
                string(phase1).c_str(), string(phase2).c_str());
      return false;
    }

    xcorr(tr1Short, *tr2, xcorrCfg.maxDelay, xcorr_lag2, xcorr_coeff2);
  }

  if (!std::isfinite(xcorr_coeff) && !std::isfinite(xcorr_coeff2))
  {
    logDebugF("Skipping cross-correlation: no meaningful coefficient for phase "
              "pair phase1='%s', phase2='%s'",
              string(phase1).c_str(), string(phase2).c_str());
    return false;
  }
  else if (std::isfinite(xcorr_coeff2))
  {
    if (!std::isfinite(xcorr_coeff) ||
        (std::abs(xcorr_coeff2) > std::abs(xcorr_coeff)))
    {
      xcorr_coeff = xcorr_coeff2;
      xcorr_lag   = xcorr_lag2;
    }
  }

  coeffOut = std::abs(xcorr_coeff);
  lagOut   = xcorr_lag;

  return true;
}

/*
 * Compute cross-correlation between two traces centered around their
 * respective picks. The cross-correlation will be performed between the
 * short trace against the longest trace, in the interval longest trace
 * middle point -/+ 'maxDelay' seconds.
 * `lagOut` will store the correction (seconds) to be applied to the pick time
 * difference (tr1 - tr2) to obtain the highest correlation coefficient (in
 * absolute value), stored in 'coeffOut'
 */
void DD::xcorr(const Trace &tr1,
               const Trace &tr2,
               double maxDelay,
               double &lagOut,
               double &coeffOut)
{
  const double freq = tr1.samplingFrequency();

  // check longest/shortest trace
  const bool swap        = tr1.sampleCount() > tr2.sampleCount();
  const Trace &trShorter = swap ? tr2 : tr1;
  const Trace &trLonger  = swap ? tr1 : tr2;

  const double *dataS = trShorter.data();
  const double *dataL = trLonger.data();
  const int sizeS     = trShorter.sampleCount();
  const int sizeL     = trLonger.sampleCount();

  // force to cross-correlate withing data boundaries
  int availableData = (sizeL - sizeS) / 2;
  int maxDelaySmps  = maxDelay * freq;
  if (maxDelaySmps > availableData) maxDelaySmps = availableData;

  crossCorrelation(dataS, sizeS, (dataL + availableData - maxDelaySmps),
                   (sizeS + maxDelaySmps * 2), lagOut, coeffOut);

  lagOut -= maxDelaySmps; // the reference is the middle of the long trace
  lagOut /= freq;         // samples to secs

  if (swap) lagOut = -lagOut;
}

shared_ptr<const Trace> DD::getWaveform(Waveform::Processor &wfProc,
                                        const TimeWindow &tw,
                                        const Event &ev,
                                        const Phase &ph,
                                        const string &component)
{
  const Station &sta = _bgCat.getStations().at(ph.stationId);

  Catalog::Phase phCopy(ph);
  Transform trans;
  if (component.empty())
  {
    return nullptr;
  }
  else if (component == "T")
  {
    trans = Transform::TRANSVERSAL;
  }
  else if (component == "R")
  {
    trans = Transform::RADIAL;
  }
  else if (component == "L2")
  {
    trans = Transform::L2;
  }
  else
  {
    trans              = Transform::NONE;
    phCopy.channelCode = getBandAndInstrumentCodes(ph.channelCode) + component;
  }

  return wfProc.get(tw, phCopy, ev, sta, _cfg.wfFilter.filterStr,
                    _cfg.wfFilter.resampleFreq, trans);
}

void DD::logXCorrSummary(const unordered_map<unsigned, Neighbours> &cluster,
                         const XcorrOptions &xcorrOpt,
                         const XCorrCache &xcorr)
{
  struct Counters
  {
    unsigned performed      = 0;
    unsigned failed         = 0;
    unsigned aboveThreshold = 0;
    vector<double> coeff;
    vector<double> lag;
  };

  map<string, Counters> pCountByStation; // key station id
  map<string, Counters> sCountByStation; // key station id

  for (const auto &kv : cluster)
  {
    const Neighbours &neighbours = kv.second;
    for (const auto &t : neighbours.phases())
    {
      const string &stationId  = std::get<0>(t);
      const string &phase      = std::get<1>(t);
      const unsigned neighEvId = std::get<2>(t);

      Catalog::Phase::Type phaseType = Phase::Type::NO;
      // P phase
      auto itpp =
          find(_cfg.validPphases.begin(), _cfg.validPphases.end(), phase);
      if (itpp != _cfg.validPphases.end())
      {
        phaseType = Phase::Type::P;
      }
      else
      {
        // S phase
        auto itsp =
            find(_cfg.validSphases.begin(), _cfg.validSphases.end(), phase);
        if (itsp != _cfg.validSphases.end())
        {
          phaseType = Phase::Type::S;
        }
      }

      if (phaseType == Phase::Type::NO)
      {
        logWarning("Unknown phase in xcorr stats");
        continue;
      }

      if (!xcorr.has(neighbours.referenceId(), neighEvId, stationId, phase))
      {
        continue;
      }
      const auto &e =
          xcorr.get(neighbours.referenceId(), neighEvId, stationId, phase);
      Counters &counters = (phaseType == Phase::Type::P)
                               ? pCountByStation[stationId]
                               : sCountByStation[stationId];

      if (e.valid)
      {
        counters.performed++;
        counters.coeff.push_back(e.coeff);
        counters.lag.push_back(e.lag);

        const auto &xcorrCfg = xcorrOpt.phase.at(phaseType);
        if (e.coeff >= xcorrCfg.minCoef)
        {
          counters.aboveThreshold++;
        }
      }
      else
      {
        counters.failed++;
      }
    }
  }

  //
  // P phase
  //
  vector<string> logPLines;
  unsigned performed      = 0;
  unsigned failed         = 0;
  unsigned aboveThreshold = 0;
  for (const auto &kv : pCountByStation)
  {
    const string &stationId = kv.first;
    const Counters &c       = kv.second;

    double min, max, q1, q2, q3;
    compute5numberSummary(c.coeff, min, max, q1, q2, q3);

    unsigned tot = c.performed + c.failed;

    logPLines.push_back(
        strf("%-11s %9u  %6.1f%%     %6.1f%%              %4.2f  %4.2f  %4.2f  "
             "%4.2f  %4.2f",
             stationId.c_str(), tot, (c.failed * 100.0 / tot),
             (c.aboveThreshold * 100.0 / tot), min, q1, q2, q3, max));

    performed += c.performed;
    failed += c.failed;
    aboveThreshold += c.aboveThreshold;
  }
  unsigned tot = performed + failed;
  logInfoF("Total P phases CC attempted %u failures %.1f%% (%u) "
           "Above coeff threshold (%.2f) %.1f%% (%u)",
           tot, (failed * 100.0 / tot), failed,
           xcorrOpt.phase.at(Phase::Type::P).minCoef,
           (aboveThreshold * 100.0 / tot), aboveThreshold);
  //
  // S phase
  //
  vector<string> logSLines;
  performed      = 0;
  failed         = 0;
  aboveThreshold = 0;
  for (const auto &kv : sCountByStation)
  {
    const string &stationId = kv.first;
    const Counters &c       = kv.second;

    double min, max, q1, q2, q3;
    compute5numberSummary(c.coeff, min, max, q1, q2, q3);

    unsigned tot = c.performed + c.failed;

    logSLines.push_back(
        strf("%-11s %9u  %6.1f%%     %6.1f%%              %4.2f  %4.2f  %4.2f  "
             "%4.2f  %4.2f",
             stationId.c_str(), tot, (c.failed * 100.0 / tot),
             (c.aboveThreshold * 100.0 / tot), min, q1, q2, q3, max));

    performed += c.performed;
    failed += c.failed;
    aboveThreshold += c.aboveThreshold;
  }
  tot = performed + failed;
  logInfoF("Total S phases CC attempted %u failures %.1f%% (%u) "
           "Above coeff threshold (%.2f) %.1f%% (%u)",
           tot, (failed * 100.0 / tot), failed,
           xcorrOpt.phase.at(Phase::Type::S).minCoef,
           (aboveThreshold * 100.0 / tot), aboveThreshold);

  logInfoF("(P phases)   Total CC   Failures   Coeff>%.2f            Min  1stQ "
           "Median 3rdQ  Max",
           xcorrOpt.phase.at(Phase::Type::P).minCoef);
  for (const string &line : logPLines)
  {
    logInfo(line);
  }

  logInfoF("(S phases)   Total CC   Failures   Coeff>%.2f            Min  1stQ "
           "Median 3rdQ  Max",
           xcorrOpt.phase.at(Phase::Type::S).minCoef);
  for (const string &line : logSLines)
  {
    logInfo(line);
  }
}

} // namespace HDD

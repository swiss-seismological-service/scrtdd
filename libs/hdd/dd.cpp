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
    os << HDD::strf("%u,%u,%s,%s,%s,%c,%s,%g,%g,%g", e.evId1, e.evId2,
                    sta.networkCode.c_str(), sta.stationCode.c_str(),
                    sta.locationCode.c_str(), e.phase,
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
    : _cfg(cfg), _srcCat(catalog),
      _bgCat(Catalog::filterPhasesAndSetWeights(_srcCat,
                                                Phase::Source::CATALOG,
                                                _cfg.validPphases,
                                                _cfg.validSphases,
                                                _cfg.pickUncertaintyClasses)),
      _ttt(std::move(ttt)), _wf(std::move(wf))
{
  disableCatalogWaveformDiskCache();
  disableAllWaveformDiskCache();
}

void DD::disableCatalogWaveformDiskCache()
{
  _useCatalogWaveformDiskCache = false;
  createWaveformCache();
}

void DD::enableCatalogWaveformDiskCache(const std::string &cacheDir)
{
  _useCatalogWaveformDiskCache = true;
  _cacheDir                    = cacheDir;
  if (!pathExists(_cacheDir))
  {
    if (!createDirectories(_cacheDir))
    {
      string msg = "Unable to create cache directory: " + _cacheDir;
      throw Exception(msg);
    }
  }
  createWaveformCache();
}

void DD::disableAllWaveformDiskCache() { _waveformCacheAll = false; }

void DD::enableAllWaveformDiskCache(const std::string &tmpCacheDir)
{
  _waveformCacheAll = true;
  _tmpCacheDir      = tmpCacheDir;

  if (!pathExists(_tmpCacheDir))
  {
    if (!createDirectories(_tmpCacheDir))
    {
      string msg = "Unable to create cache directory: " + _tmpCacheDir;
      throw Exception(msg);
    }
  }
}

void DD::createWaveformCache()
{
  _wfAccess.loader    = make_shared<Waveform::BasicLoader>(_wf);
  _wfAccess.diskCache = nullptr;
  _wfAccess.memCache  = nullptr;

  shared_ptr<Waveform::Loader> currLdr = _wfAccess.loader;

  if (_useCatalogWaveformDiskCache)
  {
    _wfAccess.diskCache =
        make_shared<Waveform::DiskCachedLoader>(_wf, currLdr, _cacheDir);
    currLdr = make_shared<Waveform::ExtraLenLoader>(_wfAccess.diskCache,
                                                    _cfg.diskTraceMinLen);
  }

  shared_ptr<Waveform::Processor> currProc =
      make_shared<Waveform::BasicProcessor>(_wf, currLdr,
                                            _cfg.wfFilter.extraTraceLen);

  // Using the MemCachedProc mechanism of wrapping Waveform::Processor it is
  // possible to add additional Waveform::Processor(S) with specialized
  // operations between currProc and memCache. For example there use to be a SNR
  // filter there

  _wfAccess.memCache = make_shared<Waveform::MemCachedProc>(currProc);
}

void DD::replaceWaveformLoader(const shared_ptr<Waveform::Loader> &baseLdr)
{
  if (_useCatalogWaveformDiskCache)
  {
    _wfAccess.diskCache->setAuxLoader(baseLdr);
  }
  else
  {
    _wfAccess.memCache->setAuxProcessor(make_shared<Waveform::BasicProcessor>(
        _wf, baseLdr, _cfg.wfFilter.extraTraceLen));
  }
}

string DD::generateWorkingSubDir(const string &prefix) const
{
  UTCTime now = UTCClock::now();
  int year, month, day, hour, min, sec, usec;
  UTCClock::toDate(now, year, month, day, hour, min, sec, usec);

  UniformRandomer ran(0, 9999);
  ran.setSeed(durToSec(now.time_since_epoch()));

  // preifx_creationTime_randomNumber
  string id = strf("%s_%04d%02d%02d%02d%02d%02d_%04zu", prefix.c_str(), year,
                   month, day, hour, min, sec, ran.next());
  return id;
}

string DD::generateWorkingSubDir(const Event &ev) const
{
  int year, month, day, hour, min, sec, usec;
  UTCClock::toDate(ev.time, year, month, day, hour, min, sec, usec);
  // singleevent_eventTime_latitude_longitude
  string prefix =
      strf("singleevent_%04d%02d%02d%02d%02d%02d_%05d_%06d", year, month, day,
           hour, min, sec, int(ev.latitude * 1000), int(ev.longitude * 1000));
  return generateWorkingSubDir(prefix);
}

void DD::unloadWaveforms() { createWaveformCache(); }

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

void DD::preloadWaveforms(const XcorrOptions &xcorrOpt)
{
  //
  // preload waveforms, store them on disk and cache them in memory
  // (already processed).
  // For better performance we want to load waveforms in batch
  // (batchloader)
  //
  logInfoF("Loading catalog waveform data (%zu events to load)",
           _bgCat.getEvents().size());

  auto forEventWaveforms =
      [this,
       &xcorrOpt](const Event &event, unsigned &numPhases, unsigned &numSPhases,
                  std::function<bool(const TimeWindow &, const Event &,
                                     const Phase &, const string &)> func) {
        auto eqlrng = _bgCat.getPhases().equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
          const Phase &phase = it->second;
          TimeWindow tw      = xcorrTimeWindowLong(xcorrOpt, phase);

          for (const string &component : xcorrComponents(xcorrOpt, phase))
          {
            if (func(tw, event, phase, component)) break;
          }

          numPhases++;
          if (phase.procInfo.type == Phase::Type::S) numSPhases++;
        }
      };

  unsigned numPhases = 0, numSPhases = 0, numEvents = 0;
  WfCounters wfcount{};

  for (const auto &kv : _bgCat.getEvents())
  {
    const Event &event = kv.second;

    logDebugF("Loading event %u waveforms...", event.id);

    // Use a BatchLoader for the current event
    shared_ptr<Waveform::BatchLoader> batchLoader =
        make_shared<Waveform::BatchLoader>(_wf);
    replaceWaveformLoader(batchLoader);

    // The following call won't try to load the waveforms. Instead it will
    // records in BatchLoader the required waveforms that will be loaded
    // in batch later
    auto requestWaveform = [this](const TimeWindow &tw, const Event &ev,
                                  const Phase &ph, const string &component) {
      getWaveform(*_wfAccess.memCache, tw, ev, ph, component);
      return false;
    };

    forEventWaveforms(event, numPhases, numSPhases, requestWaveform);

    // This will actually download the waveworms
    batchLoader->load();

    // now process and store the waveforms into memory
    auto cacheWaveform = [this](const TimeWindow &tw, const Event &ev,
                                const Phase &ph, const string &component) {
      return getWaveform(*_wfAccess.memCache, tw, ev, ph, component) != nullptr;
    };

    unsigned dummy;
    forEventWaveforms(event, dummy, dummy, cacheWaveform);

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
  replaceWaveformLoader(make_shared<Waveform::BasicLoader>(_wf));

  wfcount.update(_wfAccess.loader.get());
  wfcount.update(_wfAccess.diskCache.get());

  logInfoF("Finished preloading catalog waveform data: total events %zu total "
           "phases %u (P %.f%%, S %.f%%). Waveforms downloaded %u, not "
           "available %u, loaded from disk cache %u",
           _bgCat.getEvents().size(), numPhases,
           ((numPhases - numSPhases) * 100. / numPhases),
           (numSPhases * 100. / numPhases), wfcount.downloaded,
           wfcount.no_avail, wfcount.disk_cached);
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
      TimeWindow tw      = xcorrTimeWindowLong(xcorrOpt, phase);

      for (const string &component : xcorrComponents(xcorrOpt, phase))
      {
        // memCache->getAuxProcessor() -> we don't want to cache the
        // whole catalog waveforms into memory
        shared_ptr<const Trace> tr =
            getWaveform(*_wfAccess.memCache->getAuxProcessor(), tw, event,
                        phase, component);

        if (!tr) continue;

        string channelCode =
            getBandAndInstrumentCodes(phase.channelCode) + component;

        const string wfPath = waveformPath(basePath, event, phase, channelCode);

        logInfoF("Writing %s", wfPath.c_str());
        try
        {

          _wf->writeTrace(*tr, wfPath);
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
          _bgCat, clustOpt.minESdist, clustOpt.maxESdist,
          clustOpt.minEStoIEratio, clustOpt.minDTperEvt, clustOpt.maxDTperEvt,
          clustOpt.minNumNeigh, clustOpt.maxNumNeigh, clustOpt.numEllipsoids,
          clustOpt.maxEllipsoidSize, true);

  // Organize the neighbours by not connected clusters
  return clusterizeNeighbouringEvents(neighboursByEvent);
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
  logInfo("Starting DD relocator in multiple events mode");

  Catalog catToReloc(_bgCat);

  logInfoF("The catalog contains %zu events", catToReloc.getEvents().size());

  if (saveProcessing)
  {
    // prepare a folder for processing files
    if (processingDataDir.empty())
    {
      do
      {
        processingDataDir = generateWorkingSubDir("multievent");
      } while (pathExists(processingDataDir));
    }

    if (!pathExists(processingDataDir))
    {
      if (!createDirectories(processingDataDir))
      {
        string msg = "Unable to create directory: " + processingDataDir;
        throw Exception(msg);
      }
    }
    logInfoF("Processing data dir %s", processingDataDir.c_str());
  }

  // prepare file logger
  string logFile = joinPath(processingDataDir, "relocation.log");
  if (saveProcessing)
  {
    addFileLogger(logFile, Level::info);
  }

  if (clusters.empty())
  {
    // find Neighbours for each event in the catalog
    unordered_map<unsigned, Neighbours> neighboursByEvent =
        selectNeighbouringEventsCatalog(
            catToReloc, clustOpt.minESdist, clustOpt.maxESdist,
            clustOpt.minEStoIEratio, clustOpt.minDTperEvt, clustOpt.maxDTperEvt,
            clustOpt.minNumNeigh, clustOpt.maxNumNeigh, clustOpt.numEllipsoids,
            clustOpt.maxEllipsoidSize, true);

    // Organize the neighbours by non-connected clusters
    clusters = clusterizeNeighbouringEvents(neighboursByEvent);
  }

  logInfoF("Found %zu event clusters with the following number of events:",
           clusters.size());
  for (const auto &cluster : clusters)
  {
    logInfoF(" %zu events", cluster.size());
  }

  if (saveProcessing)
  {
    catToReloc.writeToFile(joinPath(processingDataDir, "input-event.csv"),
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

    if (saveProcessing)
    {
      string prefix = strf("cluster-%u", clusterId);
      Neighbours::writeToFile(cluster, catToReloc, prefix + "-pair.csv");
      Catalog catToDump;
      for (const auto &kv : cluster)
      {
        catToDump.add(kv.first, catToReloc, true);
      }
      catToDump.writeToFile(
          joinPath(processingDataDir, (prefix + "-event.csv")),
          joinPath(processingDataDir, (prefix + "-phase.csv")));
    }

    // Perform cross-correlation
    XCorrCache xcorr =
        buildXCorrCache(catToReloc, cluster, xcorrOpt, xcorrData);

    // the actual relocation
    std::vector<Solver::DoubleDifference> startDDs, finalDDs;
    Catalog relocatedCluster =
        relocate(catToReloc, cluster, xcorrOpt, solverOpt, false, xcorr,
                 startDDs, finalDDs);

    relocatedCatalog.add(relocatedCluster, true);

    if (saveProcessing)
    {
      string prefix = strf("cluster-%u", clusterId);
      writeDoubleDifferenceToFile(startDDs, catToReloc,
                                  prefix + "initial-double-difference");
      writeDoubleDifferenceToFile(finalDDs, catToReloc,
                                  prefix + "final-double-difference");
      relocatedCluster.writeToFile(
          joinPath(processingDataDir, (prefix + "-relocated-event.csv")),
          joinPath(processingDataDir, (prefix + "-relocated-phase.csv")));
    }
    clusterId++;

    // make sure xcorrData stores all the cross-correlation results
    xcorrData.add(xcorr);
  }

  if (saveProcessing)
  {
    relocatedCatalog.writeToFile(
        joinPath(processingDataDir, "relocated-event.csv"),
        joinPath(processingDataDir, "relocated-phase.csv"));
    xcorrData.writeToFile(_bgCat, joinPath(processingDataDir, "xcorr.csv"));
  }

  removeFileLogger(logFile);

  return relocatedCatalog;
}

Catalog DD::relocateSingleEvent(const Catalog &singleEvent,
                                bool isManual,
                                const ClusteringOptions &clustOpt1,
                                const ClusteringOptions &clustOpt2,
                                const XcorrOptions &xcorrOpt,
                                const SolverOptions &solverOpt,
                                bool saveProcessing,
                                string processingDataDir)
{
  const Catalog &bgCat = _bgCat;

  // there must be only one event in the catalog, the origin to relocate
  const Event &evToRelocate = singleEvent.getEvents().begin()->second;
  const auto &evToRelocatePhases =
      singleEvent.getPhases().equal_range(evToRelocate.id);

  logInfoF("Starting DD relocator in single event mode: event %s lat %.6f lon "
           "%.6f depth %.4f mag %.2f time %s #phases %ld",
           string(evToRelocate).c_str(), evToRelocate.latitude,
           evToRelocate.longitude, evToRelocate.depth, evToRelocate.magnitude,
           UTCClock::toString(evToRelocate.time).c_str(),
           std::distance(evToRelocatePhases.first, evToRelocatePhases.second));

  if (saveProcessing)
  {
    // prepare a folder for processing files
    if (processingDataDir.empty())
    {
      do
      {
        processingDataDir = generateWorkingSubDir("multievent");
      } while (pathExists(processingDataDir));
    }

    if (!pathExists(processingDataDir))
    {
      if (!createDirectories(processingDataDir))
      {
        string msg = "Unable to create directory: " + processingDataDir;
        throw Exception(msg);
      }
    }
    logInfoF("Processing data dir %s", processingDataDir.c_str());
  }

  // prepare file logger
  string logFile = joinPath(processingDataDir, "relocation.log");
  if (saveProcessing)
  {
    addFileLogger(logFile, Level::info);
  }

  logInfo("Performing step 1: initial location refinement (no "
          "cross-correlation)");

  string eventWorkingDir = joinPath(processingDataDir, "step1");

  Catalog evToRelocateCat = Catalog::filterPhasesAndSetWeights(
      singleEvent,
      (isManual ? Phase::Source::RT_EVENT_MANUAL
                : Phase::Source::RT_EVENT_AUTOMATIC),
      _cfg.validPphases, _cfg.validSphases, _cfg.pickUncertaintyClasses);

  XcorrOptions xcorrOptDisabled(xcorrOpt);
  xcorrOptDisabled.enable = false;

  Catalog relocatedEvCat = relocateEventSingleStep(
      bgCat, evToRelocateCat, clustOpt1, xcorrOptDisabled, solverOpt,
      saveProcessing, eventWorkingDir);

  if (!relocatedEvCat.empty())
  {
    const Event &ev = relocatedEvCat.getEvents().begin()->second;
    logInfoF("Step 1 relocation successful, new location: "
             "lat %.6f lon %.6f depth %.4f time %s",
             ev.latitude, ev.longitude, ev.depth,
             UTCClock::toString(ev.time).c_str());

    evToRelocateCat = std::move(relocatedEvCat);
  }
  else
  {
    logError("Failed to perform step 1 origin relocation");
  }

  logInfo("Performing step 2: relocate again, possibly with cross-correlation");

  eventWorkingDir = joinPath(processingDataDir, "step2");

  Catalog relocatedEvWithXcorr =
      relocateEventSingleStep(bgCat, evToRelocateCat, clustOpt2, xcorrOpt,
                              solverOpt, saveProcessing, eventWorkingDir);

  if (!relocatedEvWithXcorr.empty())
  {
    Event ev = relocatedEvWithXcorr.getEvents().begin()->second;
    logInfoF("Step 2 relocation successful, new location: "
             "lat %.6f lon %.6f depth %.4f time %s",
             ev.latitude, ev.longitude, ev.depth,
             UTCClock::toString(ev.time).c_str());

    // update the "origin change information" taking into consideration
    // the first relocation step, too
    if (!relocatedEvCat.empty())
    {
      const Event &prevRelocEv = relocatedEvCat.getEvents().begin()->second;
      if (prevRelocEv.relocInfo.isRelocated)
      {
        ev.relocInfo.startRms = prevRelocEv.relocInfo.startRms;
        relocatedEvWithXcorr.updateEvent(ev);
      }
    }
  }
  else
  {
    logError("Failed to perform step 2 origin relocation");
  }

  if (relocatedEvWithXcorr.empty())
  {
    throw Exception("Failed origin relocation");
  }
  removeFileLogger(logFile);

  return relocatedEvWithXcorr;
}

Catalog DD::relocateEventSingleStep(const Catalog &bgCat,
                                    const Catalog &evToRelocateCat,
                                    const ClusteringOptions &clustOpt,
                                    const XcorrOptions &xcorrOpt,
                                    const SolverOptions &solverOpt,
                                    bool saveProcessing,
                                    string processingDataDir)
{
  if (saveProcessing)
  {
    if (!createDirectories(processingDataDir))
    {
      string msg = "Unable to create working directory: " + processingDataDir;
      throw Exception(msg);
    }
    logInfoF("Working dir %s", processingDataDir.c_str());

    evToRelocateCat.writeToFile(
        joinPath(processingDataDir, "single-event.csv"),
        joinPath(processingDataDir, "single-event-phase.csv"),
        joinPath(processingDataDir, "single-event-station.csv"));
  }

  Catalog relocatedEvCat;

  try
  {
    // extract event to relocate
    const Event &evToRelocate = evToRelocateCat.getEvents().begin()->second;

    //
    // select neighbouring events
    //
    Neighbours neighbours = selectNeighbouringEvents(
        bgCat, evToRelocate, evToRelocateCat, clustOpt.minESdist,
        clustOpt.maxESdist, clustOpt.minEStoIEratio, clustOpt.minDTperEvt,
        clustOpt.maxDTperEvt, clustOpt.minNumNeigh, clustOpt.maxNumNeigh,
        clustOpt.numEllipsoids, clustOpt.maxEllipsoidSize, false);

    logInfoF("Found %zu neighbouring events", neighbours.ids().size());

    //
    // prepare catalog to relocate
    //
    Catalog catalog = neighbours.toCatalog(bgCat);
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
      Neighbours::writeToFile(cluster, catalog, "pair.csv");
    }

    // Perform cross-correlation
    XCorrCache xcorr = buildXCorrCache(catalog, cluster, xcorrOpt);

    // the actual relocation
    std::vector<Solver::DoubleDifference> startDDs, finalDDs;
    relocatedEvCat = relocate(catalog, cluster, xcorrOpt, solverOpt, true,
                              xcorr, startDDs, finalDDs);

    if (saveProcessing)
    {
      writeDoubleDifferenceToFile(startDDs, catalog,
                                  "initial-double-difference");
      writeDoubleDifferenceToFile(finalDDs, catalog, "final-double-difference");
      relocatedEvCat.writeToFile(
          joinPath(processingDataDir, "relocated-event.csv"),
          joinPath(processingDataDir, "relocated-phase.csv"));
    }
  }
  catch (Exception &e)
  {
    logError(e.what());
  }

  return relocatedEvCat;
}

Catalog DD::relocate(const Catalog &catalog,
                     const unordered_map<unsigned, Neighbours> &cluster,
                     const XcorrOptions &xcorrOpt,
                     const SolverOptions &solverOpt,
                     bool keepNeighboursFixed,
                     const XCorrCache &xcorr,
                     std::vector<Solver::DoubleDifference> startDDs,
                     std::vector<Solver::DoubleDifference> finalDDs) const
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
        logDebugF("Event %u relocated, but it is now outside the travel time "
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
    computeEventResiduals(solver, currentCatalog, solverOpt, cluster,
                          isFirstIteration);

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
  // copy event because we'll update it
  const Event &refEv = catalog.getEvents().at(neighbours.referenceId());
  const unordered_set<unsigned> neighboursIds = neighbours.ids();

  // Detect an event that moved to a location outside the ttt range
  bool tttAttempted = false, tttAvailable = false;

  //
  // Loop through reference event phases
  // Instead of looping though the Neighbours phases as expected, we loop though
  // each reference event phase and check if there is an entry in Neighbours
  // This approach allows to search only once in the catalog and then reuse the
  // refEv phases for all neighbours (small optimization in case of catalogs
  // with millions of phases)
  //
  auto eqlrng = catalog.getPhases().equal_range(refEv.id);
  for (auto it = eqlrng.first; it != eqlrng.second; ++it)
  {
    const Phase &refPhase  = it->second;
    const Station &station = catalog.getStations().at(refPhase.stationId);
    char phaseTypeAsChar   = static_cast<char>(refPhase.procInfo.type);
    const auto &xcorrCfg   = xcorrOpt.phase.at(refPhase.procInfo.type);

    //
    // loop through neighbouring events and look for the matching phase
    //
    for (unsigned neighEvId : neighboursIds)
    {
      if (!neighbours.has(neighEvId, refPhase.stationId,
                          refPhase.procInfo.type))
        continue;

      const Event &event = catalog.getEvents().at(neighEvId);

      const Phase &phase =
          catalog
              .searchPhase(event.id, refPhase.stationId, refPhase.procInfo.type)
              ->second;
      //
      // compute travel times for both event and `refEv`
      //
      UTCTime::duration ref_travel_time = refPhase.time - refEv.time;
      if (ref_travel_time.count() < 0)
      {
        logDebugF("Ignoring phase %s with negative travel time",
                  string(refPhase).c_str());
        continue;
      }

      UTCTime::duration travel_time = phase.time - event.time;
      if (travel_time.count() < 0)
      {
        logDebugF("Ignoring phase %s with negative travel time",
                  string(phase).c_str());
        continue;
      }

      tttAttempted = true; // at least one attempt at loading ttt
                           //
      if (!addObservationParams(solver, *_ttt, refEv, station, refPhase,
                                true) ||
          !addObservationParams(solver, *_ttt, event, station, phase,
                                !keepNeighboursFixed))
      {
        logDebugF("Skipping observation (ev %u-%u sta %s phase %c)", refEv.id,
                  event.id, station.id.c_str(), phaseTypeAsChar);
        continue;
      }

      tttAvailable = true; // at least once successfully ttt query

      //
      // compute absolute travel time differences to the solver
      //
      duration<double> diffTime = ref_travel_time - travel_time;
      double weight =
          usePickUncertainties
              ? ((refPhase.procInfo.classWeight + phase.procInfo.classWeight) /
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

        if (xcorr.has(refEv.id, event.id, refPhase.stationId,
                      refPhase.procInfo.type))
        {
          const auto &e = xcorr.get(refEv.id, event.id, refPhase.stationId,
                                    refPhase.procInfo.type);
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

      solver.addObservation(refEv.id, event.id, refPhase.stationId,
                            phaseTypeAsChar, diffTime.count(), weight,
                            xcorrUsed, xcorrCoeff);
    }
  }

  return !tttAttempted || tttAvailable;
}

bool DD::addObservationParams(Solver &solver,
                              TravelTimeTable &ttt,
                              const Catalog::Event &event,
                              const Catalog::Station &station,
                              const Catalog::Phase &phase,
                              bool computeEvChanges) const
{
  char phaseType = static_cast<char>(phase.procInfo.type);

  bool dummy1;
  double dummy2;
  //
  // check if this data has been already added to the solver
  //
  if (!solver.getObservationParams(event.id, station.id, phaseType, dummy1,
                                   dummy2, dummy2, dummy2, dummy2, dummy2))
  {
    //
    // Compute travel time information
    //
    double travelTime, takeOfAngleAzim, takeOfAngleDip, velocityAtSrc;
    try
    {
      ttt.compute(event, station, string(1, phaseType), travelTime,
                  takeOfAngleAzim, takeOfAngleDip, velocityAtSrc);
    }
    catch (Exception &e)
    {
      logWarningF(
          "Travel Time Table error: %s (Event lat %.6f lon %.6f depth %.6f "
          "Station lat %.6f lon %.6f elevation %.f )",
          e.what(), event.latitude, event.longitude, event.depth,
          station.latitude, station.longitude, station.elevation);
      return false;
    }
    double travelTimeResidual = travelTime - durToSec(phase.time - event.time);

    //
    // Populate the solver
    //
    solver.addEvent(event.id, event.latitude, event.longitude, event.depth);
    solver.addStation(station.id, station.latitude, station.longitude,
                      station.elevation);
    solver.addObservationParams(
        event.id, station.id, phaseType, computeEvChanges, travelTime,
        travelTimeResidual, takeOfAngleAzim, takeOfAngleDip, velocityAtSrc);
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
      char phaseTypeAsChar   = static_cast<char>(phase.procInfo.type);

      if (isFirstIteration)
      {
        phase.relocInfo.startResidual = 0;
      }

      phase.relocInfo.finalResidual = 0;
      phase.relocInfo.weight =
          solverOpt.usePickUncertainties ? phase.procInfo.classWeight : 1;

      bool dummy1;
      double residual, dummy2;

      if (solver.getObservationParams(event.id, station.id, phaseTypeAsChar,
                                      dummy1, dummy2, residual, dummy2, dummy2,
                                      dummy2))
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
    double deltaLat, deltaLon, deltaDepth, deltaTime;
    if (!solver.getEventChanges(event.id, deltaLat, deltaLon, deltaDepth,
                                deltaTime))
    {
      continue;
    }

    double newDepth = event.depth + deltaDepth;
    UTCTime newTime = event.time + secToDur(deltaTime);

    //
    // converte delta lat lon [ km -> degree ]
    //
    double distance = std::sqrt(square(deltaLat) + square(deltaLon));
    double azimuth  = std::atan2(deltaLon, deltaLat);

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
      char phaseTypeAsChar   = static_cast<char>(phase.procInfo.type);

      if (isFirstIteration)
      {
        phase.relocInfo.isRelocated = false;
      }

      if (!solver.isEventPhaseUsed(event.id, station.id, phaseTypeAsChar))
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
DD::buildXCorrCache(Catalog &catalog,
                    const unordered_map<unsigned, Neighbours> &cluster,
                    const XcorrOptions &xcorrOpt,
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

    buildXcorrDiffTTimePairs(catalog, neighbours, refEv, xcorrOpt, precomputed,
                             xcorr);

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
void DD::buildXcorrDiffTTimePairs(Catalog &catalog,
                                  const Neighbours &neighbours,
                                  const Event &refEv,
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
  // make use of the background catalog cache (_wfAccess.memCache) to
  // load them
  //
  shared_ptr<Waveform::Processor> tmpMemCache(
      preloadNonCatalogWaveforms(catalog, neighbours, refEv, xcorrOpt));

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
    shared_ptr<Waveform::Processor> refPrc;
    if (refPhase.procInfo.source == Phase::Source::CATALOG)
    {
      refPrc = _wfAccess.memCache;
    }
    else if (refPhase.procInfo.source == Phase::Source::RT_EVENT_MANUAL)
    {
      refPrc = tmpMemCache;
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
      {
        continue;
      }

      if (neighbours.has(neighEvId, refPhase.stationId, refPhase.procInfo.type))
      {
        const Phase &phase = catalog
                                 .searchPhase(event.id, refPhase.stationId,
                                              refPhase.procInfo.type)
                                 ->second;

        // In single-event mode `refPhase` is real-time and `phase` is from
        // the catalog. In multi-event mode both are from the catalog.
        if (phase.procInfo.source != Phase::Source::CATALOG)
        {
          throw Exception("Internal logic error: phase is not from catalog");
        }

        // skip cross-correlation if we already have the result in cache
        if (!xcorr.has(refEv.id, event.id, refPhase.stationId,
                       refPhase.procInfo.type))
        {
          try
          {
            // fetch the cross-correlation results from the precomputed values
            const XCorrCache::Entry &e = precomputed.get(
                refEv.id, event.id, refPhase.stationId, refPhase.procInfo.type);
            xcorr.add(refEv.id, event.id, refPhase.stationId,
                      refPhase.procInfo.type, e.valid, e.coeff, e.lag);
          }
          catch (...)
          {
            // Do the actual cross-correlation
            double coeff, lag;
            string component;
            bool valid =
                xcorrPhases(xcorrOpt, refEv, refPhase, *refPrc, event, phase,
                            *_wfAccess.memCache, coeff, lag, component);

            // cache cross-correlation results
            xcorr.add(refEv.id, event.id, refPhase.stationId,
                      refPhase.procInfo.type, valid, coeff, lag);
          }
        }
      }
    }
  }
}

shared_ptr<Waveform::Processor>
DD::preloadNonCatalogWaveforms(Catalog &catalog,
                               const Neighbours &neighbours,
                               const Event &refEv,
                               const XcorrOptions &xcorrOpt)
{
  //
  // For single-event relocation in real-time we want to load the waveforms
  // in batch otherwise the seedlink server gets stuck and becomes
  // unresponsive due to the multiple connections requests, one for each event
  // phase
  //
  shared_ptr<Waveform::BatchLoader> batchLoader(new Waveform::BatchLoader(_wf));

  shared_ptr<Waveform::Loader> loader = batchLoader;
  shared_ptr<Waveform::DiskCachedLoader> diskLoader;
  if (_useCatalogWaveformDiskCache && _waveformCacheAll)
  {
    diskLoader.reset(
        new Waveform::DiskCachedLoader(_wf, batchLoader, _tmpCacheDir));
    loader.reset(
        new Waveform::ExtraLenLoader(diskLoader, _cfg.diskTraceMinLen));
  }

  shared_ptr<Waveform::Processor> proc(
      new Waveform::BasicProcessor(_wf, loader, _cfg.wfFilter.extraTraceLen));

  //
  // loop through reference event phases
  //
  bool onlyCatalogPhases = true;
  auto eqlrngRef         = catalog.getPhases().equal_range(refEv.id);
  for (auto itRef = eqlrngRef.first; itRef != eqlrngRef.second; ++itRef)
  {
    const Phase &refPhase  = itRef->second;
    const Station &station = catalog.getStations().at(refPhase.stationId);
    const auto components  = xcorrComponents(xcorrOpt, refPhase);

    // We deal only with real-time event data
    if (refPhase.procInfo.source == Phase::Source::CATALOG)
    {
      continue;
    }
    onlyCatalogPhases = false;

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

      if (neighbours.has(neighEvId, refPhase.stationId, refPhase.procInfo.type))
      {
        //
        // For each match load the reference event phase waveforms
        //
        TimeWindow tw = xcorrTimeWindowLong(xcorrOpt, refPhase);

        for (const string &component : components)
        {
          // This doesn't really load the trace but force the request to reach
          // the loader, which will load all the traces later
          getWaveform(*proc, tw, refEv, refPhase, component);
        }
      }
    }
  }

  if (onlyCatalogPhases)
  {
    return nullptr;
  }

  // This will actually dowanload the waveworms all at once
  batchLoader->load();

  // print counters
  WfCounters wfcount{};
  wfcount.update(batchLoader.get());
  wfcount.update(diskLoader.get());
  logInfoF("Event %s: waveforms downloaded %u, not available %u, loaded from "
           "disk cache %u",
           string(refEv).c_str(), wfcount.downloaded, wfcount.no_avail,
           wfcount.disk_cached);

  return shared_ptr<Waveform::Processor>(new Waveform::MemCachedProc(proc));
}

TimeWindow DD::xcorrTimeWindowLong(const XcorrOptions &xcorrOpt,
                                   const Phase &phase) const
{
  const auto &xcorrCfg = xcorrOpt.phase.at(phase.procInfo.type);
  TimeWindow tw        = xcorrTimeWindowShort(xcorrOpt, phase);
  tw.setStartTime(tw.startTime() - secToDur(xcorrCfg.maxDelay));
  tw.setEndTime(tw.endTime() + secToDur(xcorrCfg.maxDelay));
  return tw;
}

TimeWindow DD::xcorrTimeWindowShort(const XcorrOptions &xcorrOpt,
                                    const Phase &phase) const
{
  const auto &xcorrCfg = xcorrOpt.phase.at(phase.procInfo.type);
  UTCTime::duration shortDuration =
      secToDur(xcorrCfg.endOffset - xcorrCfg.startOffset);
  UTCTime::duration shortTimeCorrection = secToDur(xcorrCfg.startOffset);
  return TimeWindow(phase.time + shortTimeCorrection, shortDuration);
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
  if (phase1.procInfo.type != phase2.procInfo.type)
  {
    logErrorF("Skipping cross-correlation: mismatching phases (%s and %s)",
              string(phase1).c_str(), string(phase2).c_str());
    return false;
  }

  //
  // Try to use the same channels for the cross-correlation. In case the two
  // phases differ, do not change the catalog phase channels.
  //
  const string channelCodeRoot1 = getBandAndInstrumentCodes(phase1.channelCode);
  const string channelCodeRoot2 = getBandAndInstrumentCodes(phase2.channelCode);

  if (channelCodeRoot1 != channelCodeRoot2)
  {
    bool channelsAreCompatible = false;
    for (const pair<string, string> &p : _cfg.compatibleChannels)
    {
      if ((p.first == channelCodeRoot1 && p.second == channelCodeRoot2) ||
          (p.second == channelCodeRoot1 && p.first == channelCodeRoot2))
      {
        channelsAreCompatible = true;
        break;
      }
    }

    if (!channelsAreCompatible)
    {
      logDebugF("Skipping cross-correlation: incompatible channels %s and %s "
                "(%s and %s)",
                channelCodeRoot1.c_str(), channelCodeRoot2.c_str(),
                string(phase1).c_str(), string(phase2).c_str());
      return false;
    }
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

  TimeWindow tw1 = xcorrTimeWindowLong(xcorrOpt, phase1);
  TimeWindow tw2 = xcorrTimeWindowLong(xcorrOpt, phase2);

  // Load the long `tr1`, because we want to cache the long version. Then
  // we'll trim it.
  shared_ptr<const Trace> tr1 =
      getWaveform(ph1Cache, tw1, event1, phase1, component);
  if (!tr1)
  {
    logDebugF("Skipping cross-correlation: no waveform SNR for phase %s",
              string(phase1).c_str());
    return false;
  }

  // Load the long `tr2`, because we want to cache the long version. Then
  // we'll trim it.
  shared_ptr<const Trace> tr2 =
      getWaveform(ph2Cache, tw2, event2, phase2, component);
  if (!tr2)
  {
    logDebugF("Skipping cross-correlation: no waveform for phase %s",
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
    TimeWindow tw2Short = xcorrTimeWindowShort(xcorrOpt, phase2);
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
    TimeWindow tw1Short = xcorrTimeWindowShort(xcorrOpt, phase1);
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
  else if (component == "H")
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
      const string &stationId              = std::get<0>(t);
      const Catalog::Phase::Type phaseType = std::get<1>(t);
      const unsigned neighEvId             = std::get<2>(t);

      if (!xcorr.has(neighbours.referenceId(), neighEvId, stationId, phaseType))
      {
        continue;
      }
      const auto &e =
          xcorr.get(neighbours.referenceId(), neighEvId, stationId, phaseType);
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

  unsigned performed      = 0;
  unsigned failed         = 0;
  unsigned aboveThreshold = 0;
  logInfo("Station P  #CC   fail   above-threshold     corr-coeff five-number "
          "summary");
  for (const auto &kv : pCountByStation)
  {
    const string &stationId = kv.first;
    const Counters &c       = kv.second;

    double min, max, q1, q2, q3;
    compute5numberSummary(c.coeff, min, max, q1, q2, q3);

    unsigned tot = c.performed + c.failed;

    logInfoF("%-12s %10u %.1f%% %.1f%% %.2f %.2f %.2f %.2f %.2f",
             stationId.c_str(), tot, (c.failed * 100.0 / tot),
             (c.aboveThreshold * 100.0 / tot), min, q1, q2, q3, max);
    performed += c.performed;
    failed += c.failed;
    aboveThreshold += c.aboveThreshold;
  }
  unsigned tot = performed + failed;
  logInfoF("Total %10u %.1f%% %.1f%%", tot, (failed * 100.0 / tot),
           (aboveThreshold * 100.0 / tot));

  performed      = 0;
  failed         = 0;
  aboveThreshold = 0;
  logInfo("Station S  #CC   fail   above-threshold     corr-coeff five-number "
          "summary");
  for (const auto &kv : sCountByStation)
  {
    const string &stationId = kv.first;
    const Counters &c       = kv.second;

    double min, max, q1, q2, q3;
    compute5numberSummary(c.coeff, min, max, q1, q2, q3);

    unsigned tot = c.performed + c.failed;

    logInfoF("%-12s %10u %.1f%% %.1f%% %.2f %.2f %.2f %.2f %.2f",
             stationId.c_str(), tot, (c.failed * 100.0 / tot),
             (c.aboveThreshold * 100.0 / tot), min, q1, q2, q3, max);
    performed += c.performed;
    failed += c.failed;
    aboveThreshold += c.aboveThreshold;
  }
  tot = performed + failed;
  logInfoF("Total %10u %.1f%% %.1f%%", tot, (failed * 100.0 / tot),
           (aboveThreshold * 100.0 / tot));
}

} // namespace HDD

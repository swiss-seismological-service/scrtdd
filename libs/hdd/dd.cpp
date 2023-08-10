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
using std::chrono::duration;
using Event     = HDD::Catalog::Event;
using Phase     = HDD::Catalog::Phase;
using Station   = HDD::Catalog::Station;
using Transform = HDD::Waveform::Processor::Transform;
using AQ_ACTION = HDD::SolverOptions::AQ_ACTION;
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

} // namespace

namespace HDD {

DD::DD(const Catalog &catalog,
       const Config &cfg,
       const std::shared_ptr<HDD::TravelTimeTable> &ttt,
       const std::shared_ptr<HDD::Waveform::Proxy> &wf)
    : _cfg(cfg), _srcCat(catalog),
      _bgCat(Catalog::filterPhasesAndSetWeights(_srcCat,
                                                Phase::Source::CATALOG,
                                                _cfg.validPphases,
                                                _cfg.validSphases)),
      _ttt(std::move(ttt)), _wf(std::move(wf))
{
  disableSaveProcessing();
  disableCatalogWaveformDiskCache();
  disableAllWaveformDiskCache();
}

void DD::disableSaveProcessing() { _saveProcessing = false; }

void DD::enableSaveProcessing(const string &workingDir)
{
  _saveProcessing = true;
  _workingDir     = workingDir;
  if (!pathExists(_workingDir))
  {
    if (!createDirectories(_workingDir))
    {
      string msg = "Unable to create working directory: " + _workingDir;
      throw Exception(msg);
    }
  }
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
  _wfAccess.loader.reset(new Waveform::BasicLoader(_wf));
  _wfAccess.diskCache = nullptr;
  _wfAccess.extraLen  = nullptr;
  _wfAccess.snrFilter = nullptr;
  _wfAccess.memCache  = nullptr;

  shared_ptr<Waveform::Loader> currLdr = _wfAccess.loader;

  if (_useCatalogWaveformDiskCache)
  {
    _wfAccess.diskCache.reset(
        new Waveform::DiskCachedLoader(_wf, currLdr, _cacheDir));
    _wfAccess.extraLen.reset(new Waveform::ExtraLenLoader(
        _wfAccess.diskCache, _cfg.diskTraceMinLen));
    currLdr = _wfAccess.extraLen;
  }

  shared_ptr<Waveform::Processor> currProc(
      new Waveform::BasicProcessor(_wf, currLdr, _cfg.wfFilter.extraTraceLen));

  if (_cfg.snr.minSnr > 0)
  {
    _wfAccess.snrFilter.reset(new Waveform::SnrFilterPrc(
        currProc, _cfg.snr.minSnr, _cfg.snr.noiseStart, _cfg.snr.noiseEnd,
        _cfg.snr.signalStart, _cfg.snr.signalEnd));
    currProc = _wfAccess.snrFilter;
  }

  _wfAccess.memCache.reset(new Waveform::MemCachedProc(currProc));
}

void DD::replaceWaveformCacheLoader(const shared_ptr<Waveform::Loader> &baseLdr)
{
  if (_useCatalogWaveformDiskCache)
  {
    _wfAccess.diskCache->setAuxLoader(baseLdr);
  }
  else if (_cfg.snr.minSnr > 0)
  {
    _wfAccess.snrFilter->setAuxProcessor(
        shared_ptr<Waveform::Processor>(new Waveform::BasicProcessor(
            _wf, baseLdr, _cfg.wfFilter.extraTraceLen)));
  }
  else
  {
    _wfAccess.memCache->setAuxProcessor(
        shared_ptr<Waveform::Processor>(new Waveform::BasicProcessor(
            _wf, baseLdr, _cfg.wfFilter.extraTraceLen)));
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

const std::vector<std::string> DD::xcorrComponents(const Phase &phase) const
{
  const auto xcorrCfg = _cfg.xcorr.at(phase.procInfo.type);
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

void DD::preloadWaveforms()
{
  //
  // preload waveforms, store them on disk and cache them in memory
  // (already processed).
  // For better performance we want to load waveforms in batch
  // (batchloader)
  //
  logInfo("Loading catalog waveform data (%lu events to load)",
          _bgCat.getEvents().size());

  auto forEventWaveforms =
      [this](const Event &event, unsigned &numPhases, unsigned &numSPhases,
             std::function<bool(const TimeWindow &, const Event &,
                                const Phase &, const string &)> func) {
        auto eqlrng = _bgCat.getPhases().equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
          const Phase &phase = it->second;
          TimeWindow tw      = xcorrTimeWindowLong(phase);

          for (const string &component : xcorrComponents(phase))
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

    logDebug("Loading event %u waveforms...", event.id);

    shared_ptr<Waveform::BatchLoader> batchLoader(
        new Waveform::BatchLoader(_wf));
    replaceWaveformCacheLoader(batchLoader);

    // the following doesn't load anything unless the waveform is not alreday on
    // the disk cache, instead it records in BatchLoader the waveforms we want
    // to download
    auto requestWaveform = [this](const TimeWindow &tw, const Event &ev,
                                  const Phase &ph, const string &component) {
      getWaveform(*_wfAccess.memCache, tw, ev, ph, component);
      return false;
    };

    forEventWaveforms(event, numPhases, numSPhases, requestWaveform);

    // This will actually dowanload the waveworms
    batchLoader->load();

    // now process and store the waveforms into memory
    auto cacheWaveform = [this](const TimeWindow &tw, const Event &ev,
                                const Phase &ph, const string &component) {
      return getWaveform(*_wfAccess.memCache, tw, ev, ph, component) != nullptr;
    };

    unsigned discard;
    forEventWaveforms(event, discard, discard, cacheWaveform);

    wfcount.update(batchLoader.get());

    const unsigned onePercent =
        _bgCat.getEvents().size() < 100 ? 1 : (_bgCat.getEvents().size() / 100);
    if (++numEvents % onePercent == 0)
    {
      logInfo("Loaded %.1f%% of catalog phase waveforms",
              (numEvents * 100.0 / _bgCat.getEvents().size()));
    }
  }

  // restore proper loader
  replaceWaveformCacheLoader(
      shared_ptr<Waveform::Loader>(new Waveform::BasicLoader(_wf)));

  wfcount.update(_wfAccess.loader.get());
  wfcount.update(_wfAccess.diskCache.get());

  logInfo("Finished preloading catalog waveform data: total events %lu total "
          "phases %u (P %.f%%, S %.f%%). Waveforms downloaded %u, not "
          "available %u, loaded from disk cache %u",
          _bgCat.getEvents().size(), numPhases,
          ((numPhases - numSPhases) * 100. / numPhases),
          (numSPhases * 100. / numPhases), wfcount.downloaded, wfcount.no_avail,
          wfcount.disk_cached);
}

void DD::dumpWaveforms(const string &basePath)
{
  if (!basePath.empty() && !pathExists(basePath))
  {
    if (!createDirectories(basePath))
    {
      throw Exception("Unable to create waveform dump directory: " + basePath);
    }
  }

  auto waveformPath = [](const string &wfDebugDir, const Catalog::Event &ev,
                         const Catalog::Phase &ph, const string &channelCode,
                         const string &ext) -> string {
    string debugFile =
        HDD::strf("ev%u.%s.%s.%s.%s.%s.%s.mseed", ev.id, ph.networkCode.c_str(),
                  ph.stationCode.c_str(), ph.locationCode.c_str(),
                  channelCode.c_str(), ph.type.c_str(), ext.c_str());
    return HDD::joinPath(wfDebugDir, debugFile);
  };

  for (const auto &kv : _bgCat.getEvents())
  {
    const Event &event = kv.second;
    auto eqlrng        = _bgCat.getPhases().equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      const Phase &phase = it->second;
      TimeWindow tw      = xcorrTimeWindowLong(phase);

      for (const string &component : xcorrComponents(phase))
      {
        // getAuxProcessor() -> we don't really want to keep the waveforms in
        // memory
        shared_ptr<const Trace> tr =
            getWaveform(*_wfAccess.memCache->getAuxProcessor(), tw, event,
                        phase, component);

        // check low SNR
        bool lowSNR = false;
        if (!tr && _wfAccess.snrFilter)
        {
          _wfAccess.snrFilter->setEnabled(false);
          tr = getWaveform(*_wfAccess.memCache->getAuxProcessor(), tw, event,
                           phase, component);

          if (tr) lowSNR = true;
          _wfAccess.snrFilter->setEnabled(true);
        }

        if (!tr) continue;

        string channelCode =
            getBandAndInstrumentCodes(phase.channelCode) + component;
        string ext =
            (phase.procInfo.source == Catalog::Phase::Source::THEORETICAL)
                ? "xcorrDetected"
                : (phase.isManual ? "manual" : "automatic");
        if (lowSNR) ext += ".lowSNR";

        const string wfPath =
            waveformPath(basePath, event, phase, channelCode, ext);

        logInfo("Writing %s", wfPath.c_str());
        try
        {

          _wf->writeTrace(*tr, wfPath);
        }
        catch (exception &e)
        {
          logWarning("Couldn't write waveform to disk %s: %s", wfPath.c_str(),
                     e.what());
        }
      }
    }
  }
}

std::list<Catalog> DD::findClusters(const ClusteringOptions &clustOpt)
{
  // find Neighbours for each event in the catalog
  unordered_map<unsigned, unique_ptr<Neighbours>> neighboursByEvent =
      selectNeighbouringEventsCatalog(
          _bgCat, clustOpt.minWeight, clustOpt.minESdist, clustOpt.maxESdist,
          clustOpt.minEStoIEratio, clustOpt.minDTperEvt, clustOpt.maxDTperEvt,
          clustOpt.minNumNeigh, clustOpt.maxNumNeigh, clustOpt.numEllipsoids,
          clustOpt.maxEllipsoidSize, true);

  // Organize the neighbours by not connected clusters
  list<unordered_map<unsigned, unique_ptr<Neighbours>>> clusters =
      clusterizeNeighbouringEvents(neighboursByEvent);

  std::list<Catalog> catalogs;
  for (const auto &neighCluster : clusters)
  {
    Catalog cat;
    for (const auto &kv : neighCluster) cat.add(kv.first, _bgCat, true);
    catalogs.push_back(std::move(cat));
  }
  return catalogs;
}

unique_ptr<Catalog> DD::relocateMultiEvents(const ClusteringOptions &clustOpt,
                                            const SolverOptions &solverOpt,
                                            XCorrCache &precomputed)
{
  logInfo("Starting DD relocator in multiple events mode");

  Catalog catToReloc(_bgCat);

  // prepare a folder for debug files
  string catalogWorkingDir;
  do
  {
    catalogWorkingDir = generateWorkingSubDir("multievent");
    catalogWorkingDir = joinPath(_workingDir, catalogWorkingDir);
  } while (pathExists(catalogWorkingDir));

  if (_saveProcessing)
  {
    if (!pathExists(catalogWorkingDir))
    {
      if (!createDirectories(catalogWorkingDir))
      {
        string msg = "Unable to create working directory: " + catalogWorkingDir;
        throw Exception(msg);
      }
    }
    logInfo("Working dir %s", catalogWorkingDir.c_str());
  }

  // prepare file logger
  unique_ptr<Logger::File> fileLogger;
  if (_saveProcessing)
  {
    string logFile = joinPath(catalogWorkingDir, "info.log");
    fileLogger =
        Logger::toFile(logFile, {Logger::Level::info, Logger::Level::warning,
                                 Logger::Level::error});
  }

  // find Neighbours for each event in the catalog
  unordered_map<unsigned, unique_ptr<Neighbours>> neighboursByEvent =
      selectNeighbouringEventsCatalog(
          catToReloc, clustOpt.minWeight, clustOpt.minESdist,
          clustOpt.maxESdist, clustOpt.minEStoIEratio, clustOpt.minDTperEvt,
          clustOpt.maxDTperEvt, clustOpt.minNumNeigh, clustOpt.maxNumNeigh,
          clustOpt.numEllipsoids, clustOpt.maxEllipsoidSize, true);

  // Organize the neighbours by not connected clusters
  list<unordered_map<unsigned, unique_ptr<Neighbours>>> clusters =
      clusterizeNeighbouringEvents(neighboursByEvent);

  logInfo("Found %lu event clusters", clusters.size());

  if (_saveProcessing)
  {
    catToReloc.writeToFile(joinPath(catalogWorkingDir, "starting-event.csv"),
                           joinPath(catalogWorkingDir, "starting-phase.csv"),
                           joinPath(catalogWorkingDir, "starting-station.csv"));
  }

  //
  // relocate one cluster a time
  //
  unique_ptr<Catalog> relocatedCatalog(new Catalog());

  unsigned clusterId = 1;
  for (const auto &neighCluster : clusters)
  {
    logInfo("Relocating cluster %u (%lu events)", clusterId,
            neighCluster.size());

    if (_saveProcessing)
    {
      Catalog catToDump;
      for (const auto &kv : neighCluster)
        catToDump.add(kv.first, catToReloc, true);
      string prefix = strf("cluster-%u", clusterId);
      catToDump.writeToFile(
          joinPath(catalogWorkingDir, (prefix + "-event.csv")),
          joinPath(catalogWorkingDir, (prefix + "-phase.csv")),
          joinPath(catalogWorkingDir, (prefix + "-station.csv")));
    }

    // Perform cross-correlation which also detects picks around theoretical
    // arrival times. The catalog will be updated with the corresponding
    // theoretical phases.
    const XCorrCache xcorr = buildXCorrCache(
        catToReloc, neighCluster, false, clustOpt.xcorrMaxEvStaDist,
        clustOpt.xcorrMaxInterEvDist, precomputed);

    // the actual relocation
    unique_ptr<Catalog> relocatedCluster =
        relocate(catToReloc, neighCluster, solverOpt, false, xcorr);

    relocatedCatalog->add(*relocatedCluster, true);

    if (_saveProcessing)
    {
      string prefix = strf("relocated-cluster-%u", clusterId);
      relocatedCluster->writeToFile(
          joinPath(catalogWorkingDir, (prefix + "-event.csv")),
          joinPath(catalogWorkingDir, (prefix + "-phase.csv")),
          joinPath(catalogWorkingDir, (prefix + "-station.csv")));
    }
    clusterId++;
    precomputed = std::move(xcorr);
  }

  if (_saveProcessing)
  {
    relocatedCatalog->writeToFile(
        joinPath(catalogWorkingDir, "relocated-event.csv"),
        joinPath(catalogWorkingDir, "relocated-phase.csv"),
        joinPath(catalogWorkingDir, "relocated-station.csv"));
    writeXCorrToFile(precomputed, _bgCat,
                     joinPath(catalogWorkingDir, "xcorr.csv"));
  }

  _ttt->freeResources();

  return relocatedCatalog;
}

unique_ptr<Catalog> DD::relocateSingleEvent(const Catalog &singleEvent,
                                            const ClusteringOptions &clustOpt1,
                                            const ClusteringOptions &clustOpt2,
                                            const SolverOptions &solverOpt)
{
  const Catalog &bgCat = _bgCat;

  // there must be only one event in the catalog, the origin to relocate
  const Event &evToRelocate = singleEvent.getEvents().begin()->second;
  const auto &evToRelocatePhases =
      singleEvent.getPhases().equal_range(evToRelocate.id);

  logInfo("Starting DD relocator in single event mode: event %s lat %.6f lon "
          "%.6f depth %.4f mag %.2f time %s #phases %ld",
          string(evToRelocate).c_str(), evToRelocate.latitude,
          evToRelocate.longitude, evToRelocate.depth, evToRelocate.magnitude,
          UTCClock::toString(evToRelocate.time).c_str(),
          std::distance(evToRelocatePhases.first, evToRelocatePhases.second));

  // prepare a folder for debug files
  string baseWorkingDir;
  if (_saveProcessing)
  {
    do
    {
      baseWorkingDir = generateWorkingSubDir(evToRelocate);
      baseWorkingDir = joinPath(_workingDir, baseWorkingDir);
    } while (pathExists(baseWorkingDir));

    if (!createDirectories(baseWorkingDir))
    {
      string msg = "Unable to create working directory: " + baseWorkingDir;
      throw Exception(msg);
    }
    logInfo("Working dir %s", baseWorkingDir.c_str());
  }

  // prepare file logger
  unique_ptr<Logger::File> fileLogger;
  if (_saveProcessing)
  {
    string logFile = joinPath(baseWorkingDir, "info.log");
    fileLogger =
        Logger::toFile(logFile, {Logger::Level::info, Logger::Level::warning,
                                 Logger::Level::error});
  }

  logInfo("Performing step 1: initial location refinement (no "
          "cross-correlation)");

  string eventWorkingDir = joinPath(baseWorkingDir, "step1");

  Catalog evToRelocateCat =
      Catalog::filterPhasesAndSetWeights(singleEvent, Phase::Source::RT_EVENT,
                                         _cfg.validPphases, _cfg.validSphases);

  unique_ptr<Catalog> relocatedEvCat = relocateEventSingleStep(
      bgCat, evToRelocateCat, eventWorkingDir, clustOpt1, solverOpt, false);

  if (relocatedEvCat)
  {
    const Event &ev = relocatedEvCat->getEvents().begin()->second;
    logInfo("Step 1 relocation successful, new location: "
            "lat %.6f lon %.6f depth %.4f time %s",
            ev.latitude, ev.longitude, ev.depth,
            UTCClock::toString(ev.time).c_str());
    logInfo("Relocation report: %s", relocationReport(*relocatedEvCat).c_str());

    evToRelocateCat = std::move(*relocatedEvCat);
  }
  else
  {
    logError("Failed to perform step 1 origin relocation");
  }

  logInfo("Performing step 2: relocation with cross-correlation");

  eventWorkingDir = joinPath(baseWorkingDir, "step2");

  unique_ptr<Catalog> relocatedEvWithXcorr = relocateEventSingleStep(
      bgCat, evToRelocateCat, eventWorkingDir, clustOpt2, solverOpt, true);

  if (relocatedEvWithXcorr)
  {
    Event ev = relocatedEvWithXcorr->getEvents().begin()->second;
    logInfo("Step 2 relocation successful, new location: "
            "lat %.6f lon %.6f depth %.4f time %s",
            ev.latitude, ev.longitude, ev.depth,
            UTCClock::toString(ev.time).c_str());
    logInfo("Relocation report: %s",
            relocationReport(*relocatedEvWithXcorr).c_str());

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

    logInfo("Total Changes: location=%.3f[km] depth=%.3f[km] "
            "time=%.4f[sec] Rms=%.4f[sec] (before/after %.4f/%.4f)",
            ev.relocInfo.locChange, ev.relocInfo.depthChange,
            ev.relocInfo.timeChange,
            (ev.relocInfo.finalRms - ev.relocInfo.startRms),
            ev.relocInfo.startRms, ev.relocInfo.finalRms);
  }
  else
  {
    logError("Failed to perform step 2 origin relocation");
  }

  _ttt->freeResources();

  if (!relocatedEvWithXcorr) throw Exception("Failed origin relocation");

  return relocatedEvWithXcorr;
}

unique_ptr<Catalog>
DD::relocateEventSingleStep(const Catalog &bgCat,
                            const Catalog &evToRelocateCat,
                            const string &workingDir,
                            const ClusteringOptions &clustOpt,
                            const SolverOptions &solverOpt,
                            bool doXcorr)
{
  if (_saveProcessing)
  {
    if (!createDirectories(workingDir))
    {
      string msg = "Unable to create working directory: " + workingDir;
      throw Exception(msg);
    }
    logInfo("Working dir %s", workingDir.c_str());

    evToRelocateCat.writeToFile(
        joinPath(workingDir, "single-event.csv"),
        joinPath(workingDir, "single-event-phase.csv"),
        joinPath(workingDir, "single-event-station.csv"));
  }

  unique_ptr<Catalog> relocatedEvCat;

  try
  {
    // extract event to relocate
    const Event &evToRelocate = evToRelocateCat.getEvents().begin()->second;

    //
    // select neighbouring events
    //
    bool keepUnmatchedPhases = doXcorr; // useful for detecting missed picks

    unique_ptr<Neighbours> neighbours = selectNeighbouringEvents(
        bgCat, evToRelocate, evToRelocateCat, clustOpt.minWeight,
        clustOpt.minESdist, clustOpt.maxESdist, clustOpt.minEStoIEratio,
        clustOpt.minDTperEvt, clustOpt.maxDTperEvt, clustOpt.minNumNeigh,
        clustOpt.maxNumNeigh, clustOpt.numEllipsoids, clustOpt.maxEllipsoidSize,
        keepUnmatchedPhases);

    logInfo("Found %zu neighbouring events", neighbours->ids.size());

    //
    // prepare catalog to relocate
    //
    unique_ptr<Catalog> catalog = neighbours->toCatalog(bgCat);
    unsigned evToRelocateNewId =
        catalog->add(evToRelocate.id, evToRelocateCat, false);
    neighbours->refEvId = evToRelocateNewId;

    if (_saveProcessing)
    {
      catalog->writeToFile(joinPath(workingDir, "starting-event.csv"),
                           joinPath(workingDir, "starting-phase.csv"),
                           joinPath(workingDir, "starting-station.csv"));
    }

    unordered_map<unsigned, unique_ptr<Neighbours>> neighCluster;
    neighCluster.emplace(neighbours->refEvId, std::move(neighbours));

    XCorrCache xcorr;
    if (doXcorr)
    {
      // Perform cross-correlation, which also detects picks around
      // theoretical arrival times. The catalog will be updated with the
      // corresponding phases.
      xcorr = buildXCorrCache(
          *catalog, neighCluster, clustOpt.xcorrDetectMissingPhases,
          clustOpt.xcorrMaxEvStaDist, clustOpt.xcorrMaxInterEvDist);
    }

    // the actual relocation
    relocatedEvCat = relocate(*catalog, neighCluster, solverOpt, true, xcorr);

    if (_saveProcessing)
    {
      relocatedEvCat->writeToFile(
          joinPath(workingDir, "relocated-event.csv"),
          joinPath(workingDir, "relocated-phase.csv"),
          joinPath(workingDir, "relocated-station.csv"));
    }
  }
  catch (exception &e)
  {
    logError("%s", e.what());
  }

  return relocatedEvCat;
}

unique_ptr<Catalog> DD::relocate(
    const Catalog &catalog,
    const unordered_map<unsigned, unique_ptr<Neighbours>> &neighCluster,
    const SolverOptions &solverOpt,
    bool keepNeighboursFixed,
    const XCorrCache &xcorr) const
{
  logInfo("Building and solving double-difference system...");

  //
  // iterate the solver computation multiple times
  //
  unique_ptr<const Catalog> finalCatalog(new Catalog(catalog));
  unordered_map<unsigned, unique_ptr<Neighbours>> finalNeighCluster;
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

    logInfo("Solving iteration %u num events %lu. Parameters: "
            "observWeight TT/CC=%.2f/%.2f dampingFactor=%.2f "
            "downWeightingByResidual=%.2f absLocConstraint=%.2f",
            iteration, neighCluster.size(), absTTDiffObsWeight, xcorrObsWeight,
            dampingFactor, downWeightingByResidual, absLocConstraint);

    // create a solver and then add observations
    Solver solver(solverOpt.type);

    //
    // Add absolute travel time/cross-correlation differences to the solver
    // (the observations)
    //
    for (const auto &kv : neighCluster)
    {
      addObservations(solver, absTTDiffObsWeight, xcorrObsWeight, *finalCatalog,
                      *kv.second, keepNeighboursFixed,
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
      logInfo("Cannot solve the double-difference system, stop here (%s)",
              e.what());
      break;
    }

    // prepare for next iteration
    obsparams = ObservationParams();

    // update event parameters
    finalCatalog = updateRelocatedEvents(
        solver, *finalCatalog, solverOpt, neighCluster, obsparams,
        std::max(absTTDiffObsWeight, xcorrObsWeight), finalNeighCluster);
  }

  // compute last bit of statistics for the relocated events
  return updateRelocatedEventsFinalStats(catalog, *finalCatalog,
                                         finalNeighCluster);
}

string DD::relocationReport(const Catalog &relocatedEv)
{
  const Event &event = relocatedEv.getEvents().begin()->second;
  if (!event.relocInfo.isRelocated) return "Event not relocated";

  return strf("Origin changes: location=%.3f[km] depth=%.3f[km] time=%.4f[sec] "
              "Rms change [sec]: %.4f (before/after %.4f/%.4f) "
              "Neighbours=%u Used Phases: P=%u S=%u "
              "Stations distance [km]: min=%.1f median=%.1f max=%.1f "
              "DD observations: %u (CC P/S %u/%u TT P/S %u/%u) "
              "DD residuals [msec]: before=%.f+/-%.1f after=%.f+/-%.1f",
              event.relocInfo.locChange, event.relocInfo.depthChange,
              event.relocInfo.timeChange,
              (event.relocInfo.finalRms - event.relocInfo.startRms),
              event.relocInfo.startRms, event.relocInfo.finalRms,
              event.relocInfo.numNeighbours, event.relocInfo.phases.usedP,
              event.relocInfo.phases.usedS,
              event.relocInfo.phases.stationDistMin,
              event.relocInfo.phases.stationDistMedian,
              event.relocInfo.phases.stationDistMax,
              (event.relocInfo.dd.numCCp + event.relocInfo.dd.numCCs +
               event.relocInfo.dd.numTTp + event.relocInfo.dd.numTTs),
              event.relocInfo.dd.numCCp, event.relocInfo.dd.numCCs,
              event.relocInfo.dd.numTTp, event.relocInfo.dd.numTTs,
              event.relocInfo.dd.startResidualMedian * 1000,
              event.relocInfo.dd.startResidualMAD * 1000,
              event.relocInfo.dd.finalResidualMedian * 1000,
              event.relocInfo.dd.finalResidualMAD * 1000);
}

/*
 * Add both the absolute travel time differences and the differential
 * travel times from the cross-correlation for pairs of earthquakes to the
 * solver.
 */
void DD::addObservations(Solver &solver,
                         double absTTDiffObsWeight,
                         double xcorrObsWeight,
                         const Catalog &catalog,
                         const Neighbours &neighbours,
                         bool keepNeighboursFixed,
                         bool usePickUncertainty,
                         const XCorrCache &xcorr,
                         ObservationParams &obsparams) const
{
  // copy event because we'll update it
  const Event &refEv = catalog.getEvents().at(neighbours.refEvId);

  //
  // loop through reference event phases
  //
  auto eqlrng = catalog.getPhases().equal_range(refEv.id);
  for (auto it = eqlrng.first; it != eqlrng.second; ++it)
  {
    const Phase &refPhase  = it->second;
    const Station &station = catalog.getStations().at(refPhase.stationId);
    char phaseTypeAsChar   = static_cast<char>(refPhase.procInfo.type);
    const auto xcorrCfg    = _cfg.xcorr.at(refPhase.procInfo.type);

    //
    // loop through neighbouring events and look for the matching phase
    //
    for (unsigned neighEvId : neighbours.ids)
    {
      const Event &event = catalog.getEvents().at(neighEvId);

      if (!neighbours.has(neighEvId, refPhase.stationId,
                          refPhase.procInfo.type))
        continue;

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
        logDebug("Ignoring phase %s with negative travel time",
                 string(refPhase).c_str());
        continue;
      }

      UTCTime::duration travel_time = phase.time - event.time;
      if (travel_time.count() < 0)
      {
        logDebug("Ignoring phase %s with negative travel time",
                 string(phase).c_str());
        continue;
      }

      if (!obsparams.add(*_ttt, refEv, station, refPhase, true) ||
          !obsparams.add(*_ttt, event, station, phase, !keepNeighboursFixed))
      {
        logDebug("Skipping observation (ev %u-%u sta %s phase %c)", refEv.id,
                 event.id, station.id.c_str(), phaseTypeAsChar);
        continue;
      }

      //
      // compute absolute travel time differences to the solver
      //
      duration<double> diffTime = ref_travel_time - travel_time;
      double weight =
          usePickUncertainty
              ? (refPhase.procInfo.weight + phase.procInfo.weight) / 2.0
              : 1.0;
      //
      // Check if we have cross-correlation results for the current
      // event/refEv pair at station/phase and refine the differential
      // time
      //
      bool isXcorr = false;

      if (xcorr.has(refEv.id, event.id, refPhase.stationId,
                    refPhase.procInfo.type))
      {
        const auto &e = xcorr.get(refEv.id, event.id, refPhase.stationId,
                                  refPhase.procInfo.type);

        if (e.valid && e.coeff >= xcorrCfg.minCoef)
        {
          isXcorr = true;
          diffTime -= duration<double>(e.lag);
          weight *= xcorrObsWeight;
        }
      }

      if (!isXcorr)
      {
        weight *= absTTDiffObsWeight;
      }

      solver.addObservation(refEv.id, event.id, refPhase.stationId,
                            phaseTypeAsChar, diffTime.count(), weight, isXcorr);
    }
  }
}

bool DD::ObservationParams::add(HDD::TravelTimeTable &ttt,
                                const Event &event,
                                const Station &station,
                                const Phase &phase,
                                bool computeEvChanges)
{
  char phaseType = static_cast<char>(phase.procInfo.type);
  const string key =
      std::to_string(event.id) + "@" + station.id + ":" + phaseType;
  if (_entries.find(key) == _entries.end())
  {
    try
    {
      double travelTime, azimuth, takeOffAngle, velocityAtSrc;
      ttt.compute(event, station, string(1, phaseType), travelTime, azimuth,
                  takeOffAngle, velocityAtSrc);
      double ttResidual = travelTime - durToSec(phase.time - event.time);
      _entries[key]     = Entry{event,        station,       phaseType,
                            travelTime,   ttResidual,    azimuth,
                            takeOffAngle, velocityAtSrc, computeEvChanges};
    }
    catch (exception &e)
    {
      logWarning(
          "Travel Time Table error: %s (Event lat %.6f lon %.6f depth %.6f "
          "Station lat %.6f lon %.6f elevation %.f )",
          e.what(), event.latitude, event.longitude, event.depth,
          station.latitude, station.longitude, station.elevation);
      return false;
    }
  }
  return true;
}

const DD::ObservationParams::Entry &DD::ObservationParams::get(
    unsigned eventId, const string stationId, char phaseType) const
{
  const string key =
      std::to_string(eventId) + "@" + stationId + ":" + phaseType;
  return _entries.at(key);
}

void DD::ObservationParams::addToSolver(Solver &solver) const
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

unique_ptr<Catalog> DD::updateRelocatedEvents(
    const Solver &solver,
    const Catalog &catalog,
    const SolverOptions &solverOpt,
    const unordered_map<unsigned, unique_ptr<Neighbours>> &neighCluster,
    ObservationParams &obsparams,
    double pickWeightScaler,
    unordered_map<unsigned, unique_ptr<Neighbours>> &finalNeighCluster // output
) const
{
  unordered_map<string, Station> stations    = catalog.getStations();
  map<unsigned, Event> events                = catalog.getEvents();
  unordered_multimap<unsigned, Phase> phases = catalog.getPhases();
  unsigned relocatedEvs                      = 0;
  vector<double> allRms;

  //
  // loop through each event with its corresponding neighbours cluster
  //
  for (const auto &kv : neighCluster)
  {
    Event &event = events.at(kv.second->refEvId);

    //
    // get relocation changes (km and sec) computed by the solver for
    // the current event
    //
    double deltaLat, deltaLon, deltaDepth, deltaTime;
    if (!solver.getEventChanges(event.id, deltaLat, deltaLon, deltaDepth,
                                deltaTime))
    {
      event.relocInfo.isRelocated = false;
      continue;
    }

    // check for airquakes
    if (solverOpt.airQuakes.action != AQ_ACTION::NONE &&
        (event.depth + deltaDepth) <
            -(solverOpt.airQuakes.elevationThreshold / 1000.))
    {
      if (solverOpt.airQuakes.action == AQ_ACTION::RESET)
      {
        logDebug("Revert airquake event %s to previous location",
                 string(event).c_str());
        continue;
      }
      else if (solverOpt.airQuakes.action == AQ_ACTION::RESET_DEPTH)
      {
        logDebug("Revert airquake event %s to previous depth",
                 string(event).c_str());
        deltaDepth = 0;
      }
    }

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

    unique_ptr<Neighbours> finalNeighbours(new Neighbours());
    finalNeighbours->refEvId = event.id;
    bool isFirstIteration    = !event.relocInfo.isRelocated;

    //
    // update event location/time and compute statistics
    //
    relocatedEvs++;
    event.latitude  = newLat;
    event.longitude = newLon;
    event.depth += deltaDepth;
    event.time += secToDur(deltaTime);
    event.relocInfo.finalRms      = 0;
    event.relocInfo.isRelocated   = true;
    event.relocInfo.numNeighbours = 0;
    event.relocInfo.phases        = {};
    event.relocInfo.dd.numTTp     = 0;
    event.relocInfo.dd.numTTs     = 0;
    event.relocInfo.dd.numCCp     = 0;
    event.relocInfo.dd.numCCs     = 0;

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
      double meanDDResidual, meanAPrioriWeight, meanFinalWeight;

      if (!solver.getObservationParamsChanges(
              event.id, station.id, phaseTypeAsChar, startTTObs, startCCObs,
              finalTotalObs, meanAPrioriWeight, meanFinalWeight, meanDDResidual,
              neighbourIds))
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
        phase.relocInfo.startMeanDDResidual = meanDDResidual;
      phase.relocInfo.finalMeanDDResidual = meanDDResidual;
      obsResiduals.push_back(meanDDResidual);

      if (obsparams.add(*_ttt, event, station, phase, true))
      {
        double travelTime =
            obsparams.get(event.id, station.id, phaseTypeAsChar).travelTime;
        phase.relocInfo.finalTTResidual =
            travelTime - durToSec(phase.time - event.time);
        rmsCount++;
      }
      else
      {
        phase.relocInfo.finalTTResidual = 0;
      }

      event.relocInfo.finalRms += square(phase.relocInfo.finalTTResidual);
      if (phase.procInfo.type == Phase::Type::P)
      {
        event.relocInfo.phases.usedP++;
        event.relocInfo.dd.numCCp += phase.relocInfo.numCCObs;
        event.relocInfo.dd.numTTp += phase.relocInfo.numTTObs;
      }
      if (phase.procInfo.type == Phase::Type::S)
      {
        event.relocInfo.phases.usedS++;
        event.relocInfo.dd.numCCs += phase.relocInfo.numCCObs;
        event.relocInfo.dd.numTTs += phase.relocInfo.numTTObs;
      }

      for (unsigned nId : neighbourIds)
      {
        finalNeighbours->add(nId, station.id, phase.procInfo.type);
      }
    }

    if (rmsCount > 0)
    {
      event.relocInfo.finalRms = std::sqrt(event.relocInfo.finalRms / rmsCount);
      allRms.push_back(event.relocInfo.finalRms);
    }

    double residualMedian = computeMedian(obsResiduals);
    double residualMAD =
        computeMedianAbsoluteDeviation(obsResiduals, residualMedian);
    if (isFirstIteration)
    {
      event.relocInfo.dd.startResidualMedian = residualMedian;
      event.relocInfo.dd.startResidualMAD    = residualMAD;
    }
    event.relocInfo.dd.finalResidualMedian = residualMedian;
    event.relocInfo.dd.finalResidualMAD    = residualMAD;

    event.relocInfo.numNeighbours               = finalNeighbours->ids.size();
    finalNeighCluster[finalNeighbours->refEvId] = std::move(finalNeighbours);
  }

  const double allRmsMedian = computeMedian(allRms);
  const double allRmsMAD = computeMedianAbsoluteDeviation(allRms, allRmsMedian);

  logInfo("Successfully relocated %u events, RMS median %.3f [msec] median "
          "absolute deviation %.3f [msec]",
          relocatedEvs, allRmsMedian * 1000, allRmsMAD * 1000);

  return unique_ptr<Catalog>(new Catalog(stations, events, phases));
}

unique_ptr<Catalog> DD::updateRelocatedEventsFinalStats(
    const Catalog &startCatalog,
    const Catalog &finalCatalog,
    const unordered_map<unsigned, unique_ptr<Neighbours>> &neighCluster) const
{
  unique_ptr<Catalog> catalogToReturn(new Catalog());
  vector<double> allRms;
  vector<double> stationDist;

  for (const auto &kv : neighCluster)
  {
    const unique_ptr<Neighbours> &neighbours = kv.second;
    auto it = finalCatalog.getEvents().find(neighbours->refEvId);

    // If the event hasn't been relocated, remove it from the final catalog.
    if (it == finalCatalog.getEvents().end() ||
        !it->second.relocInfo.isRelocated)
    {
      continue;
    }

    unique_ptr<Catalog> tmpCat =
        finalCatalog.extractEvent(neighbours->refEvId, true);

    const Event &startEvent = startCatalog.getEvents().at(neighbours->refEvId);
    Event finalEvent        = tmpCat->getEvents().at(neighbours->refEvId);

    //
    // Compute location/time change of start/final origin.
    //
    finalEvent.relocInfo.locChange =
        computeDistance(finalEvent.latitude, finalEvent.longitude,
                        startEvent.latitude, startEvent.longitude);
    finalEvent.relocInfo.depthChange = finalEvent.depth - startEvent.depth;
    finalEvent.relocInfo.timeChange =
        durToSec(finalEvent.time - startEvent.time);

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
        double travelTime = _ttt->compute(
            startEvent, station,
            string(1, static_cast<char>(finalPhase.procInfo.type)));
        double residual =
            travelTime - durToSec(finalPhase.time - startEvent.time);
        finalPhase.relocInfo.startTTResidual = residual;
        tmpCat->updatePhase(finalPhase, false);
        finalEvent.relocInfo.startRms += residual * residual;
        rmsCount++;
      }
      catch (exception &e)
      {
        logWarning(
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
        // Check if this station is actually part of the event phases
        // (remember keepUnmatchedPhases option during clustering). Make sure
        // that the phase was used for relocation.
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

    tmpCat->updateEvent(finalEvent, false);
    catalogToReturn->add(finalEvent.id, *tmpCat, true);
  }

  const double allRmsMedian = computeMedian(allRms);
  const double allRmsMAD = computeMedianAbsoluteDeviation(allRms, allRmsMedian);

  logInfo("Events RMS before relocation: median %.3f [msec] median absolute "
          "deviation %.3f [msec]",
          allRmsMedian * 1000, allRmsMAD * 1000);

  return catalogToReturn;
}

void DD::addMissingEventPhases(const Event &refEv,
                               Catalog &refEvCatalog,
                               const Catalog &searchCatalog,
                               const Neighbours &neighbours)
{
  vector<Phase> newPhases =
      findMissingEventPhases(refEv, refEvCatalog, searchCatalog, neighbours);

  for (Phase &ph : newPhases)
  {
    refEvCatalog.updatePhase(ph, true);
    const Station &station = searchCatalog.getStations().at(ph.stationId);
    refEvCatalog.addStation(station);
  }
}

vector<Phase> DD::findMissingEventPhases(const Event &refEv,
                                         Catalog &refEvCatalog,
                                         const Catalog &searchCatalog,
                                         const Neighbours &neighbours)
{
  //
  // find stations for which the `refEv` doesn't have phases
  //
  vector<MissingStationPhase> missingPhases =
      getMissingPhases(refEv, refEvCatalog, searchCatalog);

  //
  // for each missed phase try to detect it
  //
  vector<Phase> newPhases;
  for (const MissingStationPhase &pair : missingPhases)
  {
    const Station &station      = searchCatalog.getStations().at(pair.first);
    const Phase::Type phaseType = pair.second;

    //
    // Loop through every other event and select the ones who have a manually
    // picked phase for the missing station.
    //
    vector<DD::PhasePeer> peers =
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
      double travelTime      = durToSec(phase.time - event.time);
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

vector<DD::MissingStationPhase>
DD::getMissingPhases(const Event &refEv,
                     Catalog &refEvCatalog,
                     const Catalog &searchCatalog) const
{
  const auto &refEvPhases = refEvCatalog.getPhases().equal_range(refEv.id);

  //
  // Loop through stations and find those for which the `refEv` doesn't have
  // phases.
  //
  vector<MissingStationPhase> missingPhases;
  for (const auto &kv : searchCatalog.getStations())
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

vector<DD::PhasePeer> DD::findPhasePeers(const Station &station,
                                         const Phase::Type &phaseType,
                                         const Catalog &searchCatalog,
                                         const Neighbours &neighbours) const
{
  //
  // Loop through every other event and select those manual phases of the
  // `station` we are interested in.
  //
  vector<PhasePeer> phasePeers;

  for (unsigned neighEvId : neighbours.ids)
  {
    const Event &event = searchCatalog.getEvents().at(neighEvId);

    if (neighbours.has(neighEvId, station.id, phaseType))
    {
      const Phase &phase =
          searchCatalog.searchPhase(neighEvId, station.id, phaseType)->second;

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

Phase DD::createThoreticalPhase(const Station &station,
                                const Phase::Type &phaseType,
                                const Event &refEv,
                                const vector<DD::PhasePeer> &peers,
                                double phaseVelocity)
{
  const auto xcorrCfg = _cfg.xcorr.at(phaseType);

  // store most recent `channelCode` used
  struct
  {
    string channelCode;
    UTCTime time;
  } streamInfo = {"", UTCTime()};

  for (const PhasePeer &peer : peers)
  {
    // const Event& event = peer.first;
    const Phase &phase = peer.second;
    // get the closest stream to the `refEv` (w.r.t. time)
    if (std::abs((refEv.time - phase.time).count()) <
        std::abs((refEv.time - streamInfo.time).count()))
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
      getBandAndInstrumentCodes(streamInfo.channelCode) + "X";
  refEvNewPhase.isManual      = false;
  refEvNewPhase.procInfo.type = phaseType;

  // use phase velocity to compute phase time
  double stationDistance = computeDistance(refEv, station);
  refEvNewPhase.time = refEv.time + secToDur(stationDistance / phaseVelocity);

  refEvNewPhase.lowerUncertainty = Catalog::DEFAULT_AUTOMATIC_PICK_UNCERTAINTY;
  refEvNewPhase.upperUncertainty = Catalog::DEFAULT_AUTOMATIC_PICK_UNCERTAINTY;
  refEvNewPhase.procInfo.weight  = Catalog::computePickWeight(refEvNewPhase);
  refEvNewPhase.procInfo.source  = Phase::Source::THEORETICAL;
  refEvNewPhase.type             = strf("%ct", static_cast<char>(phaseType));

  return refEvNewPhase;
}

XCorrCache DD::buildXCorrCache(
    Catalog &catalog,
    const unordered_map<unsigned, unique_ptr<Neighbours>> &neighCluster,
    bool computeTheoreticalPhases,
    double xcorrMaxEvStaDist,
    double xcorrMaxInterEvDist,
    const XCorrCache &precomputed)
{
  logInfo("Computing differential times via cross-correlation...");

  XCorrCache xcorr(precomputed);
  unsigned long performed = 0;

  for (const auto &kv : neighCluster)
  {
    const unique_ptr<Neighbours> &neighbours = kv.second;
    const Event &refEv = catalog.getEvents().at(neighbours->refEvId);

    // Compute theoretical phases for stations that have no picks. The
    // cross-correlation will be used to detect and fix the pick time.
    if (computeTheoreticalPhases)
    {
      addMissingEventPhases(refEv, catalog, catalog, *neighbours);
    }

    buildXcorrDiffTTimePairs(catalog, *neighbours, refEv, xcorrMaxEvStaDist,
                             xcorrMaxInterEvDist, xcorr);

    // Update theoretical and automatic phase pick time and uncertainties
    // based on cross-correlation results. Also, drop theoretical phases
    // wihout any good cross-correlation result.
    fixPhases(catalog, refEv, xcorr);

    const unsigned onePerThousand =
        neighCluster.size() < 1000 ? 1 : (neighCluster.size() / 1000);
    if (++performed % onePerThousand == 0)
    {
      logInfo("Cross-correlation completion %.1f%%",
              (performed * 100 / (double)neighCluster.size()));
    }
  }

  WfCounters wfcount{};
  wfcount.update(_wfAccess.loader.get());
  wfcount.update(_wfAccess.diskCache.get());

  logInfo("Catalog waveform data: waveforms downloaded %u, not available "
          "%u, loaded from disk cache %u",
          wfcount.downloaded, wfcount.no_avail, wfcount.disk_cached);

  logXCorrSummary(xcorr);

  return xcorr;
}

/*
 * Compute and store to `XCorrCache` cross-correlated differential travel
 * times for pairs of the earthquake.
 */
void DD::buildXcorrDiffTTimePairs(Catalog &catalog,
                                  const Neighbours &neighbours,
                                  const Event &refEv,
                                  double xcorrMaxEvStaDist,
                                  double xcorrMaxInterEvDist,
                                  XCorrCache &xcorr)
{
  logDebug("Computing cross-correlation differential travel times for event %s",
           string(refEv).c_str());
  //
  // Prepare the waveform loaders for single-event. Since they are
  // temporary and discarded after the relocation we don't want to
  // make use of the background catalog caching and keep them in
  // memory forever
  //
  shared_ptr<Waveform::Processor> preloaded(preloadNonCatalogWaveforms(
      catalog, neighbours, refEv, xcorrMaxEvStaDist, xcorrMaxInterEvDist));
  shared_ptr<Waveform::Processor> seWfLdrNoSnr; // loader with not SNR check
  shared_ptr<Waveform::Processor> seWfLdr; // same loader but with SNR check

  if (preloaded)
  {
    seWfLdrNoSnr.reset(new Waveform::MemCachedProc(preloaded));
    seWfLdr = seWfLdrNoSnr;
    if (_cfg.snr.minSnr > 0)
    {
      shared_ptr<Waveform::SnrFilterPrc> snrLdr(new Waveform::SnrFilterPrc(
          seWfLdrNoSnr, _cfg.snr.minSnr, _cfg.snr.noiseStart, _cfg.snr.noiseEnd,
          _cfg.snr.signalStart, _cfg.snr.signalEnd));
      seWfLdr.reset(new Waveform::MemCachedProc(snrLdr));
    }
  }

  //
  // loop through reference event phases
  //
  auto eqlrngRef = catalog.getPhases().equal_range(refEv.id);
  for (auto itRef = eqlrngRef.first; itRef != eqlrngRef.second; ++itRef)
  {
    const Phase &refPhase  = itRef->second;
    const Station &station = catalog.getStations().at(refPhase.stationId);
    const auto xcorrCfg    = _cfg.xcorr.at(refPhase.procInfo.type);

    //
    // skip stations too far away
    //
    double stationDistance = computeDistance(refEv, station);
    if (stationDistance > xcorrMaxEvStaDist && xcorrMaxEvStaDist >= 0) continue;

    //
    // select appropriate waveform loader for refPhase
    //
    shared_ptr<Waveform::Processor> refPrc;
    if (refPhase.procInfo.source == Phase::Source::CATALOG)
      refPrc = _wfAccess.memCache;
    else if (refPhase.isManual)
      refPrc = seWfLdr;
    else // not from catalog, not manual
    {
      // For "untrusted" phases the SNR is checked AFTER the cross-correlation
      // once the pick time has been adjusted using the lag computed by the
      // cross-corr against a "trusted" phase (manual or catalog).
      refPrc = seWfLdrNoSnr;
    }

    //
    // loop through neighbouring events and cross-correlate with `refPhase`
    //
    for (unsigned neighEvId : neighbours.ids)
    {
      const Event &event = catalog.getEvents().at(neighEvId);

      //
      // skip events too far away
      //
      double interEventDistance = computeDistance(refEv, event);
      if (interEventDistance > xcorrMaxInterEvDist && xcorrMaxInterEvDist >= 0)
        continue;

      if (neighbours.has(neighEvId, refPhase.stationId, refPhase.procInfo.type))
      {
        const Phase &phase = catalog
                                 .searchPhase(event.id, refPhase.stationId,
                                              refPhase.procInfo.type)
                                 ->second;

        // In single-event mode `refPhase` is real-time and `phase` is from
        // the catalog. In multi-event mode both are from the catalog.
        if (phase.procInfo.source != Phase::Source::CATALOG)
          throw Exception("Internal logic error: phase is not from catalog");

        // skipp cross-correlation if we already have the result in cache
        if (!xcorr.has(refEv.id, event.id, refPhase.stationId,
                       refPhase.procInfo.type))
        {
          // Do the actual cross-correlation
          double coeff, lag;
          string component;
          bool valid = xcorrPhases(refEv, refPhase, *refPrc, event, phase,
                                   *_wfAccess.memCache, coeff, lag, component);

          // Check the SNR (using the pick time adjusted with the
          // cross-correlation lag) of those phases, where the check hadn't
          // been performed yet.
          if (valid && coeff >= xcorrCfg.minCoef && _cfg.snr.minSnr > 0 &&
              !refPhase.isManual &&
              refPhase.procInfo.source != Phase::Source::CATALOG)
          {
            UTCTime adjustedPickTime = refPhase.time - secToDur(lag);
            TimeWindow snrWin =
                _wfAccess.snrFilter->snrTimeWindow(adjustedPickTime);

            shared_ptr<const Trace> trace =
                getWaveform(*seWfLdrNoSnr, snrWin, refEv, refPhase, component);

            if (!trace ||
                !_wfAccess.snrFilter->goodSnr(*trace, adjustedPickTime))
            {
              valid = false;
            }
          }

          // cache cross-correlation results
          xcorr.add(refEv.id, event.id, refPhase.stationId,
                    refPhase.procInfo.type, valid, coeff, lag, component);
        }
      }
    }
  }
}

shared_ptr<Waveform::Processor>
DD::preloadNonCatalogWaveforms(Catalog &catalog,
                               const Neighbours &neighbours,
                               const Event &refEv,
                               double xcorrMaxEvStaDist,
                               double xcorrMaxInterEvDist)
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

  // Load a large enough waveform to be able to check the SNR after
  // the updating of the pick time
  if (_cfg.snr.minSnr > 0)
  {
    double maxDelay = std::max({_cfg.xcorr.at(Phase::Type::P).maxDelay,
                                _cfg.xcorr.at(Phase::Type::S).maxDelay});
    double extraBefore =
        std::min({_cfg.snr.noiseStart, _cfg.snr.signalStart}) - maxDelay;
    extraBefore = extraBefore > 0 ? 0 : std::abs(extraBefore);
    extraBefore += _cfg.wfFilter.extraTraceLen;
    double extraAfter =
        std::max({_cfg.snr.noiseEnd, _cfg.snr.signalEnd}) + maxDelay;
    extraAfter = extraAfter < 0 ? 0 : extraAfter;
    extraAfter += _cfg.wfFilter.extraTraceLen;
    loader.reset(new Waveform::ExtraLenLoader(loader, extraBefore, extraAfter));
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
    const auto components  = xcorrComponents(refPhase);

    // We deal only with real-time event data
    if (refPhase.procInfo.source == Phase::Source::CATALOG) continue;
    onlyCatalogPhases = false;

    //
    // skip stations too far away
    //
    double stationDistance = computeDistance(refEv, station);
    if (stationDistance > xcorrMaxEvStaDist && xcorrMaxEvStaDist >= 0) continue;

    //
    // loop through neighbouring events and cross-correlate with `refPhase`
    //
    for (unsigned neighEvId : neighbours.ids)
    {
      const Event &event = catalog.getEvents().at(neighEvId);

      //
      // skip events too far away
      //
      double interEventDistance = computeDistance(refEv, event);
      if (interEventDistance > xcorrMaxInterEvDist && xcorrMaxInterEvDist >= 0)
        continue;

      if (neighbours.has(neighEvId, refPhase.stationId, refPhase.procInfo.type))
      {
        //
        // For each match load the reference event phase waveforms
        //
        TimeWindow tw = xcorrTimeWindowLong(refPhase);

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
  logInfo("Event %s: waveforms downloaded %u, not available %u, loaded from "
          "disk cache %u",
          string(refEv).c_str(), wfcount.downloaded, wfcount.no_avail,
          wfcount.disk_cached);
  return proc;
}

/*
 * Update theoretical and automatic phase pick times and uncertainties based
 * on cross-correlation results. Drop theoretical phases not passing the
 * cross-correlation verification.
 */
void DD::fixPhases(Catalog &catalog,
                   const Event &refEv,
                   XCorrCache &xcorr) const
{
  unsigned totP = 0, totS = 0;
  unsigned newP = 0, newS = 0;

  vector<Phase> phasesToBeRemoved;
  vector<Phase> newPhases;

  // loop through catalog phases
  auto eqlrng = catalog.getPhases().equal_range(refEv.id);
  for (auto it = eqlrng.first; it != eqlrng.second; it++)
  {
    const Phase &phase  = it->second;
    const auto xcorrCfg = _cfg.xcorr.at(phase.procInfo.type);

    if (phase.procInfo.type == Phase::Type::P) totP++;
    if (phase.procInfo.type == Phase::Type::S) totS++;

    // nothing to do if the phase is manual or from catalog
    if (phase.isManual || phase.procInfo.source == Phase::Source::CATALOG)
    {
      continue;
    }

    //
    // Load all good xcorr result for this automatic phase
    //
    map<double, std::pair<std::reference_wrapper<const Phase>,
                          std::reference_wrapper<const XCorrCache::Entry>>>
        goodXcorr; // this is sorted
    if (xcorr.has(refEv.id, phase.stationId, phase.procInfo.type))
    {
      for (const auto &kv :
           xcorr.get(refEv.id, phase.stationId, phase.procInfo.type))
      {
        unsigned neighbourId       = kv.first;
        const XCorrCache::Entry &e = kv.second;
        if (e.valid && e.coeff >= xcorrCfg.minCoef)
        {
          const Phase &neighbourPhase =
              _bgCat
                  .searchPhase(neighbourId, phase.stationId,
                               phase.procInfo.type)
                  ->second;

          goodXcorr.insert({e.lag, {neighbourPhase, e}});
        }
      }
    }

    //
    // Nothing to do if we dont't have at least 3 good xcorr results
    // (3 is kind of arbitrary, I know)
    //
    if (goodXcorr.size() < 3)
    {
      // remove thoretical phases without enough good cross-correlations
      if (phase.procInfo.source == Phase::Source::THEORETICAL)
        phasesToBeRemoved.push_back(phase);
      continue;
    }

    //
    // Use median-1, median, median+1 values to update phase time and
    // uncertainty
    //
    double mean_lag  = 0;
    double min_lag   = 0;
    double max_lag   = 0;
    double totWeight = 0;
    string component;

    size_t median = goodXcorr.size() / 2;
    auto kv       = goodXcorr.begin();
    std::advance(kv, median - 1);
    for (int i = 0; i < 3; i++, kv++)
    {
      const Phase &neighbourPhase = kv->second.first;
      const XCorrCache::Entry &e  = kv->second.second;
      double weight               = e.coeff * neighbourPhase.procInfo.weight;
      mean_lag += e.lag * weight;
      min_lag += (e.lag - neighbourPhase.lowerUncertainty) * weight;
      max_lag += (e.lag + neighbourPhase.upperUncertainty) * weight;
      totWeight += weight;
      component = e.component; // this is a little ugly
    }
    mean_lag /= totWeight;
    min_lag /= totWeight;
    max_lag /= totWeight;

    //
    // set new phase time and uncertainty
    //
    Phase newPhase(phase);
    newPhase.time -= secToDur(mean_lag);
    newPhase.lowerUncertainty = mean_lag - min_lag;
    newPhase.upperUncertainty = max_lag - mean_lag;
    newPhase.procInfo.weight  = Catalog::computePickWeight(newPhase);
    newPhase.procInfo.source  = Phase::Source::XCORR;
    newPhase.type = strf("%cx", static_cast<char>(newPhase.procInfo.type));

    if (phase.procInfo.source == Phase::Source::THEORETICAL)
    {
      newPhase.channelCode =
          getBandAndInstrumentCodes(newPhase.channelCode) + component;

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
    catalog.removePhase(ph.eventId, ph.stationId, ph.procInfo.type);
    xcorr.remove(ph.eventId, ph.stationId, ph.procInfo.type);
  }

  for (Phase &ph : newPhases)
  {
    catalog.updatePhase(ph, true);
  }

  logDebug("Event %s total phases %u (%u P and %u S): created %u (%u P "
           "and %u S) from theoretical picks",
           string(refEv).c_str(), (totP + totS), totP, totS, (newP + newS),
           newP, newS);
}

TimeWindow DD::xcorrTimeWindowLong(const Phase &phase) const
{
  const auto xcorrCfg = _cfg.xcorr.at(phase.procInfo.type);
  TimeWindow tw       = xcorrTimeWindowShort(phase);
  tw.setStartTime(tw.startTime() - secToDur(xcorrCfg.maxDelay));
  tw.setEndTime(tw.endTime() + secToDur(xcorrCfg.maxDelay));
  return tw;
}

TimeWindow DD::xcorrTimeWindowShort(const Phase &phase) const
{
  const auto xcorrCfg = _cfg.xcorr.at(phase.procInfo.type);
  UTCTime::duration shortDuration =
      secToDur(xcorrCfg.endOffset - xcorrCfg.startOffset);
  UTCTime::duration shortTimeCorrection = secToDur(xcorrCfg.startOffset);
  return TimeWindow(phase.time + shortTimeCorrection, shortDuration);
}

bool DD::xcorrPhases(const Event &event1,
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
    logError("Skipping cross-correlation: mismatching phases (%s and %s)",
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
      logDebug("Skipping cross-correlation: incompatible channels %s and %s "
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
  for (const string &component : xcorrComponents(phase1))
  {
    performed =
        xcorrPhasesOneComponent(event1, phase1, ph1Cache, event2, phase2,
                                ph2Cache, component, coeffOut, lagOut);
    if (performed)
    {
      componentOut = component;
      break;
    }
  }
  return performed;
}

bool DD::xcorrPhasesOneComponent(const Event &event1,
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

  auto xcorrCfg = _cfg.xcorr.at(phase1.procInfo.type);

  TimeWindow tw1 = xcorrTimeWindowLong(phase1);
  TimeWindow tw2 = xcorrTimeWindowLong(phase2);

  // Load the long `tr1`, because we want to cache the long version. Then
  // we'll trim it.
  shared_ptr<const Trace> tr1 =
      getWaveform(ph1Cache, tw1, event1, phase1, component);
  if (!tr1)
  {
    logDebug(
        "Skipping cross-correlation: no waveform or no good SNR for phase %s",
        string(phase1).c_str());
    return false;
  }

  // Load the long `tr2`, because we want to cache the long version. Then
  // we'll trim it.
  shared_ptr<const Trace> tr2 =
      getWaveform(ph2Cache, tw2, event2, phase2, component);
  if (!tr2)
  {
    logDebug(
        "Skipping cross-correlation: no waveform or no good SNR for phase %s",
        string(phase2).c_str());
    return false;
  }

  if (tr1->samplingFrequency() != tr2->samplingFrequency())
  {
    logWarning(
        "Skipping cross-correlation: traces have different sampling freq "
        "(%f!=%f): phase1 %s and phase 2 %s",
        tr1->samplingFrequency(), tr2->samplingFrequency(),
        string(phase1).c_str(), string(phase2).c_str());
    return false;
  }

  // Trust the manual pick on `phase2`: keep `tr2` short and cross-correlate
  // it with the larger `tr1` window.
  double xcorr_coeff = numeric_limits<double>::quiet_NaN(), xcorr_lag = 0;

  if (phase2.isManual || (!phase1.isManual && !phase2.isManual))
  {
    // Trim `tr2` to shorter length; we want to cross-correlate the short one
    // with the long one.
    Trace tr2Short(*tr2);
    TimeWindow tw2Short = xcorrTimeWindowShort(phase2);
    if (!tr2Short.slice(tw2Short))
    {
      logDebug("Skipping cross-correlation: cannot trim phase2 waveform "
               "for phase pair phase1='%s', phase2='%s'",
               string(phase1).c_str(), string(phase2).c_str());
      return false;
    }

    xcorr(*tr1, tr2Short, xcorrCfg.maxDelay, xcorr_lag, xcorr_coeff);
  }

  // Trust the manual pick on `phase1`: keep `tr1` short and cross-correlate
  // it with a larger `tr2` window.
  double xcorr_coeff2 = numeric_limits<double>::quiet_NaN(), xcorr_lag2 = 0;

  if (phase1.isManual || (!phase1.isManual && !phase2.isManual))
  {
    // Trim `tr1` to shorter length; we want to cross-correlate the short with
    // the long one.
    Trace tr1Short(*tr1);
    TimeWindow tw1Short = xcorrTimeWindowShort(phase1);
    if (!tr1Short.slice(tw1Short))
    {
      logDebug("Skipping cross-correlation: cannot trim phase1 waveform "
               "for phase pair phase1='%s', phase2='%s'",
               string(phase1).c_str(), string(phase2).c_str());
      return false;
    }

    xcorr(tr1Short, *tr2, xcorrCfg.maxDelay, xcorr_lag2, xcorr_coeff2);
  }

  if (!std::isfinite(xcorr_coeff) && !std::isfinite(xcorr_coeff2))
  {
    logDebug("Skipping cross-correlation: no meaningful coefficient for phase "
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
 * respective picks. The cross-correlation will be performed with the short
 * trace against the longest trace starting from the longest trace middle
 * point minus 'maxDelay' until middle point pluse 'maxDelay'. `delayOut` will
 * store the shift in seconds (positive or negative) from the longest trace
 * middle point at which there is the highest (absolute value) correlation
 * coefficient, stored in 'coeffOut'
 */
void DD::xcorr(const Trace &tr1,
               const Trace &tr2,
               double maxDelay,
               double &delayOut,
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
                   (sizeS + maxDelaySmps * 2), delayOut, coeffOut);

  delayOut -= maxDelaySmps; // the reference is the middle of the long trace
  delayOut /= freq;         // samples to secs
  if (swap) delayOut = -delayOut;
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

void DD::logXCorrSummary(const XCorrCache &xcorr)
{
  struct
  {
    unsigned skipped;
    unsigned performed;
    unsigned performed_s;
    unsigned performed_p;
    unsigned good_cc;
    unsigned good_cc_s;
    unsigned good_cc_p;
  } counters{};

  auto callback = [&counters, this](unsigned ev1, unsigned ev2,
                                    const std::string &stationId,
                                    const Catalog::Phase::Type &type,
                                    const XCorrCache::Entry &e) {
    if (e.valid)
    {
      const auto xcorrCfg = _cfg.xcorr.at(type);
      bool goodCoeff      = e.coeff >= xcorrCfg.minCoef;

      counters.performed++;
      if (type == Phase::Type::S)
        counters.performed_s++;
      else if (type == Phase::Type::P)
        counters.performed_p++;

      if (goodCoeff)
      {
        counters.good_cc++;
        if (type == Phase::Type::S)
          counters.good_cc_s++;
        else if (type == Phase::Type::P)
          counters.good_cc_p++;
      }
    }
    else
      counters.skipped++;
  };

  xcorr.forEach(callback);

  // each xcorr is reported twice in XCorrCache (ev1-ev2 and ev2-ev1)
  counters.skipped /= 2;
  counters.performed /= 2;
  counters.performed_s /= 2;
  counters.performed_p /= 2;
  counters.good_cc /= 2;
  counters.good_cc_s /= 2;
  counters.good_cc_p /= 2;

  logInfo("Cross-correlation performed %u (P phase %.f%%, S phase %.f%%), "
          "skipped %u (%.f%%)",
          counters.performed,
          (counters.performed_p * 100. / counters.performed),
          (counters.performed_s * 100. / counters.performed), counters.skipped,
          (counters.skipped * 100. / (counters.performed + counters.skipped)));

  logInfo("Successful cross-correlation (coefficient above threshold) %.1f%% "
          "(%u/%u). Successful P %.1f%% (%u/%u). Successful S %.1f%% (%u/%u)",
          (counters.good_cc * 100.0 / counters.performed), counters.good_cc,
          counters.performed,
          (counters.good_cc_p * 100.0 / counters.performed_p),
          counters.good_cc_p, counters.performed_p,
          (counters.good_cc_s * 100.0 / counters.performed_s),
          counters.good_cc_s, counters.performed_s);
}

namespace {

struct CallbackCounters
{
  unsigned skipped   = 0;
  unsigned performed = 0;
  vector<double> coeff;
  vector<double> lag;
};

HDD::DD::XCorrEvalStats &operator+=(HDD::DD::XCorrEvalStats &lhs,
                                    CallbackCounters const &rhs)
{
  lhs.skipped.push_back(rhs.skipped);
  lhs.performed.push_back(rhs.performed);
  lhs.coeff.insert(lhs.coeff.end(), rhs.coeff.begin(), rhs.coeff.end());
  lhs.lag.insert(lhs.lag.end(), rhs.lag.begin(), rhs.lag.end());
  return lhs;
}

} // namespace

void DD::XCorrEvalStats::summarize(unsigned &skipped,
                                   unsigned &performed,
                                   double &meanCoeff,
                                   double &meanCoeffAbsDev,
                                   double &medianCoeff,
                                   double &medianCoeffAbsDev,
                                   double &meanLag,
                                   double &meanLagAbsDev,
                                   double &medianLag,
                                   double &medianLagAbsDev) const
{
  skipped = std::accumulate(this->skipped.begin(), this->skipped.end(), 0.);
  performed =
      std::accumulate(this->performed.begin(), this->performed.end(), 0.);
  meanCoeff         = computeMean(this->coeff);
  meanCoeffAbsDev   = computeMeanAbsoluteDeviation(this->coeff, meanCoeff);
  medianCoeff       = computeMedian(this->coeff);
  medianCoeffAbsDev = computeMedianAbsoluteDeviation(this->coeff, medianCoeff);
  meanLag           = computeMean(this->lag);
  meanLagAbsDev     = computeMeanAbsoluteDeviation(this->lag, meanLag);
  medianLag         = computeMedian(this->lag);
  medianLagAbsDev   = computeMedianAbsoluteDeviation(this->lag, medianLag);
}

void DD::evalXCorr(const ClusteringOptions &clustOpt,
                   const evalXcorrCallback &cb,
                   XCorrCache &xcorr,
                   const double interEvDistStep,
                   const double staDistStep,
                   bool theoretical)
{
  XCorrEvalStats pTotStats, sTotStats;
  map<string, XCorrEvalStats> pStatsByStation;           // key station id
  map<string, XCorrEvalStats> sStatsByStation;           // key station id
  map<unsigned, XCorrEvalStats> pStatsByStaDistance;     // key distance
  map<unsigned, XCorrEvalStats> sStatsByStaDistance;     // key distance
  map<unsigned, XCorrEvalStats> pStatsByInterEvDistance; // key distance
  map<unsigned, XCorrEvalStats> sStatsByInterEvDistance; // key distance

  int processedEvents = 0;

  for (const auto &kv : _bgCat.getEvents())
  {
    const Event &event = kv.second;

    // find the neighbouring events
    unique_ptr<Neighbours> neighbours;
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

    unique_ptr<Catalog> catalog;

    if (theoretical)
    {
      // create theoretical phases for this event instead of fetching its
      // phases from the catalog
      catalog = neighbours->toCatalog(_bgCat, false);
      addMissingEventPhases(event, *catalog, _bgCat, *neighbours);
    }
    else
    {
      catalog = neighbours->toCatalog(_bgCat, true);
    }

    // Cross-correlate every neighbour phase with its corresponding event
    // theoretical phase.
    buildXcorrDiffTTimePairs(*catalog, *neighbours, event,
                             clustOpt.xcorrMaxEvStaDist,
                             clustOpt.xcorrMaxInterEvDist, xcorr);

    // Update theoretical and automatic phase pick time and uncertainties
    // based on cross-correlation results. Drop theoretical phases wihout any
    // good cross-correlation result.
    if (theoretical)
    {
      fixPhases(*catalog, event, xcorr);
    }

    //
    // Collect Stats From the XcorrCache
    //
    for (const auto &kv : neighbours->allPhases())
    {
      for (Phase::Type phaseType : kv.second)
      {
        const string stationId = kv.first;
        const auto xcorrCfg    = _cfg.xcorr.at(phaseType);

        //
        // Collect stats for current station/phaseType
        //
        CallbackCounters counters;
        auto callback = [&counters, &neighbours,
                         &xcorrCfg](unsigned ev1, unsigned ev2,
                                    const std::string &stationId,
                                    const Catalog::Phase::Type &type,
                                    const XCorrCache::Entry &e) {
          if (neighbours->has(ev2, stationId, type))
          {
            if (e.valid && (e.coeff >= xcorrCfg.minCoef))
            {
              counters.performed++;
              counters.coeff.push_back(e.coeff);
              counters.lag.push_back(e.lag);
            }
            else
              counters.skipped++;
          }
        };
        xcorr.forEach(event.id, stationId, phaseType, callback);

        //
        //  group stats by event, station, station distance
        //
        const Station &station = _bgCat.getStations().at(stationId);
        double stationDistance = computeDistance(event, station);
        unsigned stationDistanceBucket =
            unsigned(stationDistance / staDistStep);

        if (phaseType == Phase::Type::P)
        {
          pTotStats += counters;
          pStatsByStation[stationId] += counters;
          pStatsByStaDistance[stationDistanceBucket] += counters;
        }
        else if (phaseType == Phase::Type::S)
        {
          sTotStats += counters;
          sStatsByStation[stationId] += counters;
          sStatsByStaDistance[stationDistanceBucket] += counters;
        }

        //
        // group stats by inter-event distance
        //
        for (unsigned neighEvId : neighbours->ids)
        {
          if (!neighbours->has(neighEvId, stationId, phaseType)) continue;
          if (!xcorr.has(event.id, neighEvId, stationId, phaseType)) continue;

          const auto &e = xcorr.get(event.id, neighEvId, stationId, phaseType);
          CallbackCounters counters;
          if (e.valid && (e.coeff >= xcorrCfg.minCoef))
          {
            counters.performed++;
            counters.coeff.push_back(e.coeff);
            counters.lag.push_back(e.lag);
          }
          else
            counters.skipped++;

          const Event &neighbEv  = catalog->getEvents().at(neighEvId);
          double interEvDistance = computeDistance(event, neighbEv);
          const unsigned interEvDistanceBucket =
              unsigned(interEvDistance / interEvDistStep);

          if (phaseType == Phase::Type::P)
            pStatsByInterEvDistance[interEvDistanceBucket] += counters;
          else if (phaseType == Phase::Type::S)
            sStatsByInterEvDistance[interEvDistanceBucket] += counters;
        }
      }
    }

    //
    // User update
    //
    processedEvents++;
    const unsigned onePercent =
        _bgCat.getEvents().size() < 100 ? 1 : (_bgCat.getEvents().size() / 100);
    if (processedEvents % onePercent == 0)
    {
      logInfo("Processed %.1f%%...",
              (processedEvents * 100.0 / _bgCat.getEvents().size()));
    }

    if (processedEvents % (10 * onePercent) == 0)
    {
      cb(pTotStats, sTotStats, pStatsByStation, sStatsByStation,
         pStatsByStaDistance, sStatsByStaDistance, pStatsByInterEvDistance,
         sStatsByInterEvDistance, interEvDistStep, staDistStep,
         (double(processedEvents) / _bgCat.getEvents().size()));
    }
  }

  cb(pTotStats, sTotStats, pStatsByStation, sStatsByStation,
     pStatsByStaDistance, sStatsByStaDistance, pStatsByInterEvDistance,
     sStatsByInterEvDistance, interEvDistStep, staDistStep, 1.0);

  WfCounters wfcount{};
  wfcount.update(_wfAccess.loader.get());
  wfcount.update(_wfAccess.diskCache.get());
  logInfo("Catalog waveforms downloaded %u, not available "
          "%u, loaded from disk cache %u",
          wfcount.downloaded, wfcount.no_avail, wfcount.disk_cached);
}

} // namespace HDD

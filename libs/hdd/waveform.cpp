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

#include "waveform.h"
#include "log.h"
#include "timewindow.h"
#include "utils.h"

using namespace std;
using Catalog         = HDD::Catalog;
using ThreeComponents = HDD::Waveform::ThreeComponents;
using Transform       = HDD::Waveform::Processor::Transform;

namespace {

string waveformId(const HDD::TimeWindow &tw,
                  const string &networkCode,
                  const string &stationCode,
                  const string &locationCode,
                  const string &channelCode)
{
  return HDD::strf("%s.%s.%s.%s.%s.%s", networkCode.c_str(),
                   stationCode.c_str(), locationCode.c_str(),
                   channelCode.c_str(),
                   HDD::UTCClock::toString(tw.startTime()).c_str(),
                   HDD::UTCClock::toString(tw.endTime()).c_str());
}

std::string waveformId(const HDD::TimeWindow &tw, const Catalog::Phase &ph)
{
  return waveformId(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                    ph.channelCode);
}

} // namespace

namespace HDD {
namespace Waveform {

shared_ptr<const Trace> BasicLoader::get(const TimeWindow &tw,
                                         const Catalog::Phase &ph)
{
  unique_ptr<Trace> trace;

  try
  {
    trace = _wf->loadTrace(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                           ph.channelCode);
  }
  catch (exception &e)
  {
    logDebug(
        "Cannnot load trace (stream %s.%s.%s.%s from %s length %.2f sec): %s",
        ph.networkCode.c_str(), ph.stationCode.c_str(), ph.locationCode.c_str(),
        ph.channelCode.c_str(), UTCClock::toString(tw.startTime()).c_str(),
        durToSec(tw.length()), e.what());
  }

  if (!trace)
    _counters_wf_no_avail++;
  else
    _counters_wf_downloaded++;

  return trace;
}

shared_ptr<const Trace> BatchLoader::get(const TimeWindow &tw,
                                         const Catalog::Phase &ph)
{
  if (!_dataLoaded)
  {
    request(tw, ph);
    return nullptr;
  }
  else
  {
    const string wfId = waveformId(tw, ph);
    const auto it     = _traces.find(wfId);
    return it != _traces.end() ? it->second : nullptr;
  }
}

void BatchLoader::request(const TimeWindow &tw, const Catalog::Phase &ph)
{
  if (_dataLoaded)
  {
    throw Exception("Cannot request more traces after they have been loaded");
  }

  const string streamID = ph.networkCode + "." + ph.stationCode + "." +
                          ph.locationCode + "." + ph.channelCode;
  auto eqlrng = _requests.equal_range(streamID);
  for (auto it = eqlrng.first; it != eqlrng.second; ++it)
  {
    if (it->second == tw) return; // skip duplicated requests
  }
  _requests.emplace(streamID, tw);
}

void BatchLoader::load()
{
  if (_dataLoaded) return;

  const auto onTraceLoaded = [this](const std::string &streamID,
                                    const TimeWindow &tw,
                                    unique_ptr<Trace> trace) {
    const string wfId =
        waveformId(tw, trace->networkCode(), trace->stationCode(),
                   trace->locationCode(), trace->channelCode());
    _traces[wfId] = std::move(trace);
    _counters_wf_downloaded++;
  };

  const auto onTraceFailed = [this](const std::string &streamID,
                                    const TimeWindow &tw,
                                    const std::string &error) {
    _counters_wf_no_avail++;
    logDebug("Cannnot load trace (stream %s from %s length %.2f sec): %s",
             streamID.c_str(), UTCClock::toString(tw.startTime()).c_str(),
             durToSec(tw.length()), error.c_str());
  };

  if (_requests.size() > 0)
  {
    _wf->loadTraces(_requests, onTraceLoaded, onTraceFailed);

    logInfo("Fetched %u/%lu waveforms, not available %u",
            _counters_wf_downloaded, _requests.size(), _counters_wf_no_avail);

    _requests.clear();
  }
  _dataLoaded = true;
}

TimeWindow ExtraLenLoader::traceTimeWindowToLoad(const TimeWindow &neededTW,
                                                 const UTCTime &pickTime) const
{
  TimeWindow extraTW(pickTime - secToDur(_beforePickLen),
                     pickTime + secToDur(_afterPickLen));
  TimeWindow twToLoad = extraTW.merge(neededTW);

  //
  // round the start/end time to the nearest second
  //
  int year, month, day, hour, min, sec, usec;

  UTCClock::toDate(twToLoad.startTime(), year, month, day, hour, min, sec,
                   usec);
  if (usec > 0)
  {
    twToLoad.setStartTime(
        UTCClock::fromDate(year, month, day, hour, min, sec, 0));
  }

  UTCClock::toDate(twToLoad.endTime(), year, month, day, hour, min, sec, usec);
  if (usec > 0)
  {
    twToLoad.setEndTime(
        UTCClock::fromDate(year, month, day, hour, min, sec + 1, 0));
  }

  return twToLoad;
}

shared_ptr<const Trace> ExtraLenLoader::get(const TimeWindow &tw,
                                            const Catalog::Phase &ph)
{
  const TimeWindow twToLoad     = traceTimeWindowToLoad(tw, ph.time);
  shared_ptr<const Trace> trace = _auxLdr->get(twToLoad, ph);
  if (trace && twToLoad != tw)
  {
    shared_ptr<Trace> nonConstTrace(new Trace(*trace));
    if (!nonConstTrace->slice(tw))
    {
      logDebug("Error while loading phase '%s': cannot slice trace "
               "from %s length %.2f sec. Trace "
               "data from %s length %.2f sec, samples %zu sampfreq %f",
               string(ph).c_str(),
               HDD::UTCClock::toString(tw.startTime()).c_str(),
               HDD::durToSec(tw.length()),
               HDD::UTCClock::toString(trace->startTime()).c_str(),
               HDD::durToSec(trace->timeWindow().length()),
               trace->sampleCount(), trace->samplingFrequency());
      return nullptr;
    }
    trace = nonConstTrace;
  }
  return trace;
}

shared_ptr<const Trace> DiskCachedLoader::get(const TimeWindow &tw,
                                              const Catalog::Phase &ph)
{
  bool isCached;
  shared_ptr<const Trace> trace = getFromCache(
      tw, ph.networkCode, ph.stationCode, ph.locationCode, ph.channelCode);
  if (trace)
  {
    isCached = true;
    _counters_wf_cached++;
  }
  else
  {
    isCached = false;
    trace    = _auxLdr->get(tw, ph);
  }

  if (trace && !isCached)
  {
    storeInCache(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                 ph.channelCode, *trace);
  }
  return trace;
}

bool DiskCachedLoader::isCached(const TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev) const
{
  const string cacheFile =
      waveformPath(_cacheDir, tw, ph.networkCode, ph.stationCode,
                   ph.locationCode, ph.channelCode);
  return pathExists(cacheFile);
}

unique_ptr<Trace> DiskCachedLoader::getFromCache(const TimeWindow &tw,
                                                 const string &networkCode,
                                                 const string &stationCode,
                                                 const string &locationCode,
                                                 const string &channelCode)
{
  const string cacheFile = waveformPath(_cacheDir, tw, networkCode, stationCode,
                                        locationCode, channelCode);

  if (!pathExists(cacheFile)) return nullptr;

  try
  {
    return _wf->readTrace(cacheFile);
  }
  catch (exception &e)
  {
    logWarning("Couldn't load waveform %s: %s", cacheFile.c_str(), e.what());
    return nullptr;
  }
}

void DiskCachedLoader::storeInCache(const TimeWindow &tw,
                                    const string &networkCode,
                                    const string &stationCode,
                                    const string &locationCode,
                                    const string &channelCode,
                                    const Trace &trace)
{
  const string cacheFile = waveformPath(_cacheDir, tw, networkCode, stationCode,
                                        locationCode, channelCode);
  try
  {
    _wf->writeTrace(trace, cacheFile);
  }
  catch (exception &e)
  {
    logWarning("Couldn't write waveform to disk %s: %s", cacheFile.c_str(),
               e.what());
  }
}

string DiskCachedLoader::waveformPath(const string &cacheDir,
                                      const TimeWindow &tw,
                                      const string &networkCode,
                                      const string &stationCode,
                                      const string &locationCode,
                                      const string &channelCode) const
{
  string cacheFile =
      waveformId(tw, networkCode, stationCode, locationCode, channelCode) +
      ".mseed";
  return joinPath(cacheDir, cacheFile);
}

shared_ptr<const Trace> BasicProcessor::get(const TimeWindow &tw,
                                            const Catalog::Phase &ph,
                                            const Catalog::Event &ev,
                                            const Catalog::Station &sta,
                                            const std::string &filterStr,
                                            double resampleFreq,
                                            Transform trans)
{
  if (trans == Transform::NONE)
  {
    return loadAndProcess(tw, ph, ph.channelCode, filterStr, resampleFreq);
  }

  //
  // Apply the transformation
  //
  ThreeComponents comps;
  try
  {
    _wf->getComponentsInfo(ph, comps);
  }
  catch (exception &e)
  {
    logDebug("Cannot compute RT/L2 transformation of stream %s.%s.%s.%s "
             "from %s length %.2f sec: %s",
             ph.networkCode.c_str(), ph.stationCode.c_str(),
             ph.locationCode.c_str(),
             getBandAndInstrumentCodes(ph.channelCode).c_str(),
             UTCClock::toString(tw.startTime()).c_str(), durToSec(tw.length()),
             e.what());
    return nullptr;
  }

  shared_ptr<const Trace> trace;

  if (trans == Transform::TRANSVERSAL || trans == Transform::RADIAL)
  {
    shared_ptr<const Trace> trV(
        loadAndProcess(tw, ph, comps.names[ThreeComponents::Vertical],
                       filterStr, resampleFreq));
    shared_ptr<const Trace> trH1(
        loadAndProcess(tw, ph, comps.names[ThreeComponents::FirstHorizontal],
                       filterStr, resampleFreq));
    shared_ptr<const Trace> trH2(
        loadAndProcess(tw, ph, comps.names[ThreeComponents::SecondHorizontal],
                       filterStr, resampleFreq));
    if (trH1 && trH2 && trV)
      trace = transformRT(tw, ph, ev, sta, comps, *trV, *trH1, *trH2, trans);
  }
  else if (trans == Transform::L2)
  {
    shared_ptr<const Trace> trH1(
        loadAndProcess(tw, ph, comps.names[ThreeComponents::FirstHorizontal],
                       filterStr, resampleFreq));
    shared_ptr<const Trace> trH2(
        loadAndProcess(tw, ph, comps.names[ThreeComponents::SecondHorizontal],
                       filterStr, resampleFreq));
    if (trH1 && trH2) trace = transformL2(tw, ph, comps, *trH1, *trH2);
  }

  if (!trace)
  {
    logDebug("Cannot compute RT/L2 transformation of stream %s.%s.%s.%s "
             "from %s length %.2f sec",
             ph.networkCode.c_str(), ph.stationCode.c_str(),
             ph.locationCode.c_str(),
             getBandAndInstrumentCodes(ph.channelCode).c_str(),
             UTCClock::toString(tw.startTime()).c_str(), durToSec(tw.length()));
  }

  return trace;
}

std::shared_ptr<Trace>
BasicProcessor::loadAndProcess(const TimeWindow &tw,
                               Catalog::Phase ph, // copied
                               const string &channelCode,
                               const string &filterStr,
                               double resampleFreq) const
{
  ph.channelCode = channelCode;

  const TimeWindow twToLoad =
      tw.extend(secToDur(_extraTraceLen), secToDur(_extraTraceLen));

  shared_ptr<const Trace> trace = _auxLdr->get(twToLoad, ph);
  if (!trace) return nullptr;

  shared_ptr<Trace> copy(new Trace(*trace));
  try
  {
    filter(*copy, true, filterStr, resampleFreq);
  }
  catch (exception &e)
  {
    logDebug("Errow while filtering waveform: %s", e.what());
    return nullptr;
  }

  if (!copy->slice(tw))
  {
    logDebug("Error while processing phase data '%s': cannot slice trace from "
             "%s length %.2f sec. Trace data from %s length %.2f sec, samples "
             "%zu sampfreq %f",
             string(ph).c_str(),
             HDD::UTCClock::toString(tw.startTime()).c_str(),
             HDD::durToSec(tw.length()),
             HDD::UTCClock::toString(copy->startTime()).c_str(),
             HDD::durToSec(copy->timeWindow().length()), copy->sampleCount(),
             copy->samplingFrequency());
    return nullptr;
  }
  return copy;
}

void BasicProcessor::filter(Trace &trace,
                            bool demeaning,
                            const string &filterStr,
                            double resampleFreq) const
{
  double *data           = trace.data();
  const size_t data_size = trace.sampleCount();

  if (demeaning)
  {
    double mean = std::accumulate(data, data + data_size, 0.0) / data_size;
    for (size_t i = 0; i < data_size; i++) data[i] -= mean;
  }

  if (!filterStr.empty())
  {
    _wf->filter(trace, filterStr);
  }

  if (resampleFreq > 0)
  {
    resample(trace, resampleFreq);
  }
}

shared_ptr<const Trace> MemCachedProc::get(const TimeWindow &tw,
                                           const Catalog::Phase &ph,
                                           const Catalog::Event &ev,
                                           const Catalog::Station &sta,
                                           const std::string &filterStr,
                                           double resampleFreq,
                                           Transform trans)
{
  const string wfId             = waveformId(tw, ph);
  shared_ptr<const Trace> trace = getFromCache(wfId);
  bool isCached;
  if (trace)
  {
    isCached = true;
  }
  else
  {
    isCached = false;

    // Check if we have already excluded the trace because we couldn't load it
    // (-> save time).
    if (_unloadables.find(wfId) == _unloadables.end())
    {
      trace = _auxPrc->get(tw, ph, ev, sta, filterStr, resampleFreq, trans);
      if (!trace) _unloadables.insert(wfId);
    }
  }

  if (trace && !isCached)
  {
    storeInCache(wfId, trace);
  }
  return trace;
}

bool MemCachedProc::isCached(const TimeWindow &tw,
                             const Catalog::Phase &ph,
                             const Catalog::Event &ev) const
{
  const string wfId = waveformId(tw, ph);
  const auto it     = _traces.find(wfId);
  return it != _traces.end();
}

shared_ptr<const Trace> MemCachedProc::getFromCache(const string &wfId)
{
  const auto it = _traces.find(wfId);
  return it != _traces.end() ? it->second : nullptr;
}

void MemCachedProc::storeInCache(const string &wfId,
                                 const shared_ptr<const Trace> &trace)
{
  _traces[wfId] = trace;
}

TimeWindow SnrFilterPrc::snrTimeWindow(const UTCTime &pickTime) const
{
  UTCTime winStart = std::min({pickTime + secToDur(_snr.noiseStart),
                               pickTime + secToDur(_snr.signalStart)});
  UTCTime winEnd   = std::max({pickTime + secToDur(_snr.noiseEnd),
                             pickTime + secToDur(_snr.signalEnd)});
  return TimeWindow(winStart, winEnd);
}

bool SnrFilterPrc::goodSnr(const Trace &trace, const UTCTime &pickTime) const
{
  double snr = computeSnr(trace, pickTime, _snr.noiseStart, _snr.noiseEnd,
                          _snr.signalStart, _snr.signalEnd);
  return snr >= _snr.minSnr;
}

shared_ptr<const Trace> SnrFilterPrc::get(const TimeWindow &tw,
                                          const Catalog::Phase &ph,
                                          const Catalog::Event &ev,
                                          const Catalog::Station &sta,
                                          const std::string &filterStr,
                                          double resampleFreq,
                                          Transform trans)
{
  const string wfId = waveformId(tw, ph);

  const TimeWindow twToLoad = snrTimeWindow(ph.time).merge(tw);
  shared_ptr<const Trace> trace =
      _auxPrc->get(twToLoad, ph, ev, sta, filterStr, resampleFreq, trans);

  if (!trace) return nullptr;

  if (_enabled && !goodSnr(*trace, ph.time))
  {
    logDebug("Trace has too low SNR (%s)", string(ph).c_str());
    _counters_wf_snr_low++;
    return nullptr;
  }

  if (twToLoad != tw)
  {
    shared_ptr<Trace> nonConstTrace(new Trace(*trace));
    if (!nonConstTrace->slice(tw))
    {
      logDebug("Error while checking SNR for phase '%s': cannot slice trace "
               "from %s length %.2f sec. Trace "
               "data from %s length %.2f sec, samples %zu sampfreq %f",
               string(ph).c_str(),
               HDD::UTCClock::toString(tw.startTime()).c_str(),
               HDD::durToSec(tw.length()),
               HDD::UTCClock::toString(trace->startTime()).c_str(),
               HDD::durToSec(trace->timeWindow().length()),
               trace->sampleCount(), trace->samplingFrequency());
      return nullptr;
    }
    trace = nonConstTrace;
  }

  return trace;
}

unique_ptr<Trace> transformL2(const TimeWindow &tw,
                              const Catalog::Phase &ph,
                              const ThreeComponents &tc,
                              const Trace trH1,
                              const Trace trH2)
{
  string channelCodeRoot = Waveform::getBandAndInstrumentCodes(ph.channelCode);

  if (trH1.samplingFrequency() != trH2.samplingFrequency())
  {
    logDebug(
        "Cannot perform L2 transformation with incompatible horizontal traces");
    return nullptr;
  }

  //
  // Gain
  //
  double h1scaler = 1.0;
  if (tc.gain[ThreeComponents::FirstHorizontal] != 0)
  {
    h1scaler = 1.0 / tc.gain[ThreeComponents::FirstHorizontal];
  }
  double h2scaler = 1.0;
  if (tc.gain[ThreeComponents::SecondHorizontal] != 0)
  {
    h2scaler = 1.0 / tc.gain[ThreeComponents::SecondHorizontal];
  }

  //
  // Find common start/end times
  //
  UTCTime startTime = std::max(trH1.startTime(), trH2.startTime());
  UTCTime endTime   = std::min(trH1.endTime(), trH2.endTime());
  size_t sampleCount =
      std::floor(durToSec(endTime - startTime) * trH1.samplingFrequency()) + 1;

  if (sampleCount <= 0)
  {
    logDebug("Cannot perform L2 transformation: traces do not overlap");
    return nullptr;
  }

  double h1Start = std::round(trH1.index(startTime));
  if (h1Start < 0 || (h1Start + sampleCount) > trH1.sampleCount())
  {
    logDebug("Cannot perform L2 transformation: internal logic error");
    return nullptr;
  }

  double h2Start = std::round(trH2.index(startTime));
  if (h2Start < 0 || (h2Start + sampleCount) > trH2.sampleCount())
  {
    logDebug("Cannot perform L2 transformation: internal logic error");
    return nullptr;
  }

  const double *h1data = trH1.data() + size_t(h1Start);
  const double *h2data = trH2.data() + size_t(h2Start);

  //
  // Compute L2
  //
  vector<double> l2vec(sampleCount);
  double *l2data = l2vec.data();
  for (size_t i = 0; i < sampleCount; ++i)
  {
    l2data[i] =
        std::sqrt(square(h1data[i] * h1scaler) + square(h2data[i] * h2scaler));
  }

  return unique_ptr<Trace>(new Trace(
      ph.networkCode, ph.stationCode, ph.locationCode, channelCodeRoot + "H",
      startTime, trH1.samplingFrequency(), std::move(l2vec)));
}

unique_ptr<Trace> transformRT(const TimeWindow &tw,
                              const Catalog::Phase &ph,
                              const Catalog::Event &ev,
                              const Catalog::Station &sta,
                              const ThreeComponents &tc,
                              const Trace trV,
                              const Trace trH1,
                              const Trace trH2,
                              Transform trans)
{
  string channelCodeRoot = Waveform::getBandAndInstrumentCodes(ph.channelCode);

  if (trH1.samplingFrequency() != trH2.samplingFrequency() ||
      trH1.samplingFrequency() != trV.samplingFrequency())
  {
    logDebug("Cannot perform RT transformation with incompatible traces");
    return nullptr;
  }

  using Matrix3d = double[3][3];

  struct Vector3d
  {
    double x, y, z;
    Vector3d &fromAngles(double radAzimuth, double radDip)
    {
      x = cos(radDip) * sin(radAzimuth);
      y = cos(radDip) * cos(radAzimuth);
      z = sin(radDip);
      return *this;
    }
    Vector3d &normalize()
    {
      double length = sqrt(x * x + y * y + z * z);
      double scale  = 1.0 / length;
      x *= scale;
      y *= scale;
      z *= scale;
      return *this;
    }
  };

  //
  // Rotation matrix ZNE
  //
  Vector3d h2Col;
  h2Col
      .fromAngles(+degToRad(tc.azimuth[ThreeComponents::SecondHorizontal]),
                  -degToRad(tc.dip[ThreeComponents::SecondHorizontal]))
      .normalize();
  Vector3d h1Col;
  h1Col
      .fromAngles(+degToRad(tc.azimuth[ThreeComponents::FirstHorizontal]),
                  -degToRad(tc.dip[ThreeComponents::FirstHorizontal]))
      .normalize();
  Vector3d vCol;
  vCol.fromAngles(+degToRad(tc.azimuth[ThreeComponents::Vertical]),
                  -degToRad(tc.dip[ThreeComponents::Vertical]))
      .normalize();

  Matrix3d orientationZNE = {
      h2Col.x, h1Col.x, vCol.x,  h2Col.y, h1Col.y,
      vCol.y,  h2Col.z, h1Col.z, vCol.z,
  };

  //
  // Rotation matrix ZRT
  //
  double backAzimuth;
  computeDistance(ev.latitude, ev.longitude, sta.latitude, sta.longitude,
                  nullptr, &backAzimuth);
  backAzimuth += 180.0;
  Matrix3d orientationZRT = {
      cos(backAzimuth),
      -sin(backAzimuth),
      0,
      sin(backAzimuth),
      cos(backAzimuth),
      0,
      0,
      0,
      1,
  };

  //
  // Full transformation matrix
  //
  auto matrixMult = [](const Matrix3d &m1, const Matrix3d &m2, Matrix3d &out) {
    for (int r = 0; r < 3; ++r)
      for (int col = 0; col < 3; ++col)
        out[r][col] = m1[r][0] * m2[0][col] + m1[r][1] * m2[1][col] +
                      m1[r][2] * m2[2][col];
  };
  Matrix3d transformation;
  matrixMult(orientationZRT, orientationZNE, transformation);

  //
  // Gain
  //
  double h1scaler = 1.0;
  if (tc.gain[ThreeComponents::FirstHorizontal] != 0)
  {
    h1scaler = 1.0 / tc.gain[ThreeComponents::FirstHorizontal];
  }
  double h2scaler = 1.0;
  if (tc.gain[ThreeComponents::SecondHorizontal] != 0)
  {
    h2scaler = 1.0 / tc.gain[ThreeComponents::SecondHorizontal];
  }
  double vscaler = 1.0;
  if (tc.gain[ThreeComponents::Vertical] != 0)
  {
    vscaler = 1.0 / tc.gain[ThreeComponents::Vertical];
  }

  //
  // Find common start/end times
  //
  UTCTime startTime =
      std::max({trH1.startTime(), trH2.startTime(), trV.startTime()});
  UTCTime endTime = std::min({trH1.endTime(), trH2.endTime(), trV.endTime()});
  size_t sampleCount =
      std::floor(durToSec(endTime - startTime) * trH1.samplingFrequency()) + 1;

  if (sampleCount <= 0)
  {
    logDebug("Cannot perform L2 transformation: traces do not overlap");
    return nullptr;
  }

  double h1Start = std::round(trH1.index(startTime));
  if (h1Start < 0 || (h1Start + sampleCount) > trH1.sampleCount())
  {
    logDebug("Cannot perform RT transformation: internal logic error");
    return nullptr;
  }

  double h2Start = std::round(trH2.index(startTime));
  if (h2Start < 0 || (h2Start + sampleCount) > trH2.sampleCount())
  {
    logDebug("Cannot perform RT transformation: internal logic error");
    return nullptr;
  }

  double vStart = std::round(trV.index(startTime));
  if (vStart < 0 || (vStart + sampleCount) > trV.sampleCount())
  {
    logDebug("Cannot perform RT transformation: internal logic error");
    return nullptr;
  }

  const double *h1data = trH1.data() + size_t(h1Start);
  const double *h2data = trH2.data() + size_t(h2Start);
  const double *vdata  = trV.data() + size_t(vStart);

  //
  // Apply RT transformation now
  //
  vector<double> rtvec(sampleCount);
  double *rtdata = rtvec.data();
  for (size_t i = 0; i < sampleCount; ++i)
  {
    if (trans == Transform::TRANSVERSAL) // T trasversal
      rtdata[i] = transformation[0][0] * h2data[i] * h2scaler +
                  transformation[0][1] * h1data[i] * h1scaler +
                  transformation[0][2] * vdata[i] * vscaler;
    else if (trans == Transform::RADIAL) // R radial
      rtdata[i] = transformation[1][0] * h2data[i] * h2scaler +
                  transformation[1][1] * h1data[i] * h1scaler +
                  transformation[1][2] * vdata[i] * vscaler;
    else // V vertical
      rtdata[i] = transformation[2][0] * h2data[i] * h2scaler +
                  transformation[2][1] * h1data[i] * h1scaler +
                  transformation[2][2] * vdata[i] * vscaler;
  }

  string channelCode = channelCodeRoot;
  if (trans == Transform::TRANSVERSAL)
    channelCode += "T";
  else if (trans == Transform::RADIAL)
    channelCode += "R";
  else
    channelCode += "V";
  return unique_ptr<Trace>(
      new Trace(ph.networkCode, ph.stationCode, ph.locationCode, channelCode,
                trH1.startTime(), trH1.samplingFrequency(), std::move(rtvec)));
}

void resample(Trace &trace, double new_sf)
{
  if (new_sf <= 0 || trace.samplingFrequency() == new_sf) return;

  const double data_sf       = trace.samplingFrequency();
  const double resamp_factor = new_sf / data_sf;
  const double nyquist       = std::min(new_sf, data_sf) / 2.;

  const double *data     = trace.data();
  const size_t data_size = trace.sampleCount();

  vector<double> resampled(data_size * resamp_factor);

  /*
   * Compute one sample of the resampled data:
   *
   * x: new sample point location (relative to old indexes)
   *    (e.g. every other integer for 0.5x decimation)
   * fmax: low pass filter cutoff frequency. Fmax should be less
   *       than half of data_freq, and less than half of the
   *       new sample frequency (the reciprocal of the x step size).
   * win_len: width of windowed Sinc used as the low pass filter
   *          Filter quality increases with a larger window width.
   *          The wider the window, the closer fmax can approach half of
   *          data_freq or the new sample frequency
   *
   * If the x step size is rational the same Window and Sinc values
   * will be recalculated repeatedly. Therefore these values can either
   * be cached, or pre-calculated and stored in a table (polyphase
   * interpolation); or interpolated from a smaller pre-calculated table;
   * or computed from a set of low-order polynomials fitted to each
   * section or lobe between  zero-crossings of the windowed Sinc (Farrow)
   *
   * Credits: Ronald Nicholson, "Ron's Digital Signal Processing Page"
   */
  auto new_sample = [&data, &data_size, &data_sf](double x, double fmax,
                                                  int win_len) -> double {
    const double gain = 2 * fmax / data_sf; // Calc gain correction factor
    double newSmp     = 0;

    // For 1 window width
    for (int win_i = -(win_len / 2.); win_i < (win_len / 2.); win_i++)
    {
      const int j        = int(x + win_i); // input sample index
      const double win_x = j - x;
      if (j >= 0 && size_t(j) < data_size)
      {
        // calculate von Hann Window | hann(x) = sin^2(pi*x/N)
        const double hannWin = square(std::sin(M_PI * (0.5 + win_x / win_len)));

        // Scale and calculate Sinc | sinc(x) = sin(pi*x)/(pi*x); sinc(0)=1
        const double a    = M_PI * win_x * gain;
        const double sinc = (a == 0) ? 1 : std::sin(a) / a;

        newSmp += gain * hannWin * sinc * data[j];
      }
    }
    return newSmp;
  };

  double *resampled_data = resampled.data();
  for (size_t i = 0; i < resampled.size(); i++)
  {
    double x          = i / resamp_factor;
    resampled_data[i] = new_sample(x, nyquist, 21);
  }

  trace.setData(std::move(resampled));
  trace.setSamplingFrequency(new_sf);
}

double computeSnr(const Trace &tr,
                  const UTCTime &pickTime,
                  double noiseOffsetStart,
                  double noiseOffsetEnd,
                  double signalOffsetStart,
                  double signalOffsetEnd)
{
  const double *data     = tr.data();
  const size_t data_size = tr.sampleCount();

  const double noiseStartIdx =
      std::round(tr.index(pickTime + secToDur(noiseOffsetStart)));
  const double noiseEndIdx =
      std::round(tr.index(pickTime + secToDur(noiseOffsetEnd)));
  const double signalStartIdx =
      std::round(tr.index(pickTime + secToDur(signalOffsetStart)));
  const double signalEndIdx =
      std::round(tr.index(pickTime + secToDur(signalOffsetEnd)));

  if ((std::min({noiseStartIdx, noiseEndIdx, signalStartIdx, signalEndIdx}) <
       0) ||
      (std::max({noiseStartIdx, noiseEndIdx, signalStartIdx, signalEndIdx}) >=
       data_size))
  {
    logDebug(
        "Cannot compute SNR: noise/signal windows exceed waveform boundaries");
    return -1;
  }

  const size_t noiseStart  = noiseStartIdx;
  const size_t noiseEnd    = noiseEndIdx;
  const size_t signalStart = signalStartIdx;
  const size_t signalEnd   = signalEndIdx;

  double noise = 0;
  for (size_t i = noiseStart; i < noiseEnd; i++)
  {
    noise += square(data[i]);
  }
  noise /= (noiseEnd - noiseStart);

  double signal = 0;
  for (size_t i = signalStart; i < signalEnd; i++)
  {
    signal += square(data[i]);
  }
  signal /= (signalEnd - signalStart);

  return signal / noise;
}

} // namespace Waveform
} // namespace HDD

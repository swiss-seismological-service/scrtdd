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

#include "waveform.h"

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iostream>
#include <seiscomp3/client/inventory.h>
#include <seiscomp3/core/datetime.h>
#include <seiscomp3/datamodel/station.h>
#include <seiscomp3/io/recordinput.h>
#include <seiscomp3/io/records/mseedrecord.h>
#include <seiscomp3/math/filter.h>
#include <seiscomp3/math/geo.h>
#include <seiscomp3/processing/operator/ncomps.h>
#include <seiscomp3/processing/operator/transformation.h>
#include <seiscomp3/utils/files.h>

#define SEISCOMP_COMPONENT HDD
#include <seiscomp3/logging/log.h>

using namespace std;
using namespace Seiscomp;
using namespace Seiscomp::Processing;
using DataModel::ThreeComponents;
using Seiscomp::Core::stringify;
using Catalog = HDD::Catalog;

namespace {

template <class T> T nextPowerOf2(T a, T min = 1, T max = 1 << 31)
{
  int b = min;
  while (b < a)
  {
    b <<= 1;
    if (b > max) return -1;
  }
  return b;
}

string waveformDebugPath(const string &wfDebugDir,
                         const Catalog::Event &ev,
                         const Catalog::Phase &ph,
                         const std::string &ext)
{
  string debugFile =
      stringify("ev%u.%s.%s.%s.%s.%s.%s.mseed", ev.id, ph.networkCode.c_str(),
                ph.stationCode.c_str(), ph.locationCode.c_str(),
                ph.channelCode.c_str(), ph.type.c_str(), ext.c_str());
  return (boost::filesystem::path(wfDebugDir) / debugFile).string();
}

} // namespace

namespace Seiscomp {
namespace HDD {
namespace Waveform {

std::string getBandAndInstrumentCodes(const std::string &channelCode)
{
  if (channelCode.size() >= 2) return channelCode.substr(0, 2);
  return "";
}

std::string getOrientationCode(const std::string &channelCode)
{
  if (channelCode.size() == 3) return channelCode.substr(2, 3);
  return "";
}

std::string waveformId(const Core::TimeWindow &tw,
                       const std::string &networkCode,
                       const std::string &stationCode,
                       const std::string &locationCode,
                       const std::string &channelCode)
{
  return Core::stringify("%s.%s.%s.%s.%s.%s", networkCode.c_str(),
                         stationCode.c_str(), locationCode.c_str(),
                         channelCode.c_str(), tw.startTime().iso().c_str(),
                         tw.endTime().iso().c_str());
}

std::string waveformId(const HDD::Catalog::Phase &ph,
                       const Core::TimeWindow &tw)
{
  return waveformId(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                    ph.channelCode);
}

GenericRecordPtr readWaveformFromRecordStream(const string &recordStreamURL,
                                              const Core::TimeWindow &tw,
                                              const string &networkCode,
                                              const string &stationCode,
                                              const string &locationCode,
                                              const string &channelCode)
{
  IO::RecordStreamPtr rs = IO::RecordStream::Open(recordStreamURL.c_str());
  if (rs == nullptr)
  {
    string msg = "Cannot open RecordStream: " + recordStreamURL;
    throw runtime_error(msg);
  }

  rs->setTimeWindow(tw);
  rs->addStream(networkCode, stationCode, locationCode, channelCode);

  // Store each record in a RecordSequence
  IO::RecordInput inp(rs.get(), Array::DOUBLE, Record::DATA_ONLY);
  std::shared_ptr<RecordSequence> seq(new TimeWindowBuffer(tw));
  RecordPtr rec;
  while (rec = inp.next()) { seq->feed(rec.get()); }
  rs->close();

  if (seq->empty())
  {
    string msg = stringify(
        "Data could not be loaded (stream %s.%s.%s.%s from %s length %.2f sec)",
        networkCode.c_str(), stationCode.c_str(), locationCode.c_str(),
        channelCode.c_str(), tw.startTime().iso().c_str(), tw.length());
    throw runtime_error(msg);
  }

  GenericRecordPtr trace = new GenericRecord();

  if (!merge(*trace, *seq))
  {
    string msg = stringify(
        "Data records could not be merged into a single trace "
        "(%s.%s.%s.%s from %s length %.2f sec)",
        networkCode.c_str(), stationCode.c_str(), locationCode.c_str(),
        channelCode.c_str(), tw.startTime().iso().c_str(), tw.length());
    throw runtime_error(msg);
  }

  if (!trim(*trace, tw))
  {
    string msg = stringify("Incomplete trace, not enough data for requested"
                           " time window (%s.%s.%s.%s from %s length %.2f sec)",
                           networkCode.c_str(), stationCode.c_str(),
                           locationCode.c_str(), channelCode.c_str(),
                           tw.startTime().iso().c_str(), tw.length());
    throw runtime_error(msg);
  }

  return trace;
}

bool merge(GenericRecord &trace, const RecordSequence &seq)
{
  if (seq.empty()) { return false; }

  RecordCPtr first = seq.front();
  RecordCPtr last;
  double samplingFrequency = first->samplingFrequency();
  Core::TimeSpan maxAllowedGap, maxAllowedOverlap;

  maxAllowedGap     = Core::TimeSpan((double)(0.5 / samplingFrequency));
  maxAllowedOverlap = Core::TimeSpan((double)(-0.5 / samplingFrequency));

  trace.setNetworkCode(first->networkCode());
  trace.setStationCode(first->stationCode());
  trace.setLocationCode(first->locationCode());
  trace.setChannelCode(first->channelCode());

  trace.setStartTime(first->startTime());
  trace.setSamplingFrequency(samplingFrequency);

  Array::DataType datatype = first->data()->dataType();
  ArrayPtr arr = ArrayFactory::Create(datatype, datatype, 0, nullptr);

  for (const RecordCPtr &rec : seq)
  {
    if (rec->samplingFrequency() != samplingFrequency)
    {
      SEISCOMP_DEBUG(
          "%s.%s.%s.%s: record sampling frequencies are not consistent: %f != "
          "%f",
          trace.networkCode().c_str(), trace.stationCode().c_str(),
          trace.locationCode().c_str(), trace.channelCode().c_str(),
          samplingFrequency, rec->samplingFrequency());
      return false;
    }

    // Check for gaps and overlaps
    if (last)
    {
      Core::TimeSpan diff = rec->startTime() - last->endTime();
      if (diff > maxAllowedGap)
      {
        SEISCOMP_DEBUG("%s.%s.%s.%s: gap detected of %d.%06ds",
                       trace.networkCode().c_str(), trace.stationCode().c_str(),
                       trace.locationCode().c_str(),
                       trace.channelCode().c_str(), (int)diff.seconds(),
                       (int)diff.microseconds());
        return false;
      }

      if (diff < maxAllowedOverlap)
      {
        SEISCOMP_DEBUG("%s.%s.%s.%s: overlap detected of %fs",
                       trace.networkCode().c_str(), trace.stationCode().c_str(),
                       trace.locationCode().c_str(),
                       trace.channelCode().c_str(), (double)diff);
        return false;
      }
    }

    arr->append((Array *)(rec->data()));

    last = rec;
  }

  trace.setData(arr.get());

  return true;
}

bool trim(GenericRecord &trace, const Core::TimeWindow &tw)
{
  int ofs     = (int)(double(tw.startTime() - trace.startTime()) *
                  trace.samplingFrequency());
  int samples = (int)(tw.length() * trace.samplingFrequency());

  // Not enough data at start of time window
  if (ofs < 0)
  {
    SEISCOMP_DEBUG("%s: need %d more samples in past", trace.streamID().c_str(),
                   -ofs);
    return false;
  }

  // Not enough data at end of time window
  if (ofs + samples > trace.data()->size())
  {
    SEISCOMP_DEBUG("%s: need %d more samples past the end",
                   trace.streamID().c_str(),
                   trace.data()->size() - samples - ofs);
    return false;
  }

  ArrayPtr sliced = trace.data()->slice(ofs, ofs + samples);

  trace.setStartTime(tw.startTime());
  trace.setData(sliced.get());

  return true;
}

void filter(GenericRecord &trace,
            bool demeaning,
            const std::string &filterStr,
            double resampleFreq)
{
  DoubleArray *data = DoubleArray::Cast(trace.data());

  if (demeaning)
  {
    *data -= data->mean();
    trace.dataUpdated();
  }

  if (resampleFreq > 0) { resample(trace, resampleFreq, true); }

  if (!filterStr.empty())
  {
    string filterError;
    auto filter =
        Math::Filtering::InPlaceFilter<double>::Create(filterStr, &filterError);
    if (!filter)
    {
      string msg = stringify("Filter creation failed %s: %s", filterStr.c_str(),
                             filterError.c_str());
      throw runtime_error(msg);
    }
    filter->setSamplingFrequency(trace.samplingFrequency());
    filter->apply(data->size(), data->typedData());
    delete filter;
    trace.dataUpdated();
  }
}

void resample(GenericRecord &trace, double sf, bool average)
{
  if (sf <= 0) return;

  if (trace.samplingFrequency() == sf) return;

  DoubleArray *data = DoubleArray::Cast(trace.data());
  double step       = trace.samplingFrequency() / sf;

  if (trace.samplingFrequency() < sf) // upsampling
  {
    double fi = data->size() - 1;
    data->resize(data->size() / step);

    for (int i = data->size() - 1; i >= 0; i--)
    {
      (*data)[i] = (*data)[(int)fi];
      fi -= step;
    }
  }
  else // downsampling
  {
    int w     = average ? step * 0.5 + 0.5 : 0;
    int i     = 0;
    double fi = 0.0;
    int cnt   = data->size();

    if (w <= 0)
    {
      while (fi < cnt)
      {
        (*data)[i++] = (*data)[(int)fi];
        fi += step;
      }
    }
    else
    {
      while (fi < cnt)
      {
        int ci       = (int)fi;
        double scale = 1.0;
        double v     = (*data)[ci];

        for (int g = 1; g < w; ++g)
        {
          if (ci >= g)
          {
            v += (*data)[ci - g];
            scale += 1.0;
          }

          if (ci + g < cnt)
          {
            v += (*data)[ci + g];
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

void writeTrace(GenericRecordCPtr trace, const std::string &file)
{
  if (!trace) return;

  try
  {
    std::ofstream ofs(file);
    IO::MSeedRecord msRec(*trace);
    int reclen = msRec.data()->size() * msRec.data()->bytes() + 64;
    reclen     = nextPowerOf2<int>(reclen, 128,
                               1048576); // MINRECLEN 128, MAXRECLEN 1048576
    if (reclen > 0)
    {
      msRec.setOutputRecordLength(reclen);
      msRec.write(ofs);
    }
  }
  catch (exception &e)
  {
    SEISCOMP_WARNING("Couldn't write waveform to disk %s: %s", file.c_str(),
                     e.what());
  }
}

GenericRecordPtr readTrace(const std::string &file)
{
  if (!Util::fileExists(file)) return nullptr;

  try
  {
    std::ifstream ifs(file);
    IO::MSeedRecord msRec(Array::DOUBLE, Record::Hint::DATA_ONLY);
    msRec.read(ifs);
    GenericRecordPtr trace = new GenericRecord(msRec);
    trace->setData(msRec.data()->clone()); // copy data too
    return trace;
  }
  catch (exception &e)
  {
    SEISCOMP_WARNING("Couldn't load waveform %s: %s", file.c_str(), e.what());
    return nullptr;
  }
}

/*
 * Calculate the correlation series (tr1 and tr2 are already demeaned)
 *
 * delayOut is the shift in seconds (positive or negative) between tr1 and tr2
 * middle points to get the highest correlation coefficient (coeffOut) between
 * the 2 traces A delayOut of 0 is when tr1 and tr2 middle points are aligned
 */
bool xcorr(const GenericRecordCPtr &tr1,
           const GenericRecordCPtr &tr2,
           double maxDelay,
           bool qualityCheck,
           double &delayOut,
           double &coeffOut)
{
  coeffOut = std::nan("");

  if (tr1->samplingFrequency() != tr2->samplingFrequency())
  {
    SEISCOMP_INFO(
        "Cannot cross correlate traces with different sampling freq (%f!=%f)",
        tr1->samplingFrequency(), tr2->samplingFrequency());
    return false;
  }

  const double freq      = tr1->samplingFrequency();
  const int maxDelaySmps = maxDelay * freq; // secs to samples

  // check longest/shortest trace
  const bool swap             = tr1->data()->size() > tr2->data()->size();
  GenericRecordCPtr trShorter = swap ? tr2 : tr1;
  GenericRecordCPtr trLonger  = swap ? tr1 : tr2;

  const double *smpsS = DoubleArray::ConstCast(trShorter->data())->typedData();
  const double *smpsL = DoubleArray::ConstCast(trLonger->data())->typedData();
  const int smpsSsize = trShorter->data()->size();
  const int smpsLsize = trLonger->data()->size();

  //
  // for later quality check: save local maxima/minima
  //
  struct LocalMaxima
  {
    bool notDecreasing = false;
    double prevCoeff   = -1;
    vector<double> values;
    void update(double coeff)
    {
      if (!std::isfinite(coeff)) return;
      if (coeff < prevCoeff && notDecreasing) values.push_back(prevCoeff);
      notDecreasing = coeff >= prevCoeff;
      prevCoeff     = coeff;
    }
  };
  LocalMaxima localMaxs, localMins;

  auto sampleAtLong = [&, smpsSsize, smpsLsize, smpsL](int idxS, int delay) {
    int idxL = idxS + (smpsLsize - smpsSsize) / 2 + delay;
    return (idxL < 0 || idxL >= smpsLsize) ? 0 : smpsL[idxL];
  };

  // do as much computation as possible outside the main xcorr loop
  double denomL = 0, denomS = 0;
  for (int idxS = 0; idxS < smpsSsize; idxS++)
  {
    denomS += smpsS[idxS] * smpsS[idxS];
    double sampleL = sampleAtLong(idxS, -(maxDelaySmps+1));
    denomL += sampleL * sampleL;
  }

  // cross-correlation loop
  for (int delay = -maxDelaySmps; delay < maxDelaySmps; delay++)
  {
    // remove from denomL the sample that has exited the current xcorr win
    double lastSampleL = sampleAtLong(-1, delay);
    denomL -= lastSampleL * lastSampleL;
    // add to denomL the sample that has just entered the current xcorr win
    double newSampleL = sampleAtLong(smpsSsize-1, delay);
    denomL += newSampleL * newSampleL;

    // compute numerator
    double numer = 0;
    for (int idxS = 0; idxS < smpsSsize; idxS++)
      numer += smpsS[idxS] * sampleAtLong(idxS, delay);

    const double denom = std::sqrt(denomS * denomL);
    const double coeff = numer / denom;

    if (std::abs(coeff) > std::abs(coeffOut) || !std::isfinite(coeffOut))
    {
      coeffOut = coeff;
      delayOut = delay / freq; // samples to secs
    }

    // for later quality check
    localMaxs.update(coeff);
    localMins.update(-coeff);
  }

  if (swap) { delayOut = -delayOut; }

  /*
   * To avoid errors introduced by cycle skipping the differential time
   * measurement is only accepted, if all side lobe maxima CCslm of the
   * cross-correlation function fulfill the following condition:
   *
   *                CCslm < CCmax - ( (1.0-CCmax) / 2.0 )
   *
   * where CCmax corresponds to the global maximum of the cross-correlation
   * function. By discarding measurements with local maxima CCslm close to the
   * global maximum CC, the number of potential blunders due to cycle skipping
   * is significantly reduced.
   *
   * See Diehl et al. (2017): The induced earthquake sequence related to the St.
   * Gallen deep geothermal project: Fault reactivation and fluid interactions
   * imaged by microseismicity
   * */
  if (qualityCheck && std::isfinite(coeffOut))
  {
    double threshold = std::abs(coeffOut) - ((1.0 - std::abs(coeffOut)) / 2.0);
    int numMax       = 0;
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

  if (!std::isfinite(coeffOut))
  {
    coeffOut = 0;
    delayOut = 0.;
  }

  return true;
}

double computeSnr(const GenericRecordCPtr &tr,
                  const Core::Time &pickTime,
                  double noiseOffsetStart,
                  double noiseOffsetEnd,
                  double signalOffsetStart,
                  double signalOffsetEnd)
{
  const double *data = DoubleArray::ConstCast(tr->data())->typedData();
  const int size     = tr->data()->size();
  const double freq  = tr->samplingFrequency();
  const Core::Time dataStartTime = tr->startTime();

  // convert time w.r.t. guiding pick time to sample number
  auto secToSample = [&, freq, size](double sec) {
    return std::min(std::max(std::round(sec * freq), 0.), size - 1.);
  };
  const double pickOffset = (pickTime - dataStartTime).length();
  const int noiseStart    = secToSample(noiseOffsetStart + pickOffset);
  const int noiseEnd      = secToSample(noiseOffsetEnd + pickOffset);
  const int signalStart   = secToSample(signalOffsetStart + pickOffset);
  const int signalEnd     = secToSample(signalOffsetEnd + pickOffset);

  if ((std::min({noiseStart, noiseEnd, signalStart, signalEnd}) < 0) ||
      (std::max({noiseStart, noiseEnd, signalStart, signalEnd}) >= size))
  {
    SEISCOMP_ERROR(
        "Cannot compute SNR: noise/signal windows exceed waveform boundaries");
    return -1;
  }

  // Get maximum (absolute) amplitude in noise window:
  double noiseMax = -1.0;
  for (int i = noiseStart; i < noiseEnd; i++)
  { noiseMax = std::max(std::abs(data[i]), noiseMax); }

  // Get maximum (absolute) amplitude in signal window:
  double signalMax = -1.0;
  for (int i = signalStart; i < signalEnd; i++)
  { signalMax = std::max(std::abs(data[i]), signalMax); }

  return signalMax / noiseMax;
}

GenericRecordPtr Loader::readAndProjectWaveform(const Core::TimeWindow &tw,
                                                const Catalog::Phase &ph,
                                                const Catalog::Event &ev)
{
  string wfDesc = stringify("Waveform Projection for '%s'", string(ph).c_str());

  string channelCodeRoot = getBandAndInstrumentCodes(ph.channelCode);
  string component       = getOrientationCode(ph.channelCode);
  bool allComponents     = false;
  DataModel::ThreeComponents tc;
  DataModel::SensorLocation *loc = Catalog::findSensorLocation(
      ph.networkCode, ph.stationCode, ph.locationCode, tw.startTime());

  if (loc)
  {
    //
    // Check if the projection is actually required, return nullpt and
    // do not throw if that's not the case
    //

    allComponents =
        getThreeComponents(tc, loc, channelCodeRoot.c_str(), tw.startTime());

    if ((tc.comps[ThreeComponents::Vertical] &&
         tc.comps[ThreeComponents::Vertical]->code() == ph.channelCode) ||
        (tc.comps[ThreeComponents::FirstHorizontal] &&
         tc.comps[ThreeComponents::FirstHorizontal]->code() ==
             ph.channelCode) ||
        (tc.comps[ThreeComponents::SecondHorizontal] &&
         tc.comps[ThreeComponents::SecondHorizontal]->code() == ph.channelCode))
    {
      string msg = stringify("Projection not required (%s)", wfDesc.c_str());
      throw runtime_error(msg);
    }

    if (std::set<string>{"Z", "N", "E", "R", "T"}.count(component) == 0)
    {
      string msg = stringify("Unkwown projection component '%s' (%s)",
                             component.c_str(), wfDesc.c_str());
      throw runtime_error(msg);
    }
  }

  if (!loc || !allComponents)
  {
    string msg = stringify("Unable to fetch orientation information (%s)",
                           wfDesc.c_str());
    throw runtime_error(msg);
  }

  //
  // Load data and Perform the projection
  //
  SEISCOMP_DEBUG("Loading the 3 components waveforms (%s %s %s) to perform the "
                 "projection...",
                 tc.comps[ThreeComponents::Vertical]->code().c_str(),
                 tc.comps[ThreeComponents::FirstHorizontal]->code().c_str(),
                 tc.comps[ThreeComponents::SecondHorizontal]->code().c_str());

  // orientation ZNE
  Math::Matrix3d orientationZNE;
  Math::Vector3d n;
  n.fromAngles(+deg2rad(tc.comps[ThreeComponents::Vertical]->azimuth()),
               -deg2rad(tc.comps[ThreeComponents::Vertical]->dip()))
      .normalize();
  orientationZNE.setColumn(2, n);
  n.fromAngles(+deg2rad(tc.comps[ThreeComponents::FirstHorizontal]->azimuth()),
               -deg2rad(tc.comps[ThreeComponents::FirstHorizontal]->dip()))
      .normalize();
  orientationZNE.setColumn(1, n);
  n.fromAngles(+deg2rad(tc.comps[ThreeComponents::SecondHorizontal]->azimuth()),
               -deg2rad(tc.comps[ThreeComponents::SecondHorizontal]->dip()))
      .normalize();
  orientationZNE.setColumn(0, n);

  // orientation ZRT
  Math::Matrix3d orientationZRT;
  double delta, az, baz;
  Math::Geo::delazi(ev.latitude, ev.longitude, loc->latitude(),
                    loc->longitude(), &delta, &az, &baz);
  orientationZRT.loadRotateZ(deg2rad(baz + 180.0));

  // transformation matrix
  Math::Matrix3d transformation;
  map<string, string> chCodeMap;

  if (component == "Z" || component == "N" || component == "E")
  {
    transformation = orientationZNE;
    chCodeMap[channelCodeRoot + "Z"] =
        tc.comps[ThreeComponents::Vertical]->code();
    chCodeMap[channelCodeRoot + "N"] =
        tc.comps[ThreeComponents::FirstHorizontal]->code();
    chCodeMap[channelCodeRoot + "E"] =
        tc.comps[ThreeComponents::SecondHorizontal]->code();

    SEISCOMP_DEBUG("Performing ZNE projection (channelCode %s -> %s) for %s",
                   chCodeMap[ph.channelCode].c_str(), ph.channelCode.c_str(),
                   wfDesc.c_str());
  }
  else if (component == "R" || component == "T")
  {
    transformation.mult(orientationZRT, orientationZNE);
    // chCodeMap[channelCodeRoot + "Z"] =
    // tc.comps[ThreeComponents::Vertical]->code();
    chCodeMap[channelCodeRoot + "R"] =
        tc.comps[ThreeComponents::FirstHorizontal]->code();
    chCodeMap[channelCodeRoot + "T"] =
        tc.comps[ThreeComponents::SecondHorizontal]->code();

    SEISCOMP_DEBUG("Performing ZRT projection (channelCode %s -> %s) for %s",
                   chCodeMap[ph.channelCode].c_str(), ph.channelCode.c_str(),
                   wfDesc.c_str());
  }
  else
  {
    throw runtime_error("Internal logic error: this shouldn't happend");
  }

  //
  // Load the components
  //
  auto loadWaveform = [this](
                          const string &networkCode, const string &stationCode,
                          const string &locationCode, const string &channelCode,
                          const Core::TimeWindow &tw) {
    GenericRecordCPtr trace;
    if (_doCaching && !_cacheProcessed)
    {
      trace =
          getFromCache(tw, networkCode, stationCode, locationCode, channelCode);
      if (trace) _counters_wf_cached++;
    }
    if (!trace)
    {
      trace =
          readWaveformFromRecordStream(_recordStreamURL, tw, networkCode,
                                       stationCode, locationCode, channelCode);
      _counters_wf_downloaded++;
    }
    if (trace && _doCaching && !_cacheProcessed)
    {
      storeInCache(tw, networkCode, stationCode, locationCode, channelCode,
                   trace);
    }
    return trace;
  };
  GenericRecordCPtr tr1 =
      loadWaveform(ph.networkCode, ph.stationCode, ph.locationCode,
                   tc.comps[ThreeComponents::Vertical]->code(), tw);
  GenericRecordCPtr tr2 =
      loadWaveform(ph.networkCode, ph.stationCode, ph.locationCode,
                   tc.comps[ThreeComponents::FirstHorizontal]->code(), tw);
  GenericRecordCPtr tr3 =
      loadWaveform(ph.networkCode, ph.stationCode, ph.locationCode,
                   tc.comps[ThreeComponents::SecondHorizontal]->code(), tw);

  // The wrapper will direct 3 codes into the right slots using the
  // Stream configuration class and will finally use the transformation
  // operator. The advantage is that it will apply the configured gain for
  typedef Operator::StreamConfigWrapper<double, 3, Operator::Transformation>
      OpWrapper;

  // Define the final operator class:
  //  1. Send channel codes to right slots
  //  2. Align 3 channels sample wise
  //  3. Transform the resulting 3 component trace with a rotation matrix
  typedef NCompsOperator<double, 3, OpWrapper> Rotator;

  Processing::Stream streams[3];
  streams[2].init(tc.comps[ThreeComponents::Vertical]);
  streams[1].init(tc.comps[ThreeComponents::FirstHorizontal]);
  streams[0].init(tc.comps[ThreeComponents::SecondHorizontal]);
  Rotator op(
      OpWrapper(streams, Operator::Transformation<double, 3>(transformation)));

  class DataStorer
  {
  public:
    DataStorer(const string &channelCode,
               const Core::TimeWindow &tw,
               map<string, string> chCodeMap)
        : chMap(chCodeMap[channelCode], channelCode),
          _seq(new TimeWindowBuffer(tw))
    {}

    bool store(const Record *rec)
    {
      if (rec->channelCode() == chMap.first) _seq->feed(rec);
      return true;
    }

    const std::pair<string, string> chMap;
    std::shared_ptr<RecordSequence> _seq;
  };

  DataStorer projectedData(ph.channelCode, tw, chCodeMap);

  // The function that will be called after a transformed record was created
  // op.setStoreFunc(boost::bind(&RecordSequence::feed, seq, _1));
  op.setStoreFunc(boost::bind(&DataStorer::store, projectedData, _1));

  op.feed(tr1.get());
  op.feed(tr2.get());
  op.feed(tr3.get());

  std::shared_ptr<RecordSequence> seq = projectedData._seq;

  if (seq->empty())
  {
    string msg =
        stringify("No data after the projection for %s", wfDesc.c_str());
    throw runtime_error(msg);
  }

  GenericRecordPtr trace = new GenericRecord();

  if (!merge(*trace, *seq))
  {
    string msg =
        stringify("Data records could not be merged into a single trace (%s)",
                  wfDesc.c_str());
    throw runtime_error(msg);
  }

  trace->setChannelCode(ph.channelCode);

  if (!trim(*trace, tw))
  {
    string msg =
        stringify("Incomplete trace, not enough data (%s)", wfDesc.c_str());
    throw runtime_error(msg);
  }

  return trace;
}

GenericRecordCPtr Loader::get(const Core::TimeWindow &tw,
                              const Catalog::Phase &ph,
                              const Catalog::Event &ev,
                              bool demeaning,
                              const std::string &filterStr,
                              double resampleFreq)
{
  GenericRecordCPtr trace;
  bool isProcessed;
  bool isCached;

  if (_doCaching)
  {
    trace = getFromCache(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                         ph.channelCode);
    if (trace)
    {
      isCached    = true;
      isProcessed = _cacheProcessed;
      _counters_wf_cached++;
    }
  }

  // if the trace is not cached try loading it
  if (!trace)
  {
    isCached = false;
    if (_auxLdr)
    {
      // Load trace from the auxiliary loader
      trace = _auxLdr->get(tw, ph, ev, demeaning, filterStr, resampleFreq);
      isProcessed = true;
    }
    else if (!_recordStreamURL.empty())
    {
      // Load trace from the configured record stream
      try
      {
        trace = readWaveformFromRecordStream(_recordStreamURL, tw,
                                             ph.networkCode, ph.stationCode,
                                             ph.locationCode, ph.channelCode);
        _counters_wf_downloaded++;
      }
      catch (exception &e)
      {
        SEISCOMP_DEBUG("%s", e.what());
        try
        {
          // if the waveform is not available, possibly a projection
          // 123->ZNE or ZNE->ZRT is required
          trace = readAndProjectWaveform(tw, ph, ev);
          // raw traces are cached by readAndProjectWaveform
          if (!_cacheProcessed) isCached = true;
        }
        catch (exception &e)
        {
          SEISCOMP_DEBUG("%s", e.what());
        }
      }
      isProcessed = false;
    }
  }

  if (!trace) { _counters_wf_no_avail++; }
  else
  {
    // cache unprocessed trace
    if (_doCaching && !isCached && !_cacheProcessed && !isProcessed)
    {
      storeInCache(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                   ph.channelCode, trace);
    }

    // process trace
    if (!isProcessed)
    {
      GenericRecordPtr nonConstTrace(new GenericRecord(*trace));
      filter(*nonConstTrace, demeaning, filterStr, resampleFreq);
      isProcessed = true;
      trace       = nonConstTrace;
    }

    // cache processed trace
    if (_doCaching && !isCached && _cacheProcessed && isProcessed)
    {
      storeInCache(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                   ph.channelCode, trace);
    }

    // Debugging waveforms: first time we load a waveform dump it
    if (!_wfDebugDir.empty() && _doCaching && !isCached && isProcessed)
    {
      string ext = (ph.procInfo.source == Catalog::Phase::Source::THEORETICAL)
                       ? "theoretical"
                       : (ph.isManual ? "manual" : "automatic");
      writeTrace(trace, waveformDebugPath(_wfDebugDir, ev, ph, ext));
    }
  }

  return trace;
}

bool DiskCachedLoader::isCached(const Core::TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev)
{
  const std::string cacheFile =
      waveformPath(_cacheDir, tw, ph.networkCode, ph.stationCode,
                   ph.locationCode, ph.channelCode);
  return Util::fileExists(cacheFile);
}

GenericRecordCPtr
DiskCachedLoader::getFromCache(const Core::TimeWindow &tw,
                               const std::string &networkCode,
                               const std::string &stationCode,
                               const std::string &locationCode,
                               const std::string &channelCode)
{
  const std::string cacheFile = waveformPath(
      _cacheDir, tw, networkCode, stationCode, locationCode, channelCode);
  return readTrace(cacheFile);
}

void DiskCachedLoader::storeInCache(const Core::TimeWindow &tw,
                                    const std::string &networkCode,
                                    const std::string &stationCode,
                                    const std::string &locationCode,
                                    const std::string &channelCode,
                                    const GenericRecordCPtr &trace)
{
  const std::string cacheFile = waveformPath(
      _cacheDir, tw, networkCode, stationCode, locationCode, channelCode);
  writeTrace(trace, cacheFile);
}

std::string DiskCachedLoader::waveformPath(const std::string &cacheDir,
                                           const Core::TimeWindow &tw,
                                           const std::string &networkCode,
                                           const std::string &stationCode,
                                           const std::string &locationCode,
                                           const std::string &channelCode) const
{
  std::string cacheFile =
      waveformId(tw, networkCode, stationCode, locationCode, channelCode) +
      ".mseed";
  return (boost::filesystem::path(cacheDir) / cacheFile).string();
}

bool MemCachedLoader::isCached(const Core::TimeWindow &tw,
                               const Catalog::Phase &ph,
                               const Catalog::Event &ev)
{
  const string wfId = waveformId(ph, tw);
  const auto it     = _waveforms.find(wfId);
  return it != _waveforms.end();
}

GenericRecordCPtr MemCachedLoader::getFromCache(const Core::TimeWindow &tw,
                                                const std::string &networkCode,
                                                const std::string &stationCode,
                                                const std::string &locationCode,
                                                const std::string &channelCode)
{
  const string wfId =
      waveformId(tw, networkCode, stationCode, locationCode, channelCode);
  const auto it = _waveforms.find(wfId);
  return it != _waveforms.end() ? it->second : nullptr;
}

void MemCachedLoader::storeInCache(const Core::TimeWindow &tw,
                                   const std::string &networkCode,
                                   const std::string &stationCode,
                                   const std::string &locationCode,
                                   const std::string &channelCode,
                                   const GenericRecordCPtr &trace)
{
  const string wfId =
      waveformId(tw, networkCode, stationCode, locationCode, channelCode);
  _waveforms[wfId] = trace;
}

Core::TimeWindow
ExtraLenLoader::traceTimeWindowToLoad(const Core::TimeWindow &neededTW,
                                      const Core::Time &pickTime) const
{
  Core::TimeWindow twToLoad = neededTW;

  // Make sure to load at least _traceMinLenh seconds of waveform
  if (_traceMinLen > 0)
  {
    const Core::TimeSpan additionalTime(_traceMinLen / 2);

    if (twToLoad.startTime() > pickTime - additionalTime)
      twToLoad.setStartTime(pickTime - additionalTime);

    if (twToLoad.endTime() < pickTime + additionalTime)
      twToLoad.setEndTime(pickTime + additionalTime);
  }

  return twToLoad;
}

GenericRecordCPtr ExtraLenLoader::get(const Core::TimeWindow &tw,
                                      const Catalog::Phase &ph,
                                      const Catalog::Event &ev,
                                      bool demeaning,
                                      const std::string &filterStr,
                                      double resampleFreq)
{
  const Core::TimeWindow twToLoad = traceTimeWindowToLoad(tw, ph.time);
  GenericRecordCPtr trace =
      Loader::get(twToLoad, ph, ev, demeaning, filterStr, resampleFreq);
  if (trace && twToLoad != tw)
  {
    GenericRecordPtr nonConstTrace(new GenericRecord(*trace));
    if (!trim(*nonConstTrace, tw))
    {
      SEISCOMP_DEBUG("Incomplete trace, not enough data (%s)",
                     string(ph).c_str());
      return nullptr;
    }
    trace = nonConstTrace;
  }
  return trace;
}

Core::TimeWindow
SnrFilteredLoader::snrTimeWindow(const Core::Time &pickTime) const
{
  Core::Time winStart = std::min({pickTime + Core::TimeSpan(_snr.noiseStart),
                                  pickTime + Core::TimeSpan(_snr.signalStart)});
  Core::Time winEnd   = std::max({pickTime + Core::TimeSpan(_snr.noiseEnd),
                                pickTime + Core::TimeSpan(_snr.signalEnd)});
  return Core::TimeWindow(winStart, winEnd);
}

bool SnrFilteredLoader::goodSnr(const GenericRecordCPtr &trace,
                                const Core::Time &pickTime) const
{
  double snr = computeSnr(trace, pickTime, _snr.noiseStart, _snr.noiseEnd,
                          _snr.signalStart, _snr.signalEnd);
  return snr >= _snr.minSnr;
}

GenericRecordCPtr SnrFilteredLoader::get(const Core::TimeWindow &tw,
                                         const Catalog::Phase &ph,
                                         const Catalog::Event &ev,
                                         bool demeaning,
                                         const std::string &filterStr,
                                         double resampleFreq)
{
  const string wfId = waveformId(ph, tw);

  // Check if we have already excluded the trace SNR
  if (_snrExcludedWfs.count(wfId) != 0) { return nullptr; }

  const Core::TimeWindow twToLoad = tw | snrTimeWindow(ph.time);
  GenericRecordCPtr trace =
      Loader::get(twToLoad, ph, ev, demeaning, filterStr, resampleFreq);

  if (!trace) return nullptr;

  // Check if we have already validated the trace SNR, if not do it now
  if (_snrGoodWfs.count(wfId) == 0)
  {
    if (!goodSnr(trace, ph.time))
    {
      _snrExcludedWfs.insert(wfId);
      // Debugging waveforms: dump snr low traces
      if (!_wfDebugDir.empty())
      {
        SEISCOMP_DEBUG("Trace has too low SNR, discard it (%s)",
                       string(ph).c_str());
        writeTrace(trace,
                   waveformDebugPath(_wfDebugDir, ev, ph, "snr-rejected"));
      }
      _counters_wf_snr_low++;
      return nullptr;
    }
    _snrGoodWfs.insert(wfId);
  }
  if (twToLoad != tw)
  {
    GenericRecordPtr nonConstTrace(new GenericRecord(*trace));
    if (!trim(*nonConstTrace, tw))
    {
      SEISCOMP_DEBUG("Error when checking SNR, cannot trim data (%s)",
                     string(ph).c_str());
      return nullptr;
    }
    trace = nonConstTrace;
  }

  return trace;
}

} // namespace Waveform
} // namespace HDD
} // namespace Seiscomp

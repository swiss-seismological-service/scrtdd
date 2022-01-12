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

#include <fstream>
#include <iostream>
#include <mutex>

#include <seiscomp/client/inventory.h>
#include <seiscomp/core/datetime.h>
#include <seiscomp/core/genericrecord.h>
#include <seiscomp/core/recordsequence.h>
#include <seiscomp/core/timewindow.h>
#include <seiscomp/datamodel/station.h>
#include <seiscomp/io/recordinput.h>
#include <seiscomp/io/records/mseedrecord.h>
#include <seiscomp/io/recordstream.h>
#include <seiscomp/math/filter.h>
#include <seiscomp/math/geo.h>

using namespace std;
using std::chrono::duration;
using Catalog = HDD::Catalog;

using namespace Seiscomp;

namespace {

using Transform = HDD::Waveform::Processor::Transform;

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

Seiscomp::Core::Time toSC(const HDD::UTCTime &t)
{
  return Seiscomp::Core::Time(HDD::durToSec(t.time_since_epoch()));
}

HDD::UTCTime fromSC(const Seiscomp::Core::Time &t)
{
  return HDD::UTCTime() + HDD::secToDur(t.length());
}

Seiscomp::Core::TimeWindow toSC(const HDD::TimeWindow &tw)
{
  return Seiscomp::Core::TimeWindow(toSC(tw.startTime()), toSC(tw.endTime()));
}

HDD::TimeWindow fromSC(const Seiscomp::Core::TimeWindow &tw)
{
  return HDD::TimeWindow(fromSC(tw.startTime()), fromSC(tw.endTime()));
}

HDD::Trace fromSC(const Seiscomp::Record *rec)
{
  const DoubleArray *data = DoubleArray::ConstCast(rec->data());
  if (!data)
  {
    throw HDD::Exception("Internal logic error: cannot create HDD::Trace from "
                         "Seiscomp::Core::GenericRecord");
  }
  return HDD::Trace(rec->networkCode(), rec->stationCode(), rec->locationCode(),
                    rec->channelCode(), fromSC(rec->startTime()),
                    rec->samplingFrequency(), data->typedData(), data->size());
}

unique_ptr<HDD::Trace> contiguousRecord(const RecordSequence &seq,
                                        const HDD::TimeWindow &tw,
                                        double minAvailability)
{
  const Core::TimeWindow sctw = toSC(tw);
  if (seq.availability(sctw) < minAvailability)
  {
    return nullptr;
  }

  Seiscomp::GenericRecordPtr sctr = seq.contiguousRecord<double>(&sctw, false);
  unique_ptr<HDD::Trace> trace(new HDD::Trace(fromSC(sctr.get())));

  if (!trace->slice(tw))
  {
    return nullptr;
  }
  return trace;
}

unique_ptr<HDD::Trace>
readWaveformFromRecordStream(const string &recordStreamURL,
                             const HDD::TimeWindow &tw,
                             const string &networkCode,
                             const string &stationCode,
                             const string &locationCode,
                             const string &channelCode,
                             double tolerance,
                             double minAvailability,
                             bool inventoryCheck = false)
{
  const Core::TimeWindow sctw = toSC(tw);

  //
  // Checking the inventory is not stricly required since the RecordStream has
  // to perform a check on the data existence anyway. However the check can
  // save time by skipping the RecorStream request, so we do it.
  //
  Seiscomp::DataModel::Stream *stream = nullptr;
  if (inventoryCheck)
  {
    Seiscomp::DataModel::Inventory *inv =
        Seiscomp::Client::Inventory::Instance()->inventory();

    if (inv)
    {
      Seiscomp::DataModel::InventoryError error;
      stream = Seiscomp::DataModel::getStream(inv, networkCode, stationCode,
                                              locationCode, channelCode,
                                              sctw.startTime(), &error);
      // We return without even trying to load the data only if the stream was
      // found in the inventory, but it was not available at the requested point
      // in time. If the stream was not found the inventory might be incomplete
      // and we don't want to stop for that, so we proceed with the request
      if (!stream && (error == DataModel::STREAM_EPOCH_NOT_FOUND ||
                      error == DataModel::SENSOR_EPOCH_NOT_FOUND ||
                      error == DataModel::STATION_EPOCH_NOT_FOUND ||
                      error == DataModel::NETWORK_EPOCH_NOT_FOUND))
      {
        string msg = HDD::strf(
            "Cannot find '%s.%s.%s.%s' at time '%s' in the inventory: %s",
            networkCode.c_str(), stationCode.c_str(), locationCode.c_str(),
            channelCode.c_str(),
            HDD::UTCClock::toString(tw.startTime()).c_str(), error.toString());
        throw HDD::Exception(msg);
      }
    }
  }

  IO::RecordStreamPtr rs = IO::RecordStream::Open(recordStreamURL.c_str());
  if (rs == nullptr)
  {
    string msg = "Cannot open RecordStream: " + recordStreamURL;
    throw HDD::Exception(msg);
  }

  rs->setTimeWindow(sctw);
  rs->addStream(networkCode, stationCode, locationCode, channelCode);

  IO::RecordInput inp(rs.get(), Array::DOUBLE, Record::DATA_ONLY);
  TimeWindowBuffer seq(sctw, tolerance);
  RecordPtr rec;
  while (rec = inp.next())
  {
    seq.feed(rec.get());
  }
  rs->close();

  unique_ptr<HDD::Trace> trace = contiguousRecord(seq, tw, minAvailability);

  if (!trace)
  {
    string msg =
        HDD::strf("Data availability too low %.2f%%", seq.availability());
    throw HDD::Exception(msg);
  }

  // Apply gain (useful for ZRT and L2 transform)
  if (stream != nullptr && stream->gain() != 0)
  {
    double scaler = 1.0 / stream->gain();
    double *data  = trace->data();
    for (size_t i = 0; i < trace->sampleCount(); ++i) *data *= scaler;
  }

  return trace;
}

struct ThreeComponents
{
  enum Component
  {
    Vertical         = 0, /* usually Z */
    FirstHorizontal  = 1, /* usually N */
    SecondHorizontal = 2  /* usually E */
  };
  string names[3];
  double dip[3];
  double azimuth[3];
};

bool getComponentsInfo(const Catalog::Phase &ph,
                       ThreeComponents &components,
                       Transform trans)
{
  string baseErrMsg =
      HDD::strf("Unable to fetch components information for phase '%s'",
                string(ph).c_str());

  const Core::Time sctime = toSC(ph.time);
  const string channelCodeRoot =
      HDD::Waveform::getBandAndInstrumentCodes(ph.channelCode);

  Seiscomp::DataModel::Inventory *inv =
      Seiscomp::Client::Inventory::Instance()->inventory();
  if (!inv)
  {
    HDD::logDebug("%s: inventory not available", baseErrMsg.c_str());
    return false;
  }

  Seiscomp::DataModel::InventoryError error;
  Seiscomp::DataModel::SensorLocation *loc =
      Seiscomp::DataModel::getSensorLocation(
          inv, ph.networkCode, ph.stationCode, ph.locationCode, sctime, &error);

  if (!loc)
  {
    HDD::logDebug(
        "%s: unable to fetch SensorLocation information from inventory (%s)",
        baseErrMsg.c_str(), error.toString());
    return false;
  }

  DataModel::ThreeComponents tc;
  getThreeComponents(tc, loc, channelCodeRoot.c_str(), sctime);

  if (trans == Transform::TRANSVERSAL || trans == Transform::RADIAL)
  {
    if (tc.comps[ThreeComponents::Vertical] &&
        tc.comps[ThreeComponents::FirstHorizontal] &&
        tc.comps[ThreeComponents::SecondHorizontal])
    {
      components.names[ThreeComponents::Vertical] =
          tc.comps[ThreeComponents::Vertical]->code();
      components.names[ThreeComponents::FirstHorizontal] =
          tc.comps[ThreeComponents::FirstHorizontal]->code();
      components.names[ThreeComponents::SecondHorizontal] =
          tc.comps[ThreeComponents::SecondHorizontal]->code();

      components.dip[ThreeComponents::Vertical] =
          tc.comps[ThreeComponents::Vertical]->dip();
      components.dip[ThreeComponents::FirstHorizontal] =
          tc.comps[ThreeComponents::FirstHorizontal]->dip();
      components.dip[ThreeComponents::SecondHorizontal] =
          tc.comps[ThreeComponents::SecondHorizontal]->dip();

      components.azimuth[ThreeComponents::Vertical] =
          tc.comps[ThreeComponents::Vertical]->azimuth();
      components.azimuth[ThreeComponents::FirstHorizontal] =
          tc.comps[ThreeComponents::FirstHorizontal]->azimuth();
      components.azimuth[ThreeComponents::SecondHorizontal] =
          tc.comps[ThreeComponents::SecondHorizontal]->azimuth();
      return true;
    }
    HDD::logDebug("%s: three components information not found in inventory",
                  baseErrMsg.c_str());
    return false;
  }
  else if (trans == Transform::L2)
  {
    if (tc.comps[ThreeComponents::FirstHorizontal] &&
        tc.comps[ThreeComponents::SecondHorizontal])
    {
      components.names[ThreeComponents::FirstHorizontal] =
          tc.comps[ThreeComponents::FirstHorizontal]->code();
      components.names[ThreeComponents::SecondHorizontal] =
          tc.comps[ThreeComponents::SecondHorizontal]->code();
      return true;
    }
    HDD::logDebug(
        "%s: horizontal components information not found in inventory",
        baseErrMsg.c_str());
    return false;
  }

  HDD::logDebug("%s: unknown transformation", baseErrMsg.c_str());
  return false;
}

unique_ptr<HDD::Trace> transformL2(const HDD::TimeWindow &tw,
                                   const Catalog::Phase &ph,
                                   const ThreeComponents &tc,
                                   const HDD::Trace trH1,
                                   const HDD::Trace trH2)
{
  string channelCodeRoot =
      HDD::Waveform::getBandAndInstrumentCodes(ph.channelCode);

  auto tracesCompatible = [](const HDD::Trace &tr1,
                             const HDD::Trace &tr2) -> bool {
    return tr1.startTime() == tr2.startTime() &&
           tr1.samplingFrequency() == tr2.samplingFrequency() &&
           tr1.sampleCount() == tr2.sampleCount();
  };

  if (!tracesCompatible(trH1, trH2))
  {
    HDD::logDebug(
        "Cannot perform L2 transformation with incompatible horizontal traces");
    return nullptr;
  }

  vector<double> l2vec(trH1.sampleCount());
  double *l2data       = l2vec.data();
  const double *h1data = trH1.data();
  const double *h2data = trH2.data();
  for (size_t i = 0; i < l2vec.size(); ++i)
  {
    l2data[i] = std::sqrt(h1data[i] * h1data[i] + h2data[i] * h2data[i]);
  }
  unique_ptr<HDD::Trace> l2trace(new HDD::Trace(
      ph.networkCode, ph.stationCode, ph.locationCode, channelCodeRoot + "L",
      trH1.startTime(), trH1.samplingFrequency(), std::move(l2vec)));
  return l2trace;
}

unique_ptr<HDD::Trace> transformRT(const HDD::TimeWindow &tw,
                                   const Catalog::Phase &ph,
                                   const Catalog::Event &ev,
                                   const Catalog::Station &sta,
                                   const ThreeComponents &tc,
                                   const HDD::Trace trV,
                                   const HDD::Trace trH1,
                                   const HDD::Trace trH2,
                                   Transform trans)
{
  string channelCodeRoot =
      HDD::Waveform::getBandAndInstrumentCodes(ph.channelCode);

  auto tracesCompatible = [](const HDD::Trace &tr1,
                             const HDD::Trace &tr2) -> bool {
    return tr1.startTime() == tr2.startTime() &&
           tr1.samplingFrequency() == tr2.samplingFrequency() &&
           tr1.sampleCount() == tr2.sampleCount();
  };

  if (!tracesCompatible(trH1, trH2) || !tracesCompatible(trH2, trV))
  {
    HDD::logDebug("Cannot perform L2 transformation with incompatible traces");
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
      .fromAngles(+HDD::degToRad(tc.azimuth[ThreeComponents::SecondHorizontal]),
                  -HDD::degToRad(tc.dip[ThreeComponents::SecondHorizontal]))
      .normalize();
  Vector3d h1Col;
  h1Col
      .fromAngles(+HDD::degToRad(tc.azimuth[ThreeComponents::FirstHorizontal]),
                  -HDD::degToRad(tc.dip[ThreeComponents::FirstHorizontal]))
      .normalize();
  Vector3d vCol;
  vCol.fromAngles(+HDD::degToRad(tc.azimuth[ThreeComponents::Vertical]),
                  -HDD::degToRad(tc.dip[ThreeComponents::Vertical]))
      .normalize();

  Matrix3d orientationZNE = {
      h2Col.x, h1Col.x, vCol.x,  h2Col.y, h1Col.y,
      vCol.y,  h2Col.z, h1Col.z, vCol.z,
  };

  //
  // Rotation matrix ZRT
  //
  double backAzimuth;
  HDD::computeDistance(ev.latitude, ev.longitude, sta.latitude, sta.longitude,
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
  // Apply transformation now
  //
  vector<double> tvec(trH1.sampleCount());
  double *tdata        = tvec.data();
  const double *vdata  = trV.data();
  const double *h1data = trH1.data();
  const double *h2data = trH2.data();
  for (size_t i = 0; i < tvec.size(); ++i)
  {
    if (trans == Transform::TRANSVERSAL) // T trasversal
      tdata[i] = transformation[0][0] * h2data[i] +
                 transformation[0][1] * h1data[i] +
                 transformation[0][2] * vdata[i];
    else if (trans == Transform::RADIAL) // R radial
      tdata[i] = transformation[1][0] * h2data[i] +
                 transformation[1][1] * h1data[i] +
                 transformation[1][2] * vdata[i];
    else // V vertical
      tdata[i] = transformation[2][0] * h2data[i] +
                 transformation[2][1] * h1data[i] +
                 transformation[2][2] * vdata[i];
  }

  // build and return the trace
  unique_ptr<HDD::Trace> ttrace(new HDD::Trace(
      ph.networkCode, ph.stationCode, ph.locationCode, channelCodeRoot + "L",
      trH1.startTime(), trH1.samplingFrequency(), std::move(tvec)));
  return ttrace;
}

} // namespace

namespace HDD {
namespace Waveform {

string waveformId(const TimeWindow &tw,
                  const string &networkCode,
                  const string &stationCode,
                  const string &locationCode,
                  const string &channelCode)
{
  return strf("%s.%s.%s.%s.%s.%s", networkCode.c_str(), stationCode.c_str(),
              locationCode.c_str(), channelCode.c_str(),
              UTCClock::toString(tw.startTime()).c_str(),
              UTCClock::toString(tw.endTime()).c_str());
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

void filter(Trace &trace,
            bool demeaning,
            const string &filterStr,
            double resampleFreq)
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
    string filterError;
    auto filter =
        Math::Filtering::InPlaceFilter<double>::Create(filterStr, &filterError);
    if (!filter)
    {
      string msg = strf("Filter creation failed %s: %s", filterStr.c_str(),
                        filterError.c_str());
      throw Exception(msg);
    }
    filter->setSamplingFrequency(trace.samplingFrequency());
    filter->apply(data_size, data);
    delete filter;
  }

  if (resampleFreq > 0)
  {
    resample(trace, resampleFreq);
  }
}

double computeSnr(const Trace &tr,
                  const UTCTime &pickTime,
                  double noiseOffsetStart,
                  double noiseOffsetEnd,
                  double signalOffsetStart,
                  double signalOffsetEnd)
{
  const double *data          = tr.data();
  const size_t data_size      = tr.sampleCount();
  const double freq           = tr.samplingFrequency();
  const UTCTime dataStartTime = tr.startTime();

  // convert time w.r.t. guiding pick time to sample number
  auto secToSample = [&freq, &data_size](double sec) {
    double max = std::max(std::round(sec * freq), 0.);
    return std::min(max, data_size - 1.);
  };
  const double pickOffset  = durToSec(pickTime - dataStartTime);
  const size_t noiseStart  = secToSample(noiseOffsetStart + pickOffset);
  const size_t noiseEnd    = secToSample(noiseOffsetEnd + pickOffset);
  const size_t signalStart = secToSample(signalOffsetStart + pickOffset);
  const size_t signalEnd   = secToSample(signalOffsetEnd + pickOffset);

  if ((std::max({noiseStart, noiseEnd, signalStart, signalEnd}) >= data_size))
  {
    logError(
        "Cannot compute SNR: noise/signal windows exceed waveform boundaries");
    return -1;
  }

  // get maximum (absolute) amplitude in noise window
  double noiseMax = -1.0;
  for (size_t i = noiseStart; i < noiseEnd; i++)
  {
    noiseMax = std::max(std::abs(data[i]), noiseMax);
  }

  // get maximum (absolute) amplitude in signal window
  double signalMax = -1.0;
  for (size_t i = signalStart; i < signalEnd; i++)
  {
    signalMax = std::max(std::abs(data[i]), signalMax);
  }

  return signalMax / noiseMax;
}

void writeTrace(const Trace &trace, const string &file)
{
  try
  {
    ofstream ofs(file);

    GenericRecord gr(trace.networkCode(), trace.stationCode(),
                     trace.locationCode(), trace.channelCode(),
                     toSC(trace.startTime()), trace.samplingFrequency());
    gr.setData(trace.sampleCount(), trace.data(), Array::DOUBLE);

    IO::MSeedRecord msRec(gr);
    int reclen = msRec.data()->size() * msRec.data()->elementSize() + 64;
    reclen =
        nextPowerOf2<int>(reclen, 128 /*MINRECLEN*/, 1048576 /*MAXRECLEN*/);
    if (reclen > 0)
    {
      msRec.setOutputRecordLength(reclen);
      msRec.write(ofs);
    }
  }
  catch (exception &e)
  {
    logWarning("Couldn't write waveform to disk %s: %s", file.c_str(),
               e.what());
  }
}

unique_ptr<Trace> readTrace(const string &file)
{
  if (!pathExists(file)) return nullptr;

  try
  {
    ifstream ifs(file);
    IO::MSeedRecord msRec(Array::DOUBLE, Record::Hint::DATA_ONLY);
    msRec.read(ifs);
    return unique_ptr<Trace>(new Trace(fromSC(&msRec)));
  }
  catch (exception &e)
  {
    logWarning("Couldn't load waveform %s: %s", file.c_str(), e.what());
    return nullptr;
  }
}

double Loader::_tolerance{0.1};
double Loader::_minAvailability{0.95};

shared_ptr<const Trace> BasicLoader::get(const TimeWindow &tw,
                                         const Catalog::Phase &ph)
{
  unique_ptr<Trace> trace;

  if (!_recordStreamURL.empty())
  {
    try
    {
      trace = readWaveformFromRecordStream(
          _recordStreamURL, tw, ph.networkCode, ph.stationCode, ph.locationCode,
          ph.channelCode, _tolerance, _minAvailability);
    }
    catch (exception &e)
    {
      logDebug(
          "Cannnot load trace (stream %s.%s.%s.%s from %s length %.2f sec): %s",
          ph.networkCode.c_str(), ph.stationCode.c_str(),
          ph.locationCode.c_str(), ph.channelCode.c_str(),
          HDD::UTCClock::toString(tw.startTime()).c_str(),
          durToSec(tw.length()), e.what());
    }
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
    const string wfId = waveformId(tw, ph.networkCode, ph.stationCode,
                                   ph.locationCode, ph.channelCode);
    const auto it     = _waveforms.find(wfId);
    return it != _waveforms.end() ? it->second : nullptr;
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
  auto eqlrng = _request.equal_range(streamID);
  for (auto it = eqlrng.first; it != eqlrng.second; ++it)
  {
    if (it->second == tw) return; // skip duplicated requests
  }
  _request.emplace(streamID, tw);
}

void BatchLoader::load()
{
  if (_dataLoaded) return;

  unordered_multimap<string, const Core::TimeWindow> reqCopy;
  unordered_multimap<string, TimeWindowBuffer> streamBuf;
  for (const auto &kv : _request)
  {
    Core::TimeWindow tw_ = toSC(kv.second);
    reqCopy.emplace(kv.first, tw_);
    streamBuf.emplace(kv.first, TimeWindowBuffer(tw_, _tolerance));
  }

  IO::RecordStreamPtr rs;

  //
  // The reason of this loop is that RecordStream can only handle one single
  // time window per stream
  //
  while (!reqCopy.empty())
  {
    rs = IO::RecordStream::Open(_recordStreamURL.c_str());
    if (rs == nullptr)
    {
      logError("Cannot open RecordStream: %s", _recordStreamURL.c_str());
      break;
    }

    //
    // Convert multiple time windows requests close to each others to a
    // single RecordStream request on the stream
    //
    for (auto it = reqCopy.begin(), end = reqCopy.end();
         it != end;) // loop by stream
    {
      const string streamID = it->first;
      Core::TimeWindow contiguousRequest;
      auto eqlrng = reqCopy.equal_range(streamID);
      for (auto it2 = eqlrng.first;
           it2 != eqlrng.second;) // loop by stream windows
      {
        const Core::TimeWindow &tw = it2->second;
        bool requested             = false;

        if (it2 == eqlrng.first)
        {
          contiguousRequest = tw;
          requested         = true;
        }
        else if (contiguousRequest.overlaps(tw))
        {
          contiguousRequest = contiguousRequest.merge(tw);
          requested         = true;
        }
        else if (contiguousRequest.endTime() <= tw.startTime() &&
                 contiguousRequest.contiguous(tw, 60))
        {
          contiguousRequest = contiguousRequest.merge(tw);
          requested         = true;
        }
        else if (contiguousRequest.startTime() >= tw.endTime() &&
                 tw.contiguous(contiguousRequest, 60))
        {
          contiguousRequest = contiguousRequest.merge(tw);
          requested         = true;
        }

        if (requested)
          it2 = reqCopy.erase(it2);
        else
          it2++;
      }
      static const std::regex dot("\\.", std::regex::optimize);
      const vector<string> tokens(splitString(streamID, dot));
      rs->addStream(tokens.at(0), tokens.at(1), tokens.at(2), tokens.at(3),
                    contiguousRequest.startTime(), contiguousRequest.endTime());
      it = eqlrng.second;
    }

    //
    // Collect data records and feed them to the corresponding
    // TimeWindowBuffer
    //
    IO::RecordInput inp(rs.get(), Array::DOUBLE, Record::DATA_ONLY);
    RecordPtr rec;
    while (rec = inp.next())
    {
      auto eqlrng = streamBuf.equal_range(rec->streamID());
      for (auto it = eqlrng.first; it != eqlrng.second; ++it)
      {
        TimeWindowBuffer &seq = it->second;
        seq.feed(rec.get());
      }
    }
    rs->close();
  }

  //
  // Convert data records to contiguous waveforms
  //
  for (auto &kv : streamBuf)
  {
    const string &streamID = kv.first;
    TimeWindowBuffer &seq  = kv.second;
    const TimeWindow &tw   = fromSC(seq.timeWindowToStore());

    shared_ptr<const Trace> trace = contiguousRecord(seq, tw, _minAvailability);
    if (!trace)
    {
      logDebug("Cannnot load trace, data availability %.2f%%"
               "(stream %s from %s length %.2f sec)",
               seq.availability(), streamID.c_str(),
               UTCClock::toString(tw.startTime()).c_str(), tw.length());
      _counters_wf_no_avail++;
      continue;
    }
    const string wfId =
        waveformId(tw, trace->networkCode(), trace->stationCode(),
                   trace->locationCode(), trace->channelCode());
    _waveforms[wfId] = trace;
    _counters_wf_downloaded++;
  }

  logInfo("Fetched %u/%lu waveforms, not available %u", _counters_wf_downloaded,
          _request.size(), _counters_wf_no_avail);

  _request.clear();
  _dataLoaded = true;
}

TimeWindow ExtraLenLoader::traceTimeWindowToLoad(const TimeWindow &neededTW,
                                                 const UTCTime &pickTime) const
{
  TimeWindow twToLoad = neededTW;

  // Make sure to load at least `_traceMinLenh` seconds of the waveform.
  if (_beforePickLen > 0 || _afterPickLen > 0)
  {
    const UTCTime::duration additionalTimeBefore = secToDur(_beforePickLen);
    const UTCTime::duration additionalTimeAfter  = secToDur(_afterPickLen);

    if (twToLoad.startTime() > pickTime - additionalTimeBefore)
      twToLoad.setStartTime(pickTime - additionalTimeBefore);

    if (twToLoad.endTime() < pickTime + additionalTimeAfter)
      twToLoad.setEndTime(pickTime + additionalTimeAfter);
  }

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
      logDebug("Incomplete trace, not enough data (%s)", string(ph).c_str());
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
  return readTrace(cacheFile);
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
  writeTrace(trace, cacheFile);
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
  auto loadWaveform = [this, &tw, &ph, &filterStr,
                       resampleFreq](const string &channelCode) {
    Catalog::Phase copy(ph);
    copy.channelCode              = channelCode;
    shared_ptr<const Trace> trace = _auxLdr->get(tw, copy);
    if (trace) trace = process(*trace, filterStr, resampleFreq);
    return trace;
  };

  shared_ptr<const Trace> trace;

  if (trans == Transform::NONE)
  {
    trace = loadWaveform(ph.channelCode);
  }
  else
  {
    ThreeComponents comps;
    if (getComponentsInfo(ph, comps, trans))
    {
      if (trans == Transform::TRANSVERSAL || trans == Transform::RADIAL)
      {
        shared_ptr<const Trace> trV(
            loadWaveform(comps.names[ThreeComponents::Vertical]));
        shared_ptr<const Trace> trH1(
            loadWaveform(comps.names[ThreeComponents::FirstHorizontal]));
        shared_ptr<const Trace> trH2(
            loadWaveform(comps.names[ThreeComponents::SecondHorizontal]));
        if (trH1 && trH2 && trV)
          trace =
              transformRT(tw, ph, ev, sta, comps, *trV, *trH1, *trH2, trans);
      }
      else if (trans == Transform::L2)
      {
        shared_ptr<const Trace> trH1(
            loadWaveform(comps.names[ThreeComponents::FirstHorizontal]));
        shared_ptr<const Trace> trH2(
            loadWaveform(comps.names[ThreeComponents::SecondHorizontal]));
        if (trH1 && trH2) trace = transformL2(tw, ph, comps, *trH1, *trH2);
      }
    }
    if (!trace)
    {
      logDebug("Cannnot compute RT/L2 transformation of stream %s.%s.%s.%s "
               "from %s length %.2f sec",
               ph.networkCode.c_str(), ph.stationCode.c_str(),
               ph.locationCode.c_str(),
               getBandAndInstrumentCodes(ph.channelCode).c_str(),
               HDD::UTCClock::toString(tw.startTime()).c_str(),
               durToSec(tw.length()));
    }
  }

  return trace;
}

std::shared_ptr<Trace> BasicProcessor::process(const Trace &trace,
                                               const std::string &filterStr,
                                               double resampleFreq) const
{
  shared_ptr<Trace> copy(new Trace(trace));
  try
  {
    filter(*copy, true, filterStr, resampleFreq);
  }
  catch (exception &e)
  {
    logWarning("Errow while filtering waveform: %s", e.what());
    copy.reset();
  }
  return copy;
}

shared_ptr<const Trace> MemCachedProc::get(const TimeWindow &tw,
                                           const Catalog::Phase &ph,
                                           const Catalog::Event &ev,
                                           const Catalog::Station &sta,
                                           const std::string &filterStr,
                                           double resampleFreq,
                                           Transform trans)
{
  const string wfId             = waveformId(tw, ph.networkCode, ph.stationCode,
                                 ph.locationCode, ph.channelCode);
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
  const string wfId = waveformId(ph, tw);
  const auto it     = _waveforms.find(wfId);
  return it != _waveforms.end();
}

shared_ptr<const Trace> MemCachedProc::getFromCache(const string &wfId)
{
  const auto it = _waveforms.find(wfId);
  return it != _waveforms.end() ? it->second : nullptr;
}

void MemCachedProc::storeInCache(const string &wfId,
                                 const shared_ptr<const Trace> &trace)
{
  _waveforms[wfId] = trace;
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
  const string wfId = waveformId(ph, tw);

  const TimeWindow twToLoad = snrTimeWindow(ph.time).merge(tw);
  shared_ptr<const Trace> trace =
      _auxPrc->get(twToLoad, ph, ev, sta, filterStr, resampleFreq, trans);

  if (!trace) return nullptr;

  if (!goodSnr(*trace, ph.time))
  {
    logDebug("Trace has too low SNR(%s)", string(ph).c_str());
    _counters_wf_snr_low++;
    return nullptr;
  }

  if (twToLoad != tw)
  {
    shared_ptr<Trace> nonConstTrace(new Trace(*trace));
    if (!nonConstTrace->slice(tw))
    {
      logDebug("Error when checking SNR, cannot trim data (%s)",
               string(ph).c_str());
      return nullptr;
    }
    trace = nonConstTrace;
  }

  return trace;
}

} // namespace Waveform
} // namespace HDD

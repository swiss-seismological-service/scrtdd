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
#include "utils.h"

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include <cfenv>
#include <fstream>
#include <iostream>
#include <mutex>
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
using namespace Seiscomp::Processing;
using DataModel::ThreeComponents;
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
  string debugFile = HDD::strf(
      "ev%u.%s.%s.%s.%s.%s.%s.mseed", ev.id, ph.networkCode.c_str(),
      ph.stationCode.c_str(), ph.locationCode.c_str(), ph.channelCode.c_str(),
      ph.type.c_str(), ext.c_str());
  return (boost::filesystem::path(wfDebugDir) / debugFile).string();
}

} // namespace

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
  return strf("%s.%s.%s.%s.%s.%s", networkCode.c_str(), stationCode.c_str(),
              locationCode.c_str(), channelCode.c_str(),
              tw.startTime().iso().c_str(), tw.endTime().iso().c_str());
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
                                              const string &channelCode,
                                              double tolerance,
                                              double minAvailability)
{
  IO::RecordStreamPtr rs = IO::RecordStream::Open(recordStreamURL.c_str());
  if (rs == nullptr)
  {
    string msg = "Cannot open RecordStream: " + recordStreamURL;
    throw Exception(msg);
  }

  rs->setTimeWindow(tw);
  rs->addStream(networkCode, stationCode, locationCode, channelCode);

  IO::RecordInput inp(rs.get(), Array::DOUBLE, Record::DATA_ONLY);
  TimeWindowBuffer seq(tw, tolerance);
  RecordPtr rec;
  while (rec = inp.next())
  {
    seq.feed(rec.get());
  }
  rs->close();

  GenericRecordPtr trace =
      contiguousRecord(seq, tw, tolerance, minAvailability);

  if (!trace)
  {
    string msg =
        strf("Cannnot load trace, data availability %.2f%%"
             "(stream %s.%s.%s.%s from %s length %.2f sec)",
             seq.availability(), networkCode.c_str(), stationCode.c_str(),
             locationCode.c_str(), channelCode.c_str(),
             tw.startTime().iso().c_str(), tw.length());
    throw Exception(msg);
  }

  return trace;
}

GenericRecordPtr contiguousRecord(const RecordSequence &seq,
                                  const Core::TimeWindow &tw,
                                  double tolerance,
                                  double minAvailability)
{
  if (seq.availability(tw) < minAvailability)
  {
    return nullptr;
  }

  GenericRecordPtr trace(seq.contiguousRecord<double>(&tw, false));

  if (!trim(*trace, tw))
  {
    return nullptr;
  }

  return trace;
}

bool trim(GenericRecord &trace, const Core::TimeWindow &tw)
{
  if (trace.timeWindow() == tw) return true;

  int startOfs = std::floor(double(tw.startTime() - trace.startTime()) *
                            trace.samplingFrequency());
  int endOfs   = std::ceil(double(tw.endTime() - trace.startTime()) *
                         trace.samplingFrequency());

  // one sample tolerance
  if (startOfs == -1) startOfs++;
  if (endOfs == trace.data()->size() + 1) endOfs--;

  // not enough data at start of time window
  if (startOfs < 0) return false;

  // not enough data at end of time window
  if (endOfs > trace.data()->size()) return false;

  ArrayPtr sliced = trace.data()->slice(startOfs, endOfs);
  trace.setStartTime(trace.startTime() +
                     Core::TimeSpan(startOfs / trace.samplingFrequency()));
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
    filter->apply(data->size(), data->typedData());
    delete filter;
    trace.dataUpdated();
  }

  if (resampleFreq > 0)
  {
    resample(trace, resampleFreq);
  }
}

void resample(GenericRecord &trace, double new_sf)
{
  if (new_sf <= 0 || trace.samplingFrequency() == new_sf ||
      !DoubleArray::Cast(trace.data()))
    return;

  const double data_sf       = trace.samplingFrequency();
  const double resamp_factor = new_sf / data_sf;
  const double nyquist       = std::min(new_sf, data_sf) / 2.;

  const double *data  = DoubleArray::Cast(trace.data())->typedData();
  const int data_size = trace.data()->size();

  const int resampled_data_size = data_size * resamp_factor;
  double resampled_data[resampled_data_size];

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
      if (j >= 0 && j < data_size)
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

  for (int i = 0; i < resampled_data_size; i++)
  {
    double x          = i / resamp_factor;
    resampled_data[i] = new_sample(x, nyquist, 21);
  }

  DoubleArray::Cast(trace.data())->setData(resampled_data_size, resampled_data);
  trace.setSamplingFrequency(new_sf);
  trace.dataUpdated();
}

void writeTrace(GenericRecordCPtr trace, const std::string &file)
{
  if (!trace) return;

  try
  {
    std::ofstream ofs(file);
    IO::MSeedRecord msRec(*trace);
    int reclen = msRec.data()->size() * msRec.data()->elementSize() + 64;
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
    trace->setData(msRec.data()->clone()); // copy data, too
    return trace;
  }
  catch (exception &e)
  {
    SEISCOMP_WARNING("Couldn't load waveform %s: %s", file.c_str(), e.what());
    return nullptr;
  }
}

/*
 * Compute cross-correlation between two traces centered around their respective
 * picks. The cross-correlation will be performed from the longest trace middle
 * minus 'maxDelay' to the same trace middle plus 'maxDelay' (if enough data is
 * available).
 * `delayOut` will store the shift in seconds (positive or negative) from the
 * longest trace middle point at which there is the highest (absolute value)
 * correlation coefficient, stored in 'coeffOut'
 */
bool xcorr(const GenericRecordCPtr &tr1,
           const GenericRecordCPtr &tr2,
           double maxDelay,
           bool qualityCheck,
           double &delayOut,
           double &coeffOut)
{
  if (tr1->samplingFrequency() != tr2->samplingFrequency())
  {
    SEISCOMP_INFO(
        "Cannot cross correlate traces with different sampling freq (%f!=%f)",
        tr1->samplingFrequency(), tr2->samplingFrequency());
    return false;
  }
  const double freq = tr1->samplingFrequency();

  // check longest/shortest trace
  const bool swap             = tr1->data()->size() > tr2->data()->size();
  GenericRecordCPtr trShorter = swap ? tr2 : tr1;
  GenericRecordCPtr trLonger  = swap ? tr1 : tr2;

  const double *dataS = DoubleArray::ConstCast(trShorter->data())->typedData();
  const double *dataL = DoubleArray::ConstCast(trLonger->data())->typedData();
  const int sizeS     = trShorter->data()->size();
  const int sizeL     = trLonger->data()->size();

  // force to cross-correlate withing data boundaries
  int availableData = (sizeL - sizeS) / 2;
  int maxDelaySmps  = maxDelay * freq;
  if (maxDelaySmps > availableData) maxDelaySmps = availableData;

  crossCorrelation(dataS, sizeS, (dataL + availableData - maxDelaySmps),
                   (sizeS + maxDelaySmps * 2), qualityCheck, delayOut,
                   coeffOut);

  if (!std::isfinite(coeffOut))
  {
    coeffOut = 0;
    delayOut = 0.;
  }
  else
  {
    delayOut -= maxDelaySmps; // the reference is the middle of the long trace
    delayOut /= freq;         // samples to secs
    if (swap) delayOut = -delayOut;
  }
  return true;
}

void crossCorrelation(const double *dataS,
                      const int sizeS,
                      const double *dataL,
                      const int sizeL,
                      bool qualityCheck,
                      double &delayOut,
                      double &coeffOut)
{
  /*
   * for later quality check: save local maxima/minima
   */
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

  /*
   * Pearson correlation coefficient for time series X and Y of length n
   *
   *              sum((Xi-meanX) * (Yi-meanY))
   * cc = --------------------------------------------------
   *      sqrt(sum((Xi-meanX)^2)) * sqrt(sum((Yi-meanY)^2))
   *
   * Where sum(X)  is the sum of Xi for i=1 until i=n
   *
   * This can be rearranged in a form suitable for a single-pass algorithm
   * (where the mean of X and Y are not needed)
   *
   *                 n * sum(Xi*Yi) - sum(Xi) * sum(Yi)
   * cc = -----------------------------------------------------------
   *      sqrt(n*sum(Xi^2)-sum(Xi)^2) * sqrt(n*sum(Yi^2)-sum(Yi)^2))
   *
   * For cross-correlation, where we have a short trace S which is correlated
   * against a longer trace L at subsequent offset, we can pre-compute the
   * parts that involves S and re-use them at each step of the
   * cross-correlation:
   *
   *   sumS   = sum(Xi)
   *   denomS = sqrt(n*sum(Xi^2)-sum(Xi)^2)
   *
   * For the parts that involves the longer trace L alone we can compute them
   * in a rolling fashion (removing first sample of previous iteration and
   * adding the last sample of the new iteration):
   *
   *   sumL   = sum(Yi)
   *   sumL2  = sum(Yi^2)
   *   denomL = sqrt(n*sumL2-sumL^2))
   *
   * Finally, this is the equation at each step (offset) of cross-correlation:
   *
   *       n * sum(Xi*Yi) - sumS * sumL
   * cc = ------------------------------
   *             denomS * denomL
   *
   * Unfortunately, we cannot optimize sum(Xi*Yi) and this will be a inner
   * loop inside the main cross-correlation loop.
   */

  std::feclearexcept(FE_ALL_EXCEPT);

  // prepare the data before the main xcorr loop
  const int n = sizeS;
  double sumS = 0, sumS2 = 0;
  double sumL = 0, sumL2 = 0;
  for (int i = 0; i < n; i++)
  {
    sumS += dataS[i];
    sumS2 += dataS[i] * dataS[i];
    if (i >= (n - 1)) continue;
    sumL += dataL[i];
    sumL2 += dataL[i] * dataL[i];
  }
  double denomS = std::sqrt(n * sumS2 - sumS * sumS);

  // cross-correlation loop
  coeffOut           = std::nan("");
  double lastSampleL = 0;
  for (int delay = 0; delay <= (sizeL - sizeS); delay++)
  {
    // sumL/sumL2 update: remove the sample that has just exited the
    // current cross-correlation win and add the sample that has just
    // entered
    const double newSampleL = dataL[delay + n - 1];
    sumL += newSampleL - lastSampleL;
    sumL2 += newSampleL * newSampleL - lastSampleL * lastSampleL;
    lastSampleL = dataL[delay]; // prepare for next loop

    const double denomL = std::sqrt(n * sumL2 - sumL * sumL);

    double sumSL = 0;
    for (int i = 0; i < n; i++) sumSL += dataS[i] * dataL[i + delay];

    const double coeff = (n * sumSL - sumS * sumL) / (denomS * denomL);

    if (!std::isfinite(coeffOut) || std::abs(coeff) > std::abs(coeffOut))
    {
      coeffOut = coeff;
      delayOut = delay;
    }

    // for later quality check
    localMaxs.update(coeff);
    localMins.update(-coeff);
  }

  int fe = fetestexcept(FE_ALL_EXCEPT);
  if ((fe & ~FE_INEXACT) != 0) // we don't care about FE_INEXACT
  {
    SEISCOMP_WARNING("Floating point exception during cross-correlation:");
    if (fe & FE_DIVBYZERO) SEISCOMP_WARNING("FE_DIVBYZERO");
    if (fe & FE_INVALID) SEISCOMP_WARNING("FE_INVALID");
    if (fe & FE_OVERFLOW) SEISCOMP_WARNING("FE_OVERFLOW");
    if (fe & FE_UNDERFLOW) SEISCOMP_WARNING("FE_UNDERFLOW");
  }

  /*
   * To avoid errors introduced by cycle skipping, the differential time
   * measurement is only accepted if all side lobe maxima CCslm of the
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
  auto secToSample = [&freq, &size](double sec) {
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

  // get maximum (absolute) amplitude in noise window
  double noiseMax = -1.0;
  for (int i = noiseStart; i < noiseEnd; i++)
  {
    noiseMax = std::max(std::abs(data[i]), noiseMax);
  }

  // get maximum (absolute) amplitude in signal window
  double signalMax = -1.0;
  for (int i = signalStart; i < signalEnd; i++)
  {
    signalMax = std::max(std::abs(data[i]), signalMax);
  }

  return signalMax / noiseMax;
}

bool projectionRequired(const Core::TimeWindow &tw,
                        const Catalog::Phase &ph,
                        const Catalog::Event &ev,
                        DataModel::ThreeComponents &tc,
                        DataModel::SensorLocation *&loc)
{
  string channelCodeRoot = getBandAndInstrumentCodes(ph.channelCode);
  string component       = getOrientationCode(ph.channelCode);

  if (std::set<string>{"Z", "N", "E", "R", "T"}.count(component) == 0)
  {
    // unknown projection
    return false;
  }

  loc = findSensorLocation(ph.networkCode, ph.stationCode, ph.locationCode,
                           tw.startTime());

  if (loc)
  {
    bool hasThreeComponents =
        getThreeComponents(tc, loc, channelCodeRoot.c_str(), tw.startTime());

    if ((tc.comps[ThreeComponents::Vertical] &&
         (tc.comps[ThreeComponents::Vertical]->code() == ph.channelCode)) ||
        (tc.comps[ThreeComponents::FirstHorizontal] &&
         (tc.comps[ThreeComponents::FirstHorizontal]->code() ==
          ph.channelCode)) ||
        (tc.comps[ThreeComponents::SecondHorizontal] &&
         (tc.comps[ThreeComponents::SecondHorizontal]->code() ==
          ph.channelCode)))
    {
      return false;
    }

    return hasThreeComponents;
  }

  SEISCOMP_DEBUG("Unable to fetch orientation information (%s)",
                 string(ph).c_str());
  return false;
}

GenericRecordPtr projectWaveform(const Core::TimeWindow &tw,
                                 const Catalog::Phase &ph,
                                 const Catalog::Event &ev,
                                 double tolerance,
                                 double minAvailability,
                                 const GenericRecordCPtr &tr1,
                                 const GenericRecordCPtr &tr2,
                                 const GenericRecordCPtr &tr3,
                                 const DataModel::ThreeComponents &tc,
                                 const DataModel::SensorLocation *loc)
{
  if (!tr1 || !tr2 || !tr3)
  {
    throw Exception("Cannot perform projection without all 3 components data");
  }

  string channelCodeRoot = getBandAndInstrumentCodes(ph.channelCode);
  string component       = getOrientationCode(ph.channelCode);

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
                   string(ph).c_str());
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
                   string(ph).c_str());
  }
  else
  {
    throw Exception("Internal logic error: This shouldn't happend.");
  }

  // The wrapper will direct 3 codes into the right slots using the
  // Stream configuration class and will finally use the transformation
  // operator. The advantage is that it will apply the configured gain.
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

  // Set gain to 0 since we do not need the gain for cross-correlation
  // and for consistency with rest of the API (e.g.
  // readWaveformFromRecordStream).
  streams[2].gain = 0;
  streams[1].gain = 0;
  streams[0].gain = 0;

  Rotator op(
      OpWrapper(streams, Operator::Transformation<double, 3>(transformation)));

  class DataStorer
  {
  public:
    DataStorer(const string &channelCode,
               const Core::TimeWindow &tw_,
               map<string, string> chCodeMap,
               double tolerance_,
               double minAvailability_)
        : chMap(chCodeMap[channelCode], channelCode), tw(tw_),
          tolerance(tolerance_), minAvailability(minAvailability_),
          seq(tw, tolerance)
    {}

    bool store(const Record *rec)
    {
      lock_guard<mutex> lock(_mtx);
      if (rec->channelCode() == chMap.first) seq.feed(rec);
      return true;
    }

    GenericRecordPtr get()
    {
      lock_guard<mutex> lock(_mtx);
      GenericRecordPtr trace =
          contiguousRecord(seq, tw, tolerance, minAvailability);
      if (trace) trace->setChannelCode(chMap.second);
      return trace;
    }

  private:
    const std::pair<string, string> chMap;
    const Core::TimeWindow tw;
    const double tolerance;
    const double minAvailability;
    TimeWindowBuffer seq;
    std::mutex _mtx;
  };

  DataStorer projectedData(ph.channelCode, tw, chCodeMap, tolerance,
                           minAvailability);

  // configure callback (called after a transformed record was created)
  op.setStoreFunc([&projectedData](const Record *rec) -> bool {
    return projectedData.store(rec);
  });
  op.feed(tr1.get());
  op.feed(tr2.get());
  op.feed(tr3.get());

  GenericRecordPtr trace = projectedData.get();
  if (!trace)
  {
    string msg = strf("No enough data for projection (%s)", string(ph).c_str());
    throw Exception(msg);
  }
  return trace;
}

GenericRecordCPtr Loader::get(const Core::TimeWindow &tw,
                              const Catalog::Phase &ph,
                              const Catalog::Event &ev)
{
  GenericRecordCPtr trace;

  if (!_recordStreamURL.empty())
  {
    DataModel::ThreeComponents tc;
    DataModel::SensorLocation *loc;
    bool projection = projectionRequired(tw, ph, ev, tc, loc);
    if (!projection)
    {
      try
      {
        trace = readWaveformFromRecordStream(
            _recordStreamURL, tw, ph.networkCode, ph.stationCode,
            ph.locationCode, ph.channelCode, _tolerance, _minAvailability);
        _counters_wf_downloaded++;
      }
      catch (exception &e)
      {
        SEISCOMP_DEBUG("%s", e.what());
      }
    }
    else
    {
      // If the waveform is not available, possibly a projection 123->ZNE
      // or ZNE->ZRT is required.
      auto loadWaveform = [this, &tw, &ph](const string &channelCode) {
        GenericRecordCPtr trace = readWaveformFromRecordStream(
            _recordStreamURL, tw, ph.networkCode, ph.stationCode,
            ph.locationCode, channelCode, _tolerance, _minAvailability);
        _counters_wf_downloaded++;
        return trace;
      };
      try
      {
        GenericRecordCPtr tr1 =
            loadWaveform(tc.comps[ThreeComponents::Vertical]->code());
        GenericRecordCPtr tr2 =
            loadWaveform(tc.comps[ThreeComponents::FirstHorizontal]->code());
        GenericRecordCPtr tr3 =
            loadWaveform(tc.comps[ThreeComponents::SecondHorizontal]->code());
        trace = projectWaveform(tw, ph, ev, _tolerance, _minAvailability, tr1,
                                tr2, tr3, tc, loc);
      }
      catch (exception &e)
      {
        SEISCOMP_DEBUG("%s", e.what());
      }
    }
  }

  if (!trace) _counters_wf_no_avail++;

  return trace;
}

GenericRecordCPtr Loader::get(const Core::TimeWindow &tw,
                              const Catalog::Phase &ph,
                              const Catalog::Event &ev,
                              bool demeaning,
                              const std::string &filterStr,
                              double resampleFreq)
{
  GenericRecordCPtr trace = get(tw, ph, ev);
  if (trace) trace = process(trace, demeaning, filterStr, resampleFreq);
  return trace;
}

GenericRecordPtr Loader::process(const GenericRecordCPtr &trace,
                                 bool demeaning,
                                 const std::string &filterStr,
                                 double resampleFreq) const
{
  GenericRecordPtr nonConstTrace(new GenericRecord(*trace));
  try
  {
    filter(*nonConstTrace, demeaning, filterStr, resampleFreq);
    return nonConstTrace;
  }
  catch (exception &e)
  {
    SEISCOMP_WARNING("Errow while filtering waveform: %s", e.what());
    return nullptr;
  }
}

GenericRecordCPtr DiskCachedLoader::get(const Core::TimeWindow &tw,
                                        const Catalog::Phase &ph,
                                        const Catalog::Event &ev)
{
  GenericRecordCPtr trace = getFromCache(tw, ph.networkCode, ph.stationCode,
                                         ph.locationCode, ph.channelCode);
  if (trace)
  {
    _counters_wf_cached++;
    return trace;
  }

  DataModel::ThreeComponents tc;
  DataModel::SensorLocation *loc;
  bool projection = projectionRequired(tw, ph, ev, tc, loc);
  if (!projection)
  {
    trace = _auxLdr->get(tw, ph, ev);
    if (trace)
    {
      storeInCache(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                   ph.channelCode, trace);
    }
    return trace;
  }

  auto loadWaveform = [this, &tw, &ph, &ev](const string &channelCode) {
    GenericRecordCPtr trace = getFromCache(tw, ph.networkCode, ph.stationCode,
                                           ph.locationCode, channelCode);
    if (trace)
    {
      _counters_wf_cached++;
      return trace;
    }
    Catalog::Phase tmpPh = ph;
    tmpPh.channelCode    = channelCode;
    trace                = _auxLdr->get(tw, tmpPh, ev);
    if (trace)
    {
      storeInCache(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                   channelCode, trace);
    }
    return trace;
  };
  GenericRecordCPtr tr1 =
      loadWaveform(tc.comps[ThreeComponents::Vertical]->code());
  GenericRecordCPtr tr2 =
      loadWaveform(tc.comps[ThreeComponents::FirstHorizontal]->code());
  GenericRecordCPtr tr3 =
      loadWaveform(tc.comps[ThreeComponents::SecondHorizontal]->code());
  try
  {
    trace = projectWaveform(tw, ph, ev, _tolerance, _minAvailability, tr1, tr2,
                            tr3, tc, loc);
  }
  catch (exception &e)
  {
    SEISCOMP_DEBUG("%s", e.what());
  }
  return trace;
}

bool DiskCachedLoader::isCached(const Core::TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev) const
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

GenericRecordCPtr MemCachedLoader::get(const Core::TimeWindow &tw,
                                       const Catalog::Phase &ph,
                                       const Catalog::Event &ev,
                                       bool demeaning,
                                       const std::string &filterStr,
                                       double resampleFreq)
{
  bool isCached;
  GenericRecordCPtr trace = getFromCache(tw, ph.networkCode, ph.stationCode,
                                         ph.locationCode, ph.channelCode);
  if (trace)
  {
    isCached = true;
    _counters_wf_cached++;
  }
  else
  {
    isCached = false;
    trace    = _auxLdr->get(tw, ph, ev, demeaning, filterStr, resampleFreq);
  }

  if (trace && !isCached)
  {
    storeInCache(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                 ph.channelCode, trace);

    // Dump waveforms when loaded the first time for debugging
    if (!_wfDebugDir.empty())
    {
      string ext = (ph.procInfo.source == Catalog::Phase::Source::THEORETICAL)
                       ? "theoretical"
                       : (ph.isManual ? "manual" : "automatic");
      writeTrace(trace, waveformDebugPath(_wfDebugDir, ev, ph, ext));
    }
  }

  return trace;
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

  // Make sure to load at least `_traceMinLenh` seconds of the waveform.
  if (_beforePickLen > 0 || _afterPickLen > 0)
  {
    const Core::TimeSpan additionalTimeBefore(_beforePickLen);
    const Core::TimeSpan additionalTimeAfter(_afterPickLen);

    if (twToLoad.startTime() > pickTime - additionalTimeBefore)
      twToLoad.setStartTime(pickTime - additionalTimeBefore);

    if (twToLoad.endTime() < pickTime + additionalTimeAfter)
      twToLoad.setEndTime(pickTime + additionalTimeAfter);
  }

  // round the start/end time to the nearest second
  twToLoad.setStartTime(Core::Time(twToLoad.startTime().seconds(), 0));
  twToLoad.setEndTime(twToLoad.endTime().microseconds() > 0
                          ? Core::Time(twToLoad.endTime().seconds() + 1, 0)
                          : Core::Time(twToLoad.endTime().seconds(), 0));
  return twToLoad;
}

GenericRecordCPtr ExtraLenLoader::get(const Core::TimeWindow &tw,
                                      const Catalog::Phase &ph,
                                      const Catalog::Event &ev)
{
  const Core::TimeWindow twToLoad = traceTimeWindowToLoad(tw, ph.time);
  GenericRecordCPtr trace         = _auxLdr->get(twToLoad, ph, ev);
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

  // Check if we have already excluded the trace's SNR.
  if (_snrExcludedWfs.count(wfId) != 0)
  {
    return nullptr;
  }

  const Core::TimeWindow twToLoad = tw | snrTimeWindow(ph.time);
  GenericRecordCPtr trace =
      _auxLdr->get(twToLoad, ph, ev, demeaning, filterStr, resampleFreq);

  if (!trace) return nullptr;

  // Check if we have already validated the trace's SNR, if not do it now.
  if (_snrGoodWfs.count(wfId) == 0)
  {
    if (!goodSnr(trace, ph.time))
    {
      _snrExcludedWfs.insert(wfId);
      SEISCOMP_DEBUG("Trace has too low SNR(%s)", string(ph).c_str());
      // Dump SNR low traces (debugging).
      if (!_wfDebugDir.empty())
      {
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

BatchLoader::BatchLoader(const std::string &recordStream) : Loader(recordStream)
{
  _dataLoaded = false;
}

GenericRecordCPtr BatchLoader::get(const Core::TimeWindow &tw,
                                   const Catalog::Phase &ph,
                                   const Catalog::Event &ev)
{
  if (!_dataLoaded)
  {
    request(tw, ph, ev);
    return nullptr;
  }

  DataModel::ThreeComponents tc;
  DataModel::SensorLocation *loc;
  bool projection = projectionRequired(tw, ph, ev, tc, loc);

  GenericRecordCPtr trace;

  auto getTrace = [this, &tw, &ph](const string &channelCode) {
    const string wfId = waveformId(tw, ph.networkCode, ph.stationCode,
                                   ph.locationCode, channelCode);
    const auto it     = _waveforms.find(wfId);
    return it != _waveforms.end() ? it->second : nullptr;
  };

  if (!projection)
  {
    trace = getTrace(ph.channelCode);
  }
  else
  {
    GenericRecordCPtr tr1 =
        getTrace(tc.comps[ThreeComponents::Vertical]->code());
    GenericRecordCPtr tr2 =
        getTrace(tc.comps[ThreeComponents::FirstHorizontal]->code());
    GenericRecordCPtr tr3 =
        getTrace(tc.comps[ThreeComponents::SecondHorizontal]->code());
    try
    {
      trace = projectWaveform(tw, ph, ev, _tolerance, _minAvailability, tr1,
                              tr2, tr3, tc, loc);
    }
    catch (exception &e)
    {
      SEISCOMP_DEBUG("%s", e.what());
    }
  }
  return trace;
}

void BatchLoader::request(const Core::TimeWindow &tw,
                          const Catalog::Phase &ph,
                          const Catalog::Event &ev)
{
  if (_dataLoaded)
  {
    throw Exception("Cannot request more traces after they have been loaded");
  }

  DataModel::ThreeComponents tc;
  DataModel::SensorLocation *loc;
  bool projection = projectionRequired(tw, ph, ev, tc, loc);

  auto requestTrace = [this, &tw, &ph](const string &channelCode) {
    string streamID = ph.networkCode + "." + ph.stationCode + "." +
                      ph.locationCode + "." + channelCode;
    auto eqlrng = _streamMap.equal_range(streamID);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      pair<const Core::TimeWindow, TimeWindowBuffer> &pair = it->second;
      if (pair.first == tw) return; // skip duplicated requests
    }
    pair<const Core::TimeWindow, TimeWindowBuffer> pair(
        tw, TimeWindowBuffer(tw, _tolerance));
    _streamMap.emplace(streamID, pair);
  };

  if (!projection)
  {
    requestTrace(ph.channelCode);
  }
  else
  {
    requestTrace(tc.comps[ThreeComponents::Vertical]->code());
    requestTrace(tc.comps[ThreeComponents::FirstHorizontal]->code());
    requestTrace(tc.comps[ThreeComponents::SecondHorizontal]->code());
  }
}

void BatchLoader::load()
{
  if (_dataLoaded) return;

  if (!_streamMap.empty())
  {
    std::unordered_multimap<std::string, const Core::TimeWindow> tmpMap;
    for (const auto &kv : _streamMap)
    {
      const string &streamID                                     = kv.first;
      const pair<const Core::TimeWindow, TimeWindowBuffer> &pair = kv.second;
      tmpMap.emplace(streamID, pair.first);
    }

    //
    // The reason of this loop is that RecordStream can only handle one single
    // time window per stream
    //
    while (!tmpMap.empty())
    {
      _rs = IO::RecordStream::Open(_recordStreamURL.c_str());
      if (_rs == nullptr)
      {
        SEISCOMP_ERROR("Cannot open RecordStream: %s",
                       _recordStreamURL.c_str());
        break;
      }

      //
      // Convert multiple time windows requests close to each others to a
      // single RecordStream request on the stream
      //
      for (auto it = tmpMap.begin(), end = tmpMap.end();
           it != end;) // loop by stream
      {
        const string streamID = it->first;
        Core::TimeWindow contiguousRequest;
        auto eqlrng = tmpMap.equal_range(streamID);
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
            it2 = tmpMap.erase(it2);
          else
            it2++;
        }
        static const std::regex dot("\\.", std::regex::optimize);
        const std::vector<std::string> tokens(splitString(streamID, dot));
        _rs->addStream(tokens.at(0), tokens.at(1), tokens.at(2), tokens.at(3),
                       contiguousRequest.startTime(),
                       contiguousRequest.endTime());
        it = eqlrng.second;
      }

      //
      // Collect data records and feed them to the corresponding
      // TimeWindowBuffer
      //
      IO::RecordInput inp(_rs.get(), Array::DOUBLE, Record::DATA_ONLY);
      RecordPtr rec;
      while (rec = inp.next())
      {
        auto eqlrng = _streamMap.equal_range(rec->streamID());
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
          pair<const Core::TimeWindow, TimeWindowBuffer> &pair = it->second;
          TimeWindowBuffer &seq                                = pair.second;
          seq.feed(rec.get());
        }
      }
      _rs->close();
    }

    //
    // Convert data records to contiguous waveforms
    //
    for (auto &kv : _streamMap)
    {
      const string &streamID                               = kv.first;
      pair<const Core::TimeWindow, TimeWindowBuffer> &pair = kv.second;
      const Core::TimeWindow &tw                           = pair.first;
      TimeWindowBuffer &seq                                = pair.second;
      GenericRecordPtr trace =
          contiguousRecord(seq, tw, _tolerance, _minAvailability);
      if (!trace)
      {
        SEISCOMP_DEBUG("Cannnot load trace, data availability %.2f%%"
                       "(stream %s from %s length %.2f sec)",
                       seq.availability(), streamID.c_str(),
                       tw.startTime().iso().c_str(), tw.length());
        _counters_wf_no_avail++;
        continue;
      }
      const string wfId =
          waveformId(tw, trace->networkCode(), trace->stationCode(),
                     trace->locationCode(), trace->channelCode());
      _waveforms[wfId] = trace;
      _counters_wf_downloaded++;
    }
    SEISCOMP_INFO("Fetched %u/%lu waveforms, not available %u",
                  _counters_wf_downloaded, _streamMap.size(),
                  _counters_wf_no_avail);
    _streamMap.clear();
  }

  _dataLoaded = true;
  _rs         = nullptr;
}

} // namespace Waveform
} // namespace HDD

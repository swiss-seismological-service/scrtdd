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

#define BOOST_TEST_MODULE libhdd
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>

#include "hdd/dd.h"
#include <cstring>
#include <ostream>

using namespace std;
using namespace HDD;
namespace bdata = boost::unit_test::data;

namespace std {
// Required by BOOST TEST in the std namespace
std::ostream &operator<<(std::ostream &os, const HDD::Trace &tr)
{
  os << tr.streamID() << "@" << UTCClock::toString(tr.startTime()) << "~"
     << UTCClock::toString(tr.endTime());
  return os;
}
} // namespace std

namespace {

//
// Synthetic traces: a signal (sum of 1 Hz and 11 Hz sines, a short 6 Hz sine,
// a triangle) surrounded by zeros. Each trace is built with several sampling
// frequencies.
//
Trace buildSyntheticTrace1(double samplingFrequency)
{
  Trace tr("N1", "ST1", "", "EHE", UTCClock::fromDate(1981, 1, 9, 21, 56, 4, 1),
           samplingFrequency);
  vector<double> samples(tr.samplingFrequency() * 3, 0);
  for (int i = 0; i < tr.samplingFrequency(); i++)
  {
    double value =
        std::sin(2 * M_PI * i * 11 / tr.samplingFrequency()); // 11 Hz sin
    value += std::sin(2 * M_PI * i / tr.samplingFrequency()); // 1 Hz sin
    samples[samples.size() / 2 + i] = value;
  }
  tr.setData(std::move(samples));
  return tr;
}

Trace buildSyntheticTrace2(double samplingFrequency)
{
  Trace tr("N2", "ST2", "", "BHZ",
           UTCClock::fromDate(2021, 3, 19, 13, 56, 4, 2), samplingFrequency);
  vector<double> samples(tr.samplingFrequency() * 2, 0);
  for (int i = 0; i < tr.samplingFrequency() / 6; i++)
  {
    double value =
        std::sin(2 * M_PI * i * 6 / tr.samplingFrequency()); // 6 Hz sin
    samples[samples.size() / 2 + i] = value;
  }
  tr.setData(std::move(samples));
  return tr;
}

Trace buildSyntheticTrace3(double samplingFrequency)
{
  Trace tr("N3", "ST3", "", "GHZ", UTCClock::fromDate(1998, 7, 5, 6, 12, 9, 3),
           samplingFrequency);
  vector<double> samples(tr.samplingFrequency() * 3.5, 0);
  for (int i = 0; i < tr.samplingFrequency(); ++i)
  {
    double value =
        (i <= tr.samplingFrequency() / 2) ? i : tr.samplingFrequency() - i;
    samples[samples.size() / 2 + i] = value;
  }
  tr.setData(std::move(samples));
  return tr;
}

void scaleTrace(Trace &tr, double constant, double scaler)
{
  double *samples = tr.data();
  for (size_t i = 0; i < tr.sampleCount(); ++i)
  {
    samples[i] *= scaler;
    samples[i] += constant;
  }
}

void trimTrace(Trace &tr, double trimStart, double trimEnd)
{
  TimeWindow tw = tr.timeWindow().trim(secToDur(trimStart), secToDur(trimEnd));
  BOOST_REQUIRE(tr.slice(tw));
}

Trace alterTrace(const Trace &tr,
                 double trimStart,
                 double trimEnd,
                 double constant,
                 double scaler)
{
  Trace newTr(tr);
  if (trimStart > 0 || trimEnd > 0) trimTrace(newTr, trimStart, trimEnd);
  if (constant != 0 || scaler != 1) scaleTrace(newTr, constant, scaler);
  return newTr;
}

// Cross-correlate the traces and check the lag and that the correlation
// coefficient is close to 1 (or -1 for opposite polarity)
void testXCorrTraces(const Trace &tr1, const Trace &tr2, double expectedLag)
{
  double delayOut, coeffOut;
  DD::xcorr(tr1, tr2,
            (std::abs(durToSec(tr1.timeWindow().length() -
                              tr2.timeWindow().length())) + 
            2. / tr1.samplingFrequency()) / 2.0,
            delayOut, coeffOut);
  BOOST_CHECK_SMALL(expectedLag - delayOut,
                    2.0 / tr1.samplingFrequency()); // 2 samples tolerance
  BOOST_CHECK_SMALL(1.0 - std::abs(coeffOut), 0.1); // positive or negative CC
}

// Cross-correlate `tr1` with a copy of it trimmed at one end, so that the
// trace center (where the xcorr assumes the pick) moves by `timeShift`, and
// with an offset and scaling applied: the lag must be `timeShift`
void testXCorrTrace(const Trace &tr1,
                    double timeShift,
                    double constant,
                    double scaler)
{
  Trace tr2 =
      timeShift >= 0
          ? alterTrace(tr1, 0, timeShift * 2, constant, scaler)
          : alterTrace(tr1, std::abs(timeShift * 2), 0, constant, scaler);
  testXCorrTraces(tr1, tr2, timeShift);
}

// Several time shifts, offsets, scalings and polarity changes. `s` varies the
// offset and scaling between the traces
void testXCorr(const Trace &trace, int s)
{
  testXCorrTrace(trace, 0, 0, 1);
  testXCorrTrace(trace, 0, 0, -1);
  testXCorrTrace(trace, 25. / trace.samplingFrequency(), -1 * s, 100 * s);
  testXCorrTrace(trace, 25. / trace.samplingFrequency(), 1 * s, 0.01 * s);
  testXCorrTrace(trace, 10. / trace.samplingFrequency(), -10 * s, 0.1 * s);
  testXCorrTrace(trace, 10. / trace.samplingFrequency(), 10 * s, 10 * s);
  testXCorrTrace(trace, -25. / trace.samplingFrequency(), -1 * s, 100 * s);
  testXCorrTrace(trace, -25. / trace.samplingFrequency(), 1 * s, 0.01 * s);
  testXCorrTrace(trace, -10. / trace.samplingFrequency(), -10 * s, 0.1 * s);
  testXCorrTrace(trace, -10. / trace.samplingFrequency(), 10 * s, 10 * s);
}

// The same trace sampled at different frequencies and then resampled to a
// common frequency must cross-correlate with zero lag
void testResampling(const vector<Trace> &traces)
{
  // traces vector contains identical traces, but different sampling rates
  for (size_t i = 0; i < traces.size() - 1; i++)
  {
    for (size_t j = i + 1; j < traces.size(); j++)
    {
      Trace tr1(traces[i]);
      Trace tr2(traces[j]);
      Waveform::resample(tr1, 1000);
      Waveform::resample(tr2, 1000);
      testXCorrTraces(tr1, tr2, 0);

      tr1 = traces[i];
      tr2 = traces[j];
      Waveform::resample(tr1, 60);
      Waveform::resample(tr2, 60);
      testXCorrTraces(tr1, tr2, 0);

      tr1 = traces[i];
      tr2 = traces[j];
      Waveform::resample(tr1, tr2.samplingFrequency());
      testXCorrTraces(tr1, tr2, 0);

      tr1 = traces[i];
      tr2 = traces[j];
      Waveform::resample(tr2, tr1.samplingFrequency());
      testXCorrTraces(tr1, tr2, 0);
    }
  }
}

const vector<Trace> synthetic1Traces = {
    buildSyntheticTrace1(160), buildSyntheticTrace1(158),
    buildSyntheticTrace1(80), buildSyntheticTrace1(83)};

const vector<Trace> synthetic2Traces = {
    buildSyntheticTrace2(300), buildSyntheticTrace2(280),
    buildSyntheticTrace2(150), buildSyntheticTrace2(143)};

const vector<Trace> synthetic3Traces = {
    buildSyntheticTrace3(120), buildSyntheticTrace3(119),
    buildSyntheticTrace3(60), buildSyntheticTrace3(67)};

} // namespace

// Verify the cross-correlation lag and coefficient for each synthetic trace
// and sampling frequency
BOOST_DATA_TEST_CASE(test_xcorr1,
                     bdata::xrange(synthetic1Traces.size()) ^
                         bdata::make(synthetic1Traces),
                     i, trace)
{
  testXCorr(trace, i + 1);
}

BOOST_DATA_TEST_CASE(test_xcorr2,
                     bdata::xrange(synthetic2Traces.size()) ^
                         bdata::make(synthetic2Traces),
                     i, trace)
{
  testXCorr(trace, i + 1);
}

BOOST_DATA_TEST_CASE(test_xcorr3,
                     bdata::xrange(synthetic3Traces.size()) ^
                         bdata::make(synthetic3Traces),
                     i, trace)
{
  testXCorr(trace, i + 1);
}

// Verify that a trace can be sliced with any window within the trace, also
// when the window starts or ends with the trace (the trace start and end times
// are rounded to the clock resolution)
BOOST_AUTO_TEST_CASE(test_slice)
{
  for (const vector<Trace> *traces :
       {&synthetic1Traces, &synthetic2Traces, &synthetic3Traces})
  {
    for (const Trace &trace : *traces)
    {
      const double trimLen = 50. / trace.samplingFrequency(); // 50 samples

      // trim the start only, the end only, both
      for (const auto &trim : vector<pair<double, double>>{
               {trimLen, 0}, {0, trimLen}, {trimLen, trimLen}})
      {
        BOOST_TEST_MESSAGE(strf("Slice %s %.0f Hz trimming %.3f %.3f sec",
                                trace.streamID().c_str(),
                                trace.samplingFrequency(), trim.first,
                                trim.second));
        Trace tr(trace);
        const TimeWindow tw = tr.timeWindow().trim(secToDur(trim.first),
                                                   secToDur(trim.second));
        BOOST_REQUIRE(tr.slice(tw));

        // the window times are rounded too: allow 1 more sample per side
        const unsigned trimmedSides = (trim.first > 0) + (trim.second > 0);
        const size_t expected = trace.sampleCount() - 50 * trimmedSides;
        BOOST_CHECK_GE(tr.sampleCount(), expected);
        BOOST_CHECK_LE(tr.sampleCount(), expected + trimmedSides);
      }
    }
  }
}

// Verify that resampling preserves the signal timing
BOOST_AUTO_TEST_CASE(test_resampling1)
{
  testResampling(synthetic1Traces);
}

BOOST_AUTO_TEST_CASE(test_resampling2)
{
  testResampling(synthetic2Traces);
}

BOOST_AUTO_TEST_CASE(test_resampling3)
{
  testResampling(synthetic3Traces);
}


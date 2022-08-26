#define SEISCOMP_TEST_MODULE hdd
#include <seiscomp/unittest/unittests.h>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>

#include "hdd/dd.h"
#include "hdd/adapters/scwaveform.h"
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

void writeTrace(const Trace &trace, const std::string &filename)
{
  HDD::SCAdapter::WaveformProxy().writeTrace(trace, filename);
}

Trace readTrace(const std::string &filename)
{
  return *HDD::SCAdapter::WaveformProxy().readTrace(filename);
}

void filter(Trace &trace,
              bool demeaning,
              const std::string &filterStr,
              double resampleFreq)
{
  Waveform::BasicProcessor(shared_ptr<Waveform::Proxy>(
        new HDD::SCAdapter::WaveformProxy()), nullptr,  0)
    .filter(trace, demeaning, filterStr, resampleFreq);
}

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
  // writeTrace(tr, strf("syntheticTrace1_%.fHz.mseed", tr.samplingFrequency()));
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
  // writeTrace(tr, strf("syntheticTrace2_%.fHz.mseed", tr.samplingFrequency()));
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
  // writeTrace(tr, strf("syntheticTrace3_%.fHz.mseed", tr.samplingFrequency()));
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
  tr.slice(tw);
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

void testTracesEqual(const Trace &tr1, const Trace &tr2)
{
  BOOST_CHECK_EQUAL(tr1.streamID(), tr2.streamID());
  BOOST_CHECK(tr1.startTime() == tr2.startTime());
  BOOST_CHECK(tr1.endTime() == tr2.endTime());
  BOOST_CHECK_EQUAL(tr1.samplingFrequency(), tr2.samplingFrequency());
  BOOST_CHECK_EQUAL(tr1.sampleCount(), tr2.sampleCount());
  BOOST_CHECK(tr1.elementSize() == tr2.elementSize());
  BOOST_CHECK(
      std::memcmp(tr1.data(), tr2.data(),
                  std::min(tr1.sampleCount() * tr1.elementSize(),
                           tr2.sampleCount() * tr2.elementSize())) == 0);
}

void testReadWriteTrace(const Trace &tr)
{
  const string filename = "test_read_write.mseed";
  if (pathExists(filename)) removePath(filename);
  BOOST_REQUIRE(!pathExists(filename));
  writeTrace(tr, filename);
  Trace cmpTr = readTrace(filename);
  testTracesEqual(tr, cmpTr);
  removePath(filename);
}

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

void testXCorr(const Trace &trace, int s)
{
  testXCorrTrace(trace, 0, 0, 1);
  testXCorrTrace(trace, 0, 0, -1);
  testXCorrTrace(trace, 25. / trace.samplingFrequency(), -1 * s, 100 * s);
  testXCorrTrace(trace, 25. / trace.samplingFrequency(), 1 * s, 0.01 * s);
  testXCorrTrace(trace, 10. / trace.samplingFrequency(), -10 * s, 0.1 * s);
  testXCorrTrace(trace, 10. / trace.samplingFrequency(), 10 * s, 10 * s);
}

void testReampling(const vector<Trace> &traces)
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

      // writeTrace(tr1, strf("testReamplingTr1_%s_%.f-%.fHz.mseed",
      //            tr1.networkCode().c_str(),
      //            traces[i].samplingFrequency(),
      //            tr1.samplingFrequency()));
      // writeTrace(tr2,
      //            strf("testReamplingTr2_%s_%.f-%.fHz.mseed",
      //            tr2.networkCode().c_str(),
      //            traces[j].samplingFrequency(),
      //            tr2.samplingFrequency()));

      tr1 = traces[i];
      tr2 = traces[j];
      Waveform::resample(tr1, 60);
      Waveform::resample(tr2, 60);
      testXCorrTraces(tr1, tr2, 0);

      // writeTrace(tr1, strf("testReamplingTr1_%s_%.f-%.fHz.mseed",
      //            tr1.networkCode().c_str(),
      //            traces[i].samplingFrequency(),
      //            tr1.samplingFrequency()));
      // writeTrace(tr2,
      //            strf("testReamplingTr2_%s_%.f-%.fHz.mseed",
      //            tr2.networkCode().c_str(),
      //            traces[j].samplingFrequency(),
      //            tr2.samplingFrequency()));

      tr1 = traces[i];
      tr2 = traces[j];
      Waveform::resample(tr1, tr2.samplingFrequency());
      testXCorrTraces(tr1, tr2, 0);

      // writeTrace(tr1, strf("testReamplingTr1_%s_%.f-%.fHz.mseed",
      //            tr1.networkCode().c_str(),
      //            traces[i].samplingFrequency(),
      //            tr1.samplingFrequency()));

      tr1 = traces[i];
      tr2 = traces[j];
      Waveform::resample(tr2, tr1.samplingFrequency());
      testXCorrTraces(tr1, tr2, 0);

      // writeTrace(tr2, strf("testReamplingTr2_%s_%.f-%.fHz.mseed",
      //            tr2.networkCode().c_str(),
      //            traces[j].samplingFrequency(),
      //            tr2.samplingFrequency()));
    }
  }
}

void testFiltering(const vector<Trace> &traces)
{
  const string filterStr = "ITAPER(1)>>BW_HLP(1,1,20)";

  // traces vector contains identical traces, but different sampling rates
  for (size_t i = 0; i < traces.size() - 1; i++)
  {
    for (size_t j = i + 1; j < traces.size(); j++)
    {
      Trace tr1 = alterTrace(traces[i], 0, 0, 10, 1);
      Trace tr2 = alterTrace(traces[j], 0, 0, -10, 1);
      filter(tr1, true, filterStr, 500);
      filter(tr2, true, filterStr, 500);
      testXCorrTraces(tr1, tr2, 0);

      // writeTrace(tr1,strf("testFilteringTr1_%s_%.f-%.fHz.mseed",
      //            tr1.networkCode().c_str(),
      //            traces[i].samplingFrequency(),
      //            tr1.samplingFrequency()));
      // writeTrace(tr2, strf("testFilteringTr2_%s_%.f-%.fHz.mseed",
      //            tr2.networkCode().c_str(),
      //            traces[j].samplingFrequency(),
      //            tr2.samplingFrequency()));

      tr1 = alterTrace(traces[i], 0, 0, 10, 1);
      tr2 = alterTrace(traces[j], 0, 0, -10, 1);
      filter(tr1, true, filterStr, 50);
      filter(tr2, true, filterStr, 50);
      testXCorrTraces(tr1, tr2, 0);

      // writeTrace(tr1, strf("testFilteringTr1_%s_%.f-%.fHz.mseed",
      //            tr1.networkCode().c_str(),
      //            traces[i].samplingFrequency(),
      //            tr1.samplingFrequency()));
      // writeTrace(tr2, strf("testFilteringTr2_%s_%.f-%.fHz.mseed",
      //            tr2.networkCode().c_str(),
      //            traces[j].samplingFrequency(),
      //            tr2.samplingFrequency()));
    }
  }
}

void testSnrSynthetic(const Trace &trace)
{
  double snr = Waveform::computeSnr(
      trace, trace.startTime() + (trace.timeWindow().length() / 2), -0.5, 0, 0,
      0.5);
  BOOST_CHECK(snr > 2);

  snr = Waveform::computeSnr(
      trace, trace.startTime() + (trace.timeWindow().length() / 2), 0, 0.5, 0.5,
      1);
  BOOST_CHECK(snr < 2);
}

const vector<Trace> realTraces = {
    readTrace("./data/waveform/xcorr1.mseed"),
    readTrace("./data/waveform/xcorr2.mseed"),
    readTrace("./data/waveform/xcorr3.mseed"),
    readTrace("./data/waveform/xcorr4.mseed"),
    readTrace("./data/waveform/xcorr5.mseed"),
    readTrace("./data/waveform/xcorr6.mseed")};

const vector<Trace> badSnrTraces = {
    readTrace("./data/waveform/snr1.mseed"),
    readTrace("./data/waveform/snr2.mseed"),
    readTrace("./data/waveform/snr3.mseed"),
    readTrace("./data/waveform/snr4.mseed"),
    readTrace("./data/waveform/snr5.mseed"),
    readTrace("./data/waveform/snr6.mseed"),
    readTrace("./data/waveform/snr7.mseed")};

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

BOOST_DATA_TEST_CASE(test_read_write1, bdata::make(realTraces), trace)
{
  testReadWriteTrace(trace);
}

BOOST_DATA_TEST_CASE(test_read_write2, bdata::make(badSnrTraces), trace)
{
  testReadWriteTrace(trace);
}

BOOST_DATA_TEST_CASE(test_read_write3, bdata::make(synthetic1Traces), trace)
{
  testReadWriteTrace(trace);
}

BOOST_DATA_TEST_CASE(test_read_write4, bdata::make(synthetic2Traces), trace)
{
  testReadWriteTrace(trace);
}

BOOST_DATA_TEST_CASE(test_read_write5, bdata::make(synthetic3Traces), trace)
{
  testReadWriteTrace(trace);
}

BOOST_DATA_TEST_CASE(test_xcorr1,
                     bdata::xrange(realTraces.size()) ^ bdata::make(realTraces),
                     i,
                     trace)
{
  testXCorr(trace, i + 1);
}

BOOST_DATA_TEST_CASE(test_xcorr2,
                     bdata::xrange(synthetic1Traces.size()) ^
                         bdata::make(synthetic1Traces),
                     i, trace)
{
  testXCorr(trace, i + 1);
}

BOOST_DATA_TEST_CASE(test_xcorr3,
                     bdata::xrange(synthetic2Traces.size()) ^
                         bdata::make(synthetic2Traces),
                     i, trace)
{
  testXCorr(trace, i + 1);
}

BOOST_DATA_TEST_CASE(test_xcorr4,
                     bdata::xrange(synthetic3Traces.size()) ^
                         bdata::make(synthetic3Traces),
                     i, trace)
{
  testXCorr(trace, i + 1);
}

BOOST_AUTO_TEST_CASE(test_resampling1)
{
  testReampling(synthetic1Traces);
}

BOOST_AUTO_TEST_CASE(test_resampling2)
{
  testReampling(synthetic2Traces);
}

BOOST_AUTO_TEST_CASE(test_resampling3)
{
  testReampling(synthetic3Traces);
}

BOOST_AUTO_TEST_CASE(test_filtering1)
{
  testFiltering(synthetic1Traces);
}

BOOST_AUTO_TEST_CASE(test_filtering2)
{
  testFiltering(synthetic2Traces);
}

BOOST_DATA_TEST_CASE(test_snr_bad, bdata::make(badSnrTraces), trace)
{
  double snr = Waveform::computeSnr(
      trace, trace.startTime() + (trace.timeWindow().length() / 2), -3, -0.350,
      -0.350, 0.350);
  BOOST_CHECK(snr < 2);
}

BOOST_DATA_TEST_CASE(test_snr_ok, bdata::make(realTraces), trace)
{
  double snr = Waveform::computeSnr(
      trace, trace.startTime() + (trace.timeWindow().length() / 2), -0.7, -0.1,
      -0.1, 0.5);
  BOOST_CHECK(snr >= 2);
}

BOOST_DATA_TEST_CASE(test_snr_synthetic1, bdata::make(synthetic1Traces), trace)
{
  testSnrSynthetic(trace);
}

BOOST_DATA_TEST_CASE(test_snr_synthetic2, bdata::make(synthetic2Traces), trace)
{
  testSnrSynthetic(trace);
}

BOOST_DATA_TEST_CASE(test_snr_synthetic3, bdata::make(synthetic3Traces), trace)
{
  testSnrSynthetic(trace);
}

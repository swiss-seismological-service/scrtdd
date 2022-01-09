#define SEISCOMP_TEST_MODULE hdd
#include <seiscomp/unittest/unittests.h>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>

#include "waveform.h"
#include <cstring>
#include <seiscomp3/core/genericrecord.h>
#include <seiscomp3/core/typedarray.h>

using namespace std;
using namespace HDD;
using namespace Seiscomp;
using namespace Seiscomp::Core;
namespace bdata = boost::unit_test::data;

namespace {

GenericRecordPtr buildSyntheticTrace1(double samplingFrequency)
{
  GenericRecordPtr tr  = new GenericRecord("N1", "ST1", "", "EHE",
                                          Core::Time(1981, 1, 9, 21, 56, 4, 1),
                                          samplingFrequency, 10, Array::DOUBLE);
  DoubleArray *samples = new DoubleArray(tr->samplingFrequency() * 3);
  samples->fill(0);
  for (int i = 0; i < tr->samplingFrequency(); i++)
  {
    double value =
        std::sin(2 * M_PI * i * 11 / tr->samplingFrequency()); // 11 Hz sin
    value += std::sin(2 * M_PI * i / tr->samplingFrequency()); // 1 Hz sin
    samples->set(samples->size() / 2 + i, value);
  }
  tr->setData(samples);
  // writeTrace(tr, strf("syntheticTrace1_%.fHz.mseed", tr->samplingFrequency()));
  return tr;
}

GenericRecordPtr buildSyntheticTrace2(double samplingFrequency)
{
  GenericRecordPtr tr  = new GenericRecord("N2", "ST2", "", "BHZ",
                                          Core::Time(2021, 3, 19, 13, 56, 4, 2),
                                          samplingFrequency, 30, Array::DOUBLE);
  DoubleArray *samples = new DoubleArray(tr->samplingFrequency() * 2);
  samples->fill(0);
  for (int i = 0; i < tr->samplingFrequency() / 6; i++)
  {
    double value =
        std::sin(2 * M_PI * i * 6 / tr->samplingFrequency()); // 6 Hz sin
    samples->set(samples->size() / 2 + i, value);
  }
  tr->setData(samples);
  // writeTrace(tr, strf("syntheticTrace2_%.fHz.mseed", tr->samplingFrequency()));
  return tr;
}

GenericRecordPtr buildSyntheticTrace3(double samplingFrequency)
{
  GenericRecordPtr tr  = new GenericRecord("N3", "ST3", "", "GHZ",
                                          Core::Time(1998, 7, 5, 6, 12, 9, 3),
                                          samplingFrequency, 20, Array::DOUBLE);
  DoubleArray *samples = new DoubleArray(tr->samplingFrequency() * 3.5);
  samples->fill(0);
  for (int i = 0; i < tr->samplingFrequency(); ++i)
  {
    double value =
        (i <= tr->samplingFrequency() / 2) ? i : tr->samplingFrequency() - i;
    samples->set(samples->size() / 2 + i, value);
  }
  tr->setData(samples);
  // writeTrace(tr, strf("syntheticTrace3_%.fHz.mseed", tr->samplingFrequency()));
  return tr;
}

void scaleTrace(GenericRecordPtr &tr, double constant, double scaler)
{
  DoubleArray *samples = DoubleArray::Cast(tr->data());
  for (int i = 0; i < samples->size(); ++i)
  {
    double value = (*samples)[i];
    value *= scaler;
    value += constant;
    samples->set(i, value);
  }
  tr->dataUpdated();
}

void trimTrace(GenericRecordPtr &tr, double trimStart, double trimEnd)
{
  int start       = trimStart * tr->samplingFrequency();
  int end         = tr->data()->size() - trimEnd * tr->samplingFrequency();
  ArrayPtr sliced = tr->data()->slice(start, end);
  tr->setData(sliced.get());
}

GenericRecordPtr alterTrace(const GenericRecordCPtr &tr,
                            double trimStart,
                            double trimEnd,
                            double constant,
                            double scaler)
{
  GenericRecordPtr newTr(new GenericRecord(*tr));
  if (trimStart > 0 || trimEnd > 0) trimTrace(newTr, trimStart, trimEnd);
  if (constant != 0 || scaler != 1) scaleTrace(newTr, constant, scaler);
  return newTr;
}

void testTracesEqual(const GenericRecordCPtr &tr1, const GenericRecordCPtr &tr2)
{
  BOOST_REQUIRE(tr1 && tr2);
  BOOST_CHECK_EQUAL(tr1->streamID(), tr2->streamID());
  BOOST_CHECK(tr1->startTime() == tr2->startTime());
  BOOST_CHECK(tr1->endTime() == tr2->endTime());
  BOOST_CHECK_EQUAL(tr1->samplingFrequency(), tr2->samplingFrequency());
  BOOST_CHECK_EQUAL(tr1->timingQuality(), tr2->timingQuality());
  BOOST_CHECK_EQUAL(tr1->dataType(), tr2->dataType());
  BOOST_CHECK_EQUAL(tr1->sampleCount(), tr2->sampleCount());
  BOOST_CHECK(tr1->data()->size() == tr2->data()->size());
  BOOST_CHECK(tr1->data()->elementSize() == tr2->data()->elementSize());
  BOOST_CHECK(
      std::memcmp(tr1->data()->data(), tr2->data()->data(),
                  std::min(tr1->data()->elementSize() * tr1->data()->size(),
                           tr2->data()->elementSize() * tr2->data()->size())) ==
      0);
}

void testReadWriteTrace(const GenericRecordCPtr &tr)
{
  const string filename = "test_read_write.mseed";
  if (pathExists(filename)) removePath(filename);
  BOOST_REQUIRE(!pathExists(filename));
  HDD::Waveform::writeTrace(tr, filename);
  GenericRecordPtr cmpTr = HDD::Waveform::readTrace(filename);
  testTracesEqual(tr, cmpTr);
  removePath(filename);
}

void testXCorrTraces(const GenericRecordCPtr &tr1,
                     const GenericRecordCPtr &tr2,
                     double expectedLag)
{
  double delayOut, coeffOut;
  BOOST_CHECK(HDD::Waveform::xcorr(
      tr1, tr2,
      std::abs(tr1->timeWindow().length() - tr2->timeWindow().length()) / 2.0,
      true, delayOut, coeffOut));
  BOOST_CHECK_SMALL(expectedLag - delayOut,
                    2.0 / tr1->samplingFrequency()); // 2 samples tolerance
  BOOST_CHECK_SMALL(1.0 - std::abs(coeffOut), 0.1);  // positive or negative CC
}

void testXCorrTrace(const GenericRecordCPtr &tr1,
                    double timeShift,
                    double constant,
                    double scaler)
{
  GenericRecordPtr tr2 =
      timeShift >= 0
          ? alterTrace(tr1, 0, timeShift * 2, constant, scaler)
          : alterTrace(tr1, std::abs(timeShift * 2), 0, constant, scaler);
  testXCorrTraces(tr1, tr2, timeShift);
}

void testXCorr(const GenericRecordCPtr &trace, int s)
{
  testXCorrTrace(trace, 0, 0, 1);
  testXCorrTrace(trace, 0, 0, -1);
  testXCorrTrace(trace, 25. / trace->samplingFrequency(), -1 * s, 100 * s);
  testXCorrTrace(trace, 25. / trace->samplingFrequency(), 1 * s, 0.01 * s);
  testXCorrTrace(trace, 10. / trace->samplingFrequency(), -10 * s, 0.1 * s);
  testXCorrTrace(trace, 10. / trace->samplingFrequency(), 10 * s, 10 * s);
}

void testReampling(const vector<GenericRecordCPtr> &traces)
{
  // traces vector contains identical traces, but different sampling rates
  for (size_t i = 0; i < traces.size() - 1; i++)
  {
    for (size_t j = i + 1; j < traces.size(); j++)
    {
      GenericRecordPtr tr1(new GenericRecord(*traces[i]));
      GenericRecordPtr tr2(new GenericRecord(*traces[j]));
      HDD::Waveform::resample(*tr1, 1000);
      HDD::Waveform::resample(*tr2, 1000);
      testXCorrTraces(tr1, tr2, 0);

      // HDD::Waveform::writeTrace(tr1, strf("testReamplingTr1_%s_%.f-%.fHz.mseed",
      //                           tr1->networkCode().c_str(),
      //                           traces[i]->samplingFrequency(),
      //                           tr1->samplingFrequency()));
      // HDD::Waveform::writeTrace(tr2, strf("testReamplingTr2_%s_%.f-%.fHz.mseed",
      //                           tr2->networkCode().c_str(),
      //                           traces[j]->samplingFrequency(),
      //                           tr2->samplingFrequency()));

      tr1 = new GenericRecord(*traces[i]);
      tr2 = new GenericRecord(*traces[j]);
      HDD::Waveform::resample(*tr1, 60);
      HDD::Waveform::resample(*tr2, 60);
      testXCorrTraces(tr1, tr2, 0);

      // HDD::Waveform::writeTrace(tr1, strf("testReamplingTr1_%s_%.f-%.fHz.mseed",
      //                           tr1->networkCode().c_str(),
      //                           traces[i]->samplingFrequency(),
      //                           tr1->samplingFrequency()));
      // HDD::Waveform::writeTrace(tr2, strf("testReamplingTr2_%s_%.f-%.fHz.mseed",
      //                           tr2->networkCode().c_str(),
      //                           traces[j]->samplingFrequency(),
      //                           tr2->samplingFrequency()));

      tr1 = new GenericRecord(*traces[i]);
      tr2 = new GenericRecord(*traces[j]);
      HDD::Waveform::resample(*tr1, tr2->samplingFrequency());
      testXCorrTraces(tr1, tr2, 0);

      // HDD::Waveform::writeTrace(tr1, strf("testReamplingTr1_%s_%.f-%.fHz.mseed",
      //                           tr1->networkCode().c_str(),
      //                           traces[i]->samplingFrequency(),
      //                           tr1->samplingFrequency()));

      tr1 = new GenericRecord(*traces[i]);
      tr2 = new GenericRecord(*traces[j]);
      HDD::Waveform::resample(*tr2, tr1->samplingFrequency());
      testXCorrTraces(tr1, tr2, 0);

      // HDD::Waveform::writeTrace(tr2, strf("testReamplingTr2_%s_%.f-%.fHz.mseed",
      //                           tr2->networkCode().c_str(),
      //                           traces[j]->samplingFrequency(),
      //                           tr2->samplingFrequency()));
    }
  }
}

void testFiltering(const vector<GenericRecordCPtr> &traces)
{
  const string filterStr = "ITAPER(1)>>BW_HLP(1,1,20)";

  // traces vector contains identical traces, but different sampling rates
  for (size_t i = 0; i < traces.size() - 1; i++)
  {
    for (size_t j = i + 1; j < traces.size(); j++)
    {
      GenericRecordPtr tr1 = alterTrace(traces[i], 0, 0, 10, 1);
      GenericRecordPtr tr2 = alterTrace(traces[j], 0, 0, -10, 1);
      HDD::Waveform::filter(*tr1, true, filterStr, 500);
      HDD::Waveform::filter(*tr2, true, filterStr, 500);
      testXCorrTraces(tr1, tr2, 0);

      // HDD::Waveform::writeTrace(tr1, strf("testFilteringTr1_%s_%.f-%.fHz.mseed",
      //                           tr1->networkCode().c_str(),
      //                           traces[i]->samplingFrequency(),
      //                           tr1->samplingFrequency()));
      // HDD::Waveform::writeTrace(tr2, strf("testFilteringTr2_%s_%.f-%.fHz.mseed",
      //                           tr2->networkCode().c_str(),
      //                           traces[j]->samplingFrequency(),
      //                           tr2->samplingFrequency()));

      tr1 = alterTrace(traces[i], 0, 0, 10, 1);
      tr2 = alterTrace(traces[j], 0, 0, -10, 1); 
      HDD::Waveform::filter(*tr1, true, filterStr, 50);
      HDD::Waveform::filter(*tr2, true, filterStr, 50);
      testXCorrTraces(tr1, tr2, 0);

      // HDD::Waveform::writeTrace(tr1, strf("testFilteringTr1_%s_%.f-%.fHz.mseed",
      //                           tr1->networkCode().c_str(),
      //                           traces[i]->samplingFrequency(),
      //                           tr1->samplingFrequency()));
      // HDD::Waveform::writeTrace(tr2, strf("testFilteringTr2_%s_%.f-%.fHz.mseed",
      //                           tr2->networkCode().c_str(),
      //                           traces[j]->samplingFrequency(),
      //                           tr2->samplingFrequency()));
    }
  }
}

void testSnrSynthetic(const GenericRecordCPtr &trace)
{
  double snr = HDD::Waveform::computeSnr(
      trace, trace->startTime() + TimeSpan(trace->timeWindow().length() / 2),
      -0.5, 0, 0, 0.5);
  BOOST_CHECK(snr > 2);

  snr = HDD::Waveform::computeSnr(
      trace, trace->startTime() + TimeSpan(trace->timeWindow().length() / 2), 0,
      0.5, 0.5, 1);
  BOOST_CHECK(snr < 2);
}

vector<GenericRecordCPtr> realTraces = {
    HDD::Waveform::readTrace("./data/waveform/xcorr1.mseed"),
    HDD::Waveform::readTrace("./data/waveform/xcorr2.mseed"),
    HDD::Waveform::readTrace("./data/waveform/xcorr3.mseed"),
    HDD::Waveform::readTrace("./data/waveform/xcorr4.mseed"),
    HDD::Waveform::readTrace("./data/waveform/xcorr5.mseed"),
    HDD::Waveform::readTrace("./data/waveform/xcorr6.mseed")};

vector<GenericRecordCPtr> badSnrTraces = {
    HDD::Waveform::readTrace("./data/waveform/snr1.mseed"),
    HDD::Waveform::readTrace("./data/waveform/snr2.mseed"),
    HDD::Waveform::readTrace("./data/waveform/snr3.mseed"),
    HDD::Waveform::readTrace("./data/waveform/snr4.mseed"),
    HDD::Waveform::readTrace("./data/waveform/snr5.mseed"),
    HDD::Waveform::readTrace("./data/waveform/snr6.mseed"),
    HDD::Waveform::readTrace("./data/waveform/snr7.mseed")};

vector<GenericRecordCPtr> synthetic1Traces = {
    buildSyntheticTrace1(160), buildSyntheticTrace1(158),
    buildSyntheticTrace1(80), buildSyntheticTrace1(83)};

vector<GenericRecordCPtr> synthetic2Traces = {
    buildSyntheticTrace2(300), buildSyntheticTrace2(280),
    buildSyntheticTrace2(150), buildSyntheticTrace2(143)};

vector<GenericRecordCPtr> synthetic3Traces = {
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
  double snr = HDD::Waveform::computeSnr(
      trace, trace->startTime() + TimeSpan(trace->timeWindow().length() / 2),
      -3, -0.350, -0.350, 0.350);
  BOOST_CHECK(snr < 2);
}

BOOST_DATA_TEST_CASE(test_snr_ok, bdata::make(realTraces), trace)
{
  double snr = HDD::Waveform::computeSnr(
      trace, trace->startTime() + TimeSpan(trace->timeWindow().length() / 2),
      -3, -0.350, -0.350, 0.350);
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

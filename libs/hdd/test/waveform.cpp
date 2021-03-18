#define SEISCOMP_TEST_MODULE hdd
#include <seiscomp/unittest/unittests.h>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>

#include "waveform.h"
#include <boost/filesystem.hpp>
#include <cstring>
#include <seiscomp3/core/genericrecord.h>
#include <seiscomp3/core/strings.h>
#include <seiscomp3/core/typedarray.h>

using namespace std;
using namespace Seiscomp;
using namespace Seiscomp::Core;
using Seiscomp::Core::stringify;
using namespace HDD::Waveform;
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
    value *=
        (i <= tr->samplingFrequency() / 2) ? i : tr->samplingFrequency() - i;
    samples->set(samples->size() / 3 + i, value);
  }
  tr->setData(samples);
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
  return tr;
}

GenericRecordPtr buildSyntheticTrace3(double samplingFrequency)
{
  GenericRecordPtr tr  = new GenericRecord("N3", "ST1", "", "GHZ",
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
  return tr;
}

void scaleTrace(GenericRecordPtr &tr, double constant, double scaler)
{
  DoubleArray *samples = DoubleArray::Cast(tr->data());
  for (int i = 0; i < samples->size(); ++i)
  {
    double value = (*samples)[i];
    value += constant;
    value *= scaler;
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

struct GlobalFixture
{
  static vector<GenericRecordCPtr> realTraces;
  static vector<GenericRecordCPtr> badSnrTraces;
  static vector<GenericRecordCPtr> synthetic1Traces;
  static vector<GenericRecordCPtr> synthetic2Traces;
  static vector<GenericRecordCPtr> synthetic3Traces;
  static const size_t realTracesNum;
  static const size_t badSnrTracesNum;
  static const size_t synthetic1TracesNum;
  static const size_t synthetic2TracesNum;
  static const size_t synthetic3TracesNum;

  GlobalFixture()
  {
    for (size_t i = 1; i <= realTracesNum; i++)
    {
      GenericRecordPtr tr = readTrace(stringify("./waveform/xcorr%d.mseed", i));
      BOOST_REQUIRE(tr);
      realTraces.push_back(tr);
    }

    for (size_t i = 1; i <= badSnrTracesNum; i++)
    {
      GenericRecordPtr tr = readTrace(stringify("./waveform/snr%d.mseed", i));
      BOOST_REQUIRE(tr);
      badSnrTraces.push_back(tr);
    }

    GenericRecordPtr tr;
    tr = buildSyntheticTrace1(160);
    BOOST_REQUIRE(tr);
    synthetic1Traces.push_back(tr);
    tr = buildSyntheticTrace1(158);
    BOOST_REQUIRE(tr);
    synthetic1Traces.push_back(tr);
    tr = buildSyntheticTrace1(80);
    BOOST_REQUIRE(tr);
    synthetic1Traces.push_back(tr);
    tr = buildSyntheticTrace1(83);
    BOOST_REQUIRE(tr);
    synthetic1Traces.push_back(tr);

    tr = buildSyntheticTrace2(300);
    BOOST_REQUIRE(tr);
    synthetic2Traces.push_back(tr);
    tr = buildSyntheticTrace2(280);
    BOOST_REQUIRE(tr);
    synthetic2Traces.push_back(tr);
    tr = buildSyntheticTrace2(150);
    BOOST_REQUIRE(tr);
    synthetic2Traces.push_back(tr);
    tr = buildSyntheticTrace2(143);
    BOOST_REQUIRE(tr);
    synthetic2Traces.push_back(tr);

    tr = buildSyntheticTrace3(80);
    BOOST_REQUIRE(tr);
    synthetic3Traces.push_back(tr);
    tr = buildSyntheticTrace3(78);
    BOOST_REQUIRE(tr);
    synthetic3Traces.push_back(tr);
    tr = buildSyntheticTrace3(40);
    BOOST_REQUIRE(tr);
    synthetic3Traces.push_back(tr);
    tr = buildSyntheticTrace3(43);
    BOOST_REQUIRE(tr);
    synthetic3Traces.push_back(tr);
    /*
    for (size_t i = 0; i < synthetic1Traces.size(); i++)
      writeTrace(synthetic1Traces[i], stringify("synthetic1Trace%d.mseed", i));
    for (size_t i = 0; i < synthetic2Traces.size(); i++)
      writeTrace(synthetic2Traces[i], stringify("synthetic2Trace%d.mseed", i));
    for (size_t i = 0; i < synthetic3Traces.size(); i++)
      writeTrace(synthetic3Traces[i], stringify("synthetic3Trace%d.mseed", i));
    */

    BOOST_REQUIRE(realTraces.size() == realTracesNum);
    BOOST_REQUIRE(badSnrTraces.size() == badSnrTracesNum);
    BOOST_REQUIRE(synthetic1Traces.size() == synthetic1TracesNum);
    BOOST_REQUIRE(synthetic2Traces.size() == synthetic2TracesNum);
    BOOST_REQUIRE(synthetic3Traces.size() == synthetic3TracesNum);
  }
};

vector<GenericRecordCPtr> GlobalFixture::realTraces;
vector<GenericRecordCPtr> GlobalFixture::badSnrTraces;
vector<GenericRecordCPtr> GlobalFixture::synthetic1Traces;
vector<GenericRecordCPtr> GlobalFixture::synthetic2Traces;
vector<GenericRecordCPtr> GlobalFixture::synthetic3Traces;
const size_t GlobalFixture::realTracesNum       = 6;
const size_t GlobalFixture::badSnrTracesNum     = 7;
const size_t GlobalFixture::synthetic1TracesNum = 4;
const size_t GlobalFixture::synthetic2TracesNum = 4;
const size_t GlobalFixture::synthetic3TracesNum = 4;

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
  if (boost::filesystem::exists(filename)) boost::filesystem::remove(filename);
  BOOST_REQUIRE(!boost::filesystem::exists(filename));
  writeTrace(tr, filename);
  GenericRecordPtr cmpTr = readTrace(filename);
  testTracesEqual(tr, cmpTr);
  boost::filesystem::remove(filename);
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
  double delayOut, coeffOut;
  BOOST_CHECK(xcorr(
      tr1, tr2, (tr1->timeWindow().length() - tr2->timeWindow().length()) / 2,
      true, delayOut, coeffOut));
  BOOST_CHECK_SMALL(timeShift - delayOut,
                    2.0 / tr1->samplingFrequency()); // 2 samples tolerance
  BOOST_CHECK_SMALL(1.0 - coeffOut, 0.1);
}

} // namespace

BOOST_TEST_GLOBAL_FIXTURE(GlobalFixture);

BOOST_DATA_TEST_CASE(test_read_write1,
                     bdata::xrange(GlobalFixture::realTracesNum),
                     i)
{
  testReadWriteTrace(GlobalFixture::realTraces[i]);
}

BOOST_DATA_TEST_CASE(test_read_write2,
                     bdata::xrange(GlobalFixture::badSnrTracesNum),
                     i)
{
  testReadWriteTrace(GlobalFixture::badSnrTraces[i]);
}

BOOST_DATA_TEST_CASE(test_read_write3,
                     bdata::xrange(GlobalFixture::synthetic1TracesNum),
                     i)
{
  testReadWriteTrace(GlobalFixture::synthetic1Traces[i]);
}

BOOST_DATA_TEST_CASE(test_read_write4,
                     bdata::xrange(GlobalFixture::synthetic2TracesNum),
                     i)
{
  testReadWriteTrace(GlobalFixture::synthetic2Traces[i]);
}

BOOST_DATA_TEST_CASE(test_read_write5,
                     bdata::xrange(GlobalFixture::synthetic3TracesNum),
                     i)
{
  testReadWriteTrace(GlobalFixture::synthetic3Traces[i]);
}

BOOST_DATA_TEST_CASE(test_xcorr1,
                     bdata::xrange(GlobalFixture::realTracesNum),
                     i)
{
  testXCorrTrace(GlobalFixture::realTraces[i], 0, 0, 1);
  testXCorrTrace(GlobalFixture::realTraces[i], 0.1, 30, 40);
  testXCorrTrace(GlobalFixture::realTraces[i], -0.1, -30, 0.4);
}

BOOST_DATA_TEST_CASE(test_xcorr2,
                     bdata::xrange(GlobalFixture::synthetic1TracesNum),
                     i)
{
  testXCorrTrace(GlobalFixture::synthetic1Traces[i], 0, 0, 1);
  testXCorrTrace(GlobalFixture::synthetic1Traces[i], 0.12, 1.5, 150);
  testXCorrTrace(GlobalFixture::synthetic1Traces[i], -0.15, -1.5, 0.015);
}

BOOST_DATA_TEST_CASE(test_xcorr3,
                     bdata::xrange(GlobalFixture::synthetic2TracesNum),
                     i)
{
  testXCorrTrace(GlobalFixture::synthetic2Traces[i], 0, 0, 1);
  testXCorrTrace(GlobalFixture::synthetic2Traces[i], 0.05, -77, 0.158);
  testXCorrTrace(GlobalFixture::synthetic2Traces[i], -0.05, 77, 158);
}

BOOST_DATA_TEST_CASE(test_xcorr4,
                     bdata::xrange(GlobalFixture::synthetic3TracesNum),
                     i)
{
  testXCorrTrace(GlobalFixture::synthetic3Traces[i], 0, 0, 1);
  testXCorrTrace(GlobalFixture::synthetic3Traces[i], 0.12, -45, 99);
  testXCorrTrace(GlobalFixture::synthetic3Traces[i], -0.12, 14, 0.99);
}

BOOST_AUTO_TEST_CASE(test_resampling) { BOOST_CHECK(1); }

BOOST_AUTO_TEST_CASE(test_snr) { BOOST_CHECK(1); }

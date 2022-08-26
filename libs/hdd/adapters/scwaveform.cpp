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

#include "scwaveform.h"
#include "hdd/log.h"
#include "hdd/utils.h"
#include "scutils.h"

#include <fstream>
#include <iostream>

#include <seiscomp/client/inventory.h>
#include <seiscomp/core/datetime.h>
#include <seiscomp/core/genericrecord.h>
#include <seiscomp/core/recordsequence.h>
#include <seiscomp/core/timewindow.h>
#include <seiscomp/io/recordinput.h>
#include <seiscomp/io/records/mseedrecord.h>
#include <seiscomp/io/recordstream.h>
#include <seiscomp/math/filter.h>

using namespace std;
using namespace Seiscomp;
using Catalog   = HDD::Catalog;
using HDD3Comps = HDD::Waveform::ThreeComponents;
using SC3Comps  = Seiscomp::DataModel::ThreeComponents;

namespace {

using HDD::SCAdapter::fromSC;
using HDD::SCAdapter::toSC;

unique_ptr<HDD::Trace> contiguousRecord(const RecordSequence &seq,
                                        const HDD::TimeWindow &tw,
                                        double minAvailability)
{
  const Core::TimeWindow sctw = toSC(tw);
  if (seq.availability(sctw) < minAvailability)
  {
    string msg =
        HDD::strf("Data availability too low %.2f%%", seq.availability(sctw));
    throw HDD::Exception(msg);
  }

  Seiscomp::GenericRecordPtr rec = seq.contiguousRecord<double>(&sctw, false);

  if (!rec)
  {
    throw HDD::Exception("Internal logic error: cannot create HDD::Trace from "
                         "Seiscomp::Core::GenericRecord");
  }

  const Seiscomp::DoubleArray *data =
      Seiscomp::DoubleArray::ConstCast(rec->data());
  if (!data)
  {
    throw HDD::Exception("Internal logic error: cannot create HDD::Trace from "
                         "Seiscomp::Core::GenericRecord");
  }

  unique_ptr<HDD::Trace> trace(new HDD::Trace(
      rec->networkCode(), rec->stationCode(), rec->locationCode(),
      rec->channelCode(), fromSC(rec->startTime()), rec->samplingFrequency(),
      data->typedData(), data->size()));
  if (!trace->slice(tw))
  {
    string msg =
        HDD::strf("Cannot slice trace from %s length %.2f sec. Trace "
                  "data from %s length %.2f sec, samples %zu sampfreq %f",
                  HDD::UTCClock::toString(tw.startTime()).c_str(),
                  HDD::durToSec(tw.length()),
                  HDD::UTCClock::toString(trace->startTime()).c_str(),
                  HDD::durToSec(trace->timeWindow().length()),
                  trace->sampleCount(), trace->samplingFrequency());
    throw HDD::Exception(msg);
  }

  return trace;
}

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

} // namespace

namespace HDD {
namespace SCAdapter {

unique_ptr<HDD::Trace> WaveformProxy::loadTrace(const HDD::TimeWindow &tw,
                                                const string &networkCode,
                                                const string &stationCode,
                                                const string &locationCode,
                                                const string &channelCode)
{
  const Core::TimeWindow sctw = toSC(tw);

  IO::RecordStreamPtr rs = IO::RecordStream::Open(_recordStreamURL.c_str());
  if (rs == nullptr)
  {
    string msg = "Cannot open RecordStream: " + _recordStreamURL;
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

  return contiguousRecord(seq, tw, minAvailability);
}

void WaveformProxy::loadTraces(
    const unordered_multimap<string, const HDD::TimeWindow> &request,
    const function<void(const string &,
                        const HDD::TimeWindow &,
                        unique_ptr<HDD::Trace>)> &onTraceLoaded,
    const function<void(const string &,
                        const HDD::TimeWindow &,
                        const string &)> &onTraceFailed)
{
  // Prepare the requests in a more convenient format
  struct Request
  {
    string net, sta, loc, ch;
    Core::TimeWindow tw;
  };
  unordered_multimap<string, Request> reqCopy;
  unordered_multimap<string, TimeWindowBuffer> streamBuf;
  for (const auto &kv : request)
  {
    static const std::regex dot("\\.", std::regex::optimize);
    const vector<string> tokens(HDD::splitString(kv.first, dot));
    const Request req{tokens.at(0), tokens.at(1), tokens.at(2), tokens.at(3),
                      toSC(kv.second)};
    reqCopy.emplace(req.net + "." + req.sta, req); // group by net.sta
    streamBuf.emplace(kv.first, TimeWindowBuffer(req.tw, tolerance));
  }

  //
  // Split all the requests in multiple connections to RecordStream.
  // The reason of this split is that RecordStream can only handle one single
  // time window per stream. Also, in case of seedlink, it accepts only one
  // location and channel code per station
  //
  IO::RecordStreamPtr rs;

  while (!reqCopy.empty())
  {
    rs = IO::RecordStream::Open(_recordStreamURL.c_str());
    if (rs == nullptr)
    {
      HDD::logError("Cannot open RecordStream: %s", _recordStreamURL.c_str());
      break;
    }

    // For seedlink: we grouped the requests by net.sta because
    // no multiple requests for the same net.sta are allowed in the same
    // connection
    for (auto it = reqCopy.begin(), end = reqCopy.end();
         it != end;) // loop by net.sta
    {
      const Request req = it->second;
      Core::TimeWindow contiguousRequest;
      auto eqlrng = reqCopy.equal_range(it->first);        // net.sta
      for (auto it2 = eqlrng.first; it2 != eqlrng.second;) // loop by net.sta
      {
        const Request req2         = it2->second;
        const Core::TimeWindow &tw = req2.tw;
        bool requested             = false;

        //
        // Convert multiple time windows requests close to each others to a
        // single RecordStream request
        //
        if (req.loc == req2.loc && req.ch == req2.ch)
        {
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
        }

        if (requested)
          it2 = reqCopy.erase(it2);
        else
          it2++;
      }
      rs->addStream(req.net, req.sta, req.loc, req.ch,
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
    const string &streamID    = kv.first;
    TimeWindowBuffer &seq     = kv.second;
    const HDD::TimeWindow &tw = fromSC(seq.timeWindowToStore());

    try
    {
      unique_ptr<HDD::Trace> trace = contiguousRecord(seq, tw, minAvailability);
      onTraceLoaded(streamID, tw, std::move(trace));
    }
    catch (exception &e)
    {
      onTraceFailed(streamID, tw, e.what());
    }
  }
}

void WaveformProxy::getComponentsInfo(const Catalog::Phase &ph,
                                      HDD3Comps &components)
{
  const Core::Time sctime = toSC(ph.time);
  const string channelCodeRoot =
      HDD::Waveform::getBandAndInstrumentCodes(ph.channelCode);

  Seiscomp::DataModel::Inventory *inv =
      Seiscomp::Client::Inventory::Instance()->inventory();
  if (!inv)
  {
    throw HDD::Exception(
        "Unable to fetch components information: inventory not available");
  }

  Seiscomp::DataModel::InventoryError error;
  Seiscomp::DataModel::SensorLocation *loc =
      Seiscomp::DataModel::getSensorLocation(
          inv, ph.networkCode, ph.stationCode, ph.locationCode, sctime, &error);

  if (!loc)
  {
    string msg = HDD::strf(
        "Unable to fetch SensorLocation information from inventory: %s",
        error.toString());
    throw HDD::Exception(msg);
  }

  SC3Comps tc;
  getThreeComponents(tc, loc, channelCodeRoot.c_str(), sctime);

  if (tc.comps[SC3Comps::Vertical] && tc.comps[SC3Comps::FirstHorizontal] &&
      tc.comps[SC3Comps::SecondHorizontal])
  {
    components.names[HDD3Comps::Vertical] =
        tc.comps[SC3Comps::Vertical]->code();
    components.names[HDD3Comps::FirstHorizontal] =
        tc.comps[SC3Comps::FirstHorizontal]->code();
    components.names[HDD3Comps::SecondHorizontal] =
        tc.comps[SC3Comps::SecondHorizontal]->code();

    components.gain[HDD3Comps::Vertical] = tc.comps[SC3Comps::Vertical]->gain();
    components.gain[HDD3Comps::FirstHorizontal] =
        tc.comps[SC3Comps::FirstHorizontal]->gain();
    components.gain[HDD3Comps::SecondHorizontal] =
        tc.comps[SC3Comps::SecondHorizontal]->gain();

    components.dip[HDD3Comps::Vertical] = tc.comps[SC3Comps::Vertical]->dip();
    components.dip[HDD3Comps::FirstHorizontal] =
        tc.comps[SC3Comps::FirstHorizontal]->dip();
    components.dip[HDD3Comps::SecondHorizontal] =
        tc.comps[SC3Comps::SecondHorizontal]->dip();

    components.azimuth[HDD3Comps::Vertical] =
        tc.comps[SC3Comps::Vertical]->azimuth();
    components.azimuth[HDD3Comps::FirstHorizontal] =
        tc.comps[SC3Comps::FirstHorizontal]->azimuth();
    components.azimuth[HDD3Comps::SecondHorizontal] =
        tc.comps[SC3Comps::SecondHorizontal]->azimuth();
  }
  else
  {
    throw HDD::Exception("Sensor information found in inventory, but it "
                         "doesn't have three components");
  }
}

void WaveformProxy::filter(HDD::Trace &trace, const string &filterStr)
{
  double *data           = trace.data();
  const size_t data_size = trace.sampleCount();

  string filterError;
  auto filter =
      Math::Filtering::InPlaceFilter<double>::Create(filterStr, &filterError);
  if (!filter)
  {
    string msg = HDD::strf("Filter creation failed %s: %s", filterStr.c_str(),
                           filterError.c_str());
    throw HDD::Exception(msg);
  }
  filter->setSamplingFrequency(trace.samplingFrequency());
  filter->apply(data_size, data);
  delete filter;
}

void WaveformProxy::writeTrace(const HDD::Trace &trace, const string &file)
{
  ofstream ofs(file);

  GenericRecord gr(trace.networkCode(), trace.stationCode(),
                   trace.locationCode(), trace.channelCode(),
                   toSC(trace.startTime()), trace.samplingFrequency());
  gr.setData(trace.sampleCount(), trace.data(), Array::DOUBLE);

  IO::MSeedRecord msRec(gr);
  int reclen = msRec.data()->size() * msRec.data()->elementSize() + 64;
  reclen =
      ::nextPowerOf2<int>(reclen, 128 /*MINRECLEN*/, 1048576 /*MAXRECLEN*/);
  if (reclen > 0)
  {
    msRec.setOutputRecordLength(reclen);
    msRec.write(ofs);
  }
}

unique_ptr<HDD::Trace> WaveformProxy::readTrace(const string &file)
{
  ifstream ifs(file);
  IO::MSeedRecord msRec(Array::DOUBLE, Record::Hint::DATA_ONLY);
  msRec.read(ifs);

  const Seiscomp::DoubleArray *data =
      Seiscomp::DoubleArray::ConstCast(msRec.data());
  if (!data)
  {
    throw HDD::Exception("Internal logic error: cannot create HDD::Trace from "
                         "Seiscomp::Core::GenericRecord");
  }
  return unique_ptr<HDD::Trace>(new HDD::Trace(
      msRec.networkCode(), msRec.stationCode(), msRec.locationCode(),
      msRec.channelCode(), fromSC(msRec.startTime()), msRec.samplingFrequency(),
      data->typedData(), data->size()));
}
} // namespace SCAdapter
} // namespace HDD

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

#ifndef __HDD_TRACE_H__
#define __HDD_TRACE_H__

#include "timewindow.h"
#include "utctime.h"

#include <cmath>
#include <string>
#include <vector>

namespace HDD {

template <typename T_DATA, typename T_TIME> class GenericTrace
{
private:
  std::string _net;
  std::string _sta;
  std::string _loc;
  std::string _cha;
  T_TIME _stime;
  double _smpfreq;
  std::vector<T_DATA> _data;

public:
  GenericTrace(const std::string &net,
               const std::string &sta,
               const std::string &loc,
               const std::string &cha,
               const T_TIME &stime,
               double smpfreq,
               std::vector<T_DATA> &&data) // take ownership
      : _net(net), _sta(sta), _loc(loc), _cha(cha), _stime(stime),
        _smpfreq(smpfreq), _data(std::move(data))
  {}

  GenericTrace(const std::string &net,
               const std::string &sta,
               const std::string &loc,
               const std::string &cha,
               const T_TIME &stime,
               double smpfreq,
               const T_DATA *data,
               size_t dataSize) // copy
      : GenericTrace(net,
                     sta,
                     loc,
                     cha,
                     stime,
                     smpfreq,
                     std::vector<T_DATA>{data, data + dataSize})
  {}

  GenericTrace(const std::string &net,
               const std::string &sta,
               const std::string &loc,
               const std::string &cha,
               const T_TIME &stime,
               double smpfreq)
      : GenericTrace(net, sta, loc, cha, stime, smpfreq, std::vector<T_DATA>())
  {}

  ~GenericTrace() = default;

  GenericTrace(const GenericTrace &other) = default;
  GenericTrace &operator=(const GenericTrace &other) = default;

  GenericTrace(GenericTrace &&other) = default;
  GenericTrace &operator=(GenericTrace &&other) = default;

  std::string networkCode() const { return _net; }
  void setNetworkCode(const std::string &net) { _net = net; }

  std::string stationCode() const { return _sta; }
  void setStationCode(const std::string &sta) const { _sta = sta; }

  std::string locationCode() const { return _loc; }
  void setLocationCode(const std::string &loc) { _loc = loc; }

  std::string channelCode() const { return _cha; }
  void setChannelCode(const std::string &cha) { _cha = cha; }

  const T_TIME &startTime() const { return _stime; }
  void setStartTime(const T_TIME &time) { _stime = time; }

  T_TIME endTime() const { return _stime + secToDur(sampleCount() / _smpfreq); }
  void setEndTime(const T_TIME &time) { _stime += time - endTime(); }

  TimeWindow timeWindow() const { return TimeWindow(_stime, endTime()); }

  size_t sampleCount() const { return _data.size(); }

  size_t elementSize() const { return sizeof(T_DATA); }

  double samplingFrequency() const { return _smpfreq; }

  void setSamplingFrequency(double freq) { _smpfreq = freq; }

  std::string streamID() const
  {
    return _net + "." + _sta + "." + _loc + "." + _cha;
  }

  const T_DATA *data() const { return _data.data(); }

  T_DATA *data() { return _data.data(); }

  // take ownership
  void setData(std::vector<T_DATA> &&data) { _data = std::move(data); }

  // copy
  void setData(T_DATA *data, size_t dataSize)
  {
    setData(std::vector<T_DATA>{data, data + dataSize});
  }

  bool slice(const TimeWindow &tw)
  {
    if (timeWindow() == tw) return true;

    if (tw.startTime() < _stime) return false;
    if (tw.endTime() > endTime()) return false;

    size_t startOfs = std::floor(durToSec(tw.startTime() - _stime) * _smpfreq);
    size_t endOfs   = std::ceil(durToSec(tw.endTime() - _stime) * _smpfreq) - 1;

    if (startOfs >= sampleCount()) return false; // overflow
    if (endOfs >= sampleCount()) return false;

    std::vector<T_DATA> sliced(_data.begin() + startOfs,
                               _data.begin() + endOfs);
    setStartTime(startTime() + secToDur(startOfs / _smpfreq));
    setData(std::move(sliced));
    return true;
  }
};

using Trace = GenericTrace<double, UTCTime>;

} // namespace HDD

#endif

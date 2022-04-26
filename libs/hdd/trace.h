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

template <typename T_DATA, typename T_TIME, typename T_WIN> class GenericTrace
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

  T_TIME endTime() const
  {
    return sampleCount() == 0
               ? _stime
               : (_stime + secToDur((sampleCount() - 1) / _smpfreq));
  }
  void setEndTime(const T_TIME &time) { _stime += time - endTime(); }

  T_WIN timeWindow() const { return T_WIN(_stime, endTime()); }

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

  double index(const T_TIME time) const
  {
    return durToSec(time - _stime) * _smpfreq;
  }

  bool slice(const T_WIN &tw)
  {
    const T_WIN _tw = timeWindow();
    if (_tw == tw) return true;

    if (!_tw.contains(tw)) return false;

    double startOfs = std::floor(index(tw.startTime()));
    double endOfs   = std::ceil(index(tw.endTime()));

    if (startOfs < 0) return false;
    if (endOfs >= sampleCount()) return false;

    std::vector<T_DATA> sliced(_data.begin() + startOfs,
                               _data.begin() + endOfs + 1);
    setStartTime(startTime() + secToDur(startOfs / _smpfreq));
    setData(std::move(sliced));
    return true;
  }
};

using Trace = GenericTrace<double, UTCTime, TimeWindow>;

} // namespace HDD

#endif

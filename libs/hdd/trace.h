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

  GenericTrace(const GenericTrace &other)            = default;
  GenericTrace &operator=(const GenericTrace &other) = default;

  GenericTrace(GenericTrace &&other)            = default;
  GenericTrace &operator=(GenericTrace &&other) = default;

  std::string networkCode() const { return _net; }
  void setNetworkCode(const std::string &net) { _net = net; }

  std::string stationCode() const { return _sta; }
  void setStationCode(const std::string &sta) { _sta = sta; }

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

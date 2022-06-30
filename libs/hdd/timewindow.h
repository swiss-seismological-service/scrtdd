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

#ifndef __HDD_TIMEWINDOW_H__
#define __HDD_TIMEWINDOW_H__

#include "utctime.h"

namespace HDD {

template <class T, class Duration = typename T::duration>
class GenericTimeWindow
{
private:
  T _startTime, _endTime;

public:
  GenericTimeWindow() = default;
  GenericTimeWindow(const T &startTime, Duration length)
      : GenericTimeWindow(startTime, startTime + length)
  {}
  GenericTimeWindow(const T &startTime, const T &endTime)
      : _startTime(startTime)
  {
    setEndTime(endTime);
  }

  ~GenericTimeWindow() = default;

  GenericTimeWindow(const GenericTimeWindow &other) = default;
  GenericTimeWindow &operator=(const GenericTimeWindow &other) = default;

  GenericTimeWindow(GenericTimeWindow &&other) = default;
  GenericTimeWindow &operator=(GenericTimeWindow &&other) = default;

  bool operator==(const GenericTimeWindow &other) const
  {
    return _startTime == other._startTime && _endTime == other._endTime;
  }

  bool operator!=(const GenericTimeWindow &other) const
  {
    return _startTime != other._startTime || _endTime != other._endTime;
  }

  const T &startTime() const { return _startTime; }

  const T &endTime() const { return _endTime; }

  Duration length() const { return _endTime - _startTime; }

  void set(const T &startTime, const T &endTime)
  {
    setStartTime(startTime);
    setEndTime(endTime);
  }

  void setStartTime(const T &t)
  {
    _startTime = t;
    if (_endTime < _startTime)
    {
      _endTime = _startTime;
    }
  }

  void setEndTime(const T &t)
  {
    _endTime = t;
    if (_endTime < _startTime) _endTime = _startTime;
  }

  bool empty() const { return _startTime == _endTime; }

  bool contains(const T &t) const { return t >= _startTime && t <= _endTime; }

  bool contains(const GenericTimeWindow &other) const
  {
    return other._startTime >= _startTime && other._endTime <= _endTime;
  }

  bool overlaps(const GenericTimeWindow &other) const
  {
    return !overlap(other).empty();
  }

  bool contiguous(const GenericTimeWindow &other) const
  {
    return other._startTime == _endTime || other._endTime == _startTime;
  }

  GenericTimeWindow trim(Duration amountAtStart, Duration amountAtEnd) const
  {
    return GenericTimeWindow(_startTime + amountAtStart,
                             _endTime - amountAtEnd);
  }

  GenericTimeWindow extend(Duration amountAtStart, Duration amountAtEnd) const
  {
    return GenericTimeWindow(_startTime - amountAtStart,
                             _endTime + amountAtEnd);
  }

  GenericTimeWindow merge(const GenericTimeWindow &other) const
  {
    GenericTimeWindow tw(*this);
    if (other.startTime() < tw.startTime()) tw.setStartTime(other.startTime());
    if (other.endTime() > tw.endTime()) tw.setEndTime(other.endTime());
    return tw;
  }

  GenericTimeWindow overlap(const GenericTimeWindow &other) const
  {
    if (contains(other)) return other;
    if (other.contains(*this)) return *this;
    if (contains(other._startTime))
      return GenericTimeWindow(other._startTime, _endTime);
    if (contains(other._endTime))
      return GenericTimeWindow(_startTime, other._endTime);
    return GenericTimeWindow();
  }
};

using TimeWindow = GenericTimeWindow<UTCTime>;

} // namespace HDD

#endif

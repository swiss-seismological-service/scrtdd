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

#ifndef __HDD_UTCTIME_H__
#define __HDD_UTCTIME_H__

#include <chrono>
#include <string>

namespace HDD {

struct UTCClock
{
  typedef std::chrono::microseconds duration;
  typedef duration::rep rep;
  typedef duration::period period;
  typedef std::chrono::time_point<UTCClock, duration> time_point;
  static const bool is_steady = true;

  static time_point now();

  static time_point fromDate(int year  = 0,
                             int month = 0,
                             int day   = 0,
                             int hour  = 0,
                             int min   = 0,
                             int sec   = 0,
                             int usec  = 0);

  static void toDate(const time_point &tp,
                     int &year,
                     int &month,
                     int &day,
                     int &hour,
                     int &min,
                     int &sec,
                     int &usec);

  static time_point fromString(const std::string &s);

  static std::string toString(const time_point &tp);
};

using UTCTime = UTCClock::time_point;

inline UTCTime::duration secToDur(double sec)
{
  return std::chrono::duration_cast<UTCTime::duration>(
      std::chrono::duration<double>(sec));
}

inline double durToSec(UTCTime::duration d)
{
  return std::chrono::duration<double>(d).count();
}

} // namespace HDD

#endif

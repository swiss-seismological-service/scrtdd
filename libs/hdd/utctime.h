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

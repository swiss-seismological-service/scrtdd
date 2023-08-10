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

#include "utctime.h"
#include "utils.h"
#include <ctime>
#include <regex>

namespace chrono = std::chrono;
using chrono::duration_cast;
using chrono::time_point_cast;
using HDD::UTCClock;

namespace {

std::time_t to_time_t(const UTCClock::time_point &tp) noexcept
{
  return std::time_t(
      duration_cast<chrono::seconds>(tp.time_since_epoch()).count());
}

UTCClock::time_point from_time_t(std::time_t tt) noexcept
{
  return time_point_cast<UTCClock::duration>(
      chrono::time_point<UTCClock, chrono::seconds>(chrono::seconds(tt)));
}

} // namespace

namespace HDD {

UTCClock::time_point UTCClock::now()
{
  std::time_t now =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  return from_time_t(now);
}

UTCClock::time_point UTCClock::fromDate(
    int year, int month, int day, int hour, int min, int sec, int usec)
{
  std::tm tm     = {};
  tm.tm_year     = year - 1900;
  tm.tm_mon      = month - 1;
  tm.tm_mday     = day;
  tm.tm_hour     = hour;
  tm.tm_min      = min;
  tm.tm_sec      = sec;
  tm.tm_isdst    = -1;
  std::time_t tt = timegm(&tm);
  return from_time_t(tt) + chrono::microseconds(usec);
}

void UTCClock::toDate(const UTCClock::time_point &tp,
                      int &year,
                      int &month,
                      int &day,
                      int &hour,
                      int &min,
                      int &sec,
                      int &usec)
{
  std::time_t tt = to_time_t(tp);
  std::tm tm;
  gmtime_r(&tt, &tm);
  year  = tm.tm_year + 1900;
  month = tm.tm_mon + 1;
  day   = tm.tm_mday;
  hour  = tm.tm_hour;
  min   = tm.tm_min;
  chrono::microseconds leftover =
      tp - from_time_t(tt) + chrono::seconds(tm.tm_sec);
  sec  = duration_cast<chrono::seconds>(leftover).count();
  usec = (leftover - chrono::seconds(sec)).count();
}

UTCClock::time_point UTCClock::fromString(const std::string &s)
{
  static const std::regex re(
      R"((\d\d\d\d)-(\d\d)-(\d\d)T(\d\d):(\d\d):(\d\d)\.(\d+)Z)",
      std::regex::optimize);

  std::smatch m;
  if (!std::regex_match(s, m, re)) throw Exception("Invalid UTC string: " + s);

  int year  = std::stoi(m.str(1));
  int month = std::stoi(m.str(2));
  int day   = std::stoi(m.str(3));
  int hour  = std::stoi(m.str(4));
  int min   = std::stoi(m.str(5));
  int sec   = std::stoi(m.str(6));
  int usec  = std::stod(m.str(7));
  // This mess below is because we cannot read the usec as a floating point
  // number but we must read it as an integer otherwise we lose precision
  size_t numberOfDigits = m.str(7).size();
  while (numberOfDigits < 6)
  {
    usec *= 10;
    numberOfDigits++;
  }
  while (numberOfDigits > 6)
  {
    usec /= 10;
    numberOfDigits--;
  }

  if (month < 1 || month > 12) throw Exception("Invalid UTC string: " + s);
  if (day < 1 || day > 31) throw Exception("Invalid UTC string: " + s);
  if (hour < 0 || hour > 23) throw Exception("Invalid UTC string: " + s);
  if (min < 0 || min > 59) throw Exception("Invalid UTC string: " + s);
  if (sec < 0 || sec > 59) throw Exception("Invalid UTC string: " + s);

  return fromDate(year, month, day, hour, min, sec, usec);
}

std::string UTCClock::toString(const UTCClock::time_point &tp)
{
  int year, month, day, hour, min, sec, usec;
  toDate(tp, year, month, day, hour, min, sec, usec);
  std::string usecStr;
  if (usec == 0)
    usecStr = "0000";
  else
  {
    usecStr               = strf("%06d", usec);
    size_t numberOfDigits = usecStr.size();
    while (numberOfDigits > 0 && usecStr[numberOfDigits - 1] == '0')
      --numberOfDigits;
    usecStr = usecStr.substr(0, numberOfDigits);
  }
  return strf("%04d-%02d-%02dT%02d:%02d:%02d.%sZ", year, month, day, hour, min,
              sec, usecStr.c_str());
}

} // namespace HDD

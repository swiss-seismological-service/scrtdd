/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU Affero General Public License as published*
 * by the Free Software Foundation, either version 3 of the License, or    *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This program is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#ifndef __HDD_SCCONVERSION_H__
#define __HDD_SCCONVERSION_H__

#include "hdd/timewindow.h"
#include "hdd/utctime.h"

#include <seiscomp/core/datetime.h>
#include <seiscomp/core/genericrecord.h>
#include <seiscomp/core/timewindow.h>
#include <seiscomp/core/typedarray.h>

namespace HDDSCAdapter {

inline Seiscomp::Core::Time toSC(const HDD::UTCTime &t)
{
  return Seiscomp::Core::Time(HDD::durToSec(t.time_since_epoch()));
}

inline HDD::UTCTime fromSC(const Seiscomp::Core::Time &t)
{
  return HDD::UTCTime() + HDD::secToDur(t.length());
}

inline Seiscomp::Core::TimeWindow toSC(const HDD::TimeWindow &tw)
{
  return Seiscomp::Core::TimeWindow(toSC(tw.startTime()), toSC(tw.endTime()));
}

inline HDD::TimeWindow fromSC(const Seiscomp::Core::TimeWindow &tw)
{
  return HDD::TimeWindow(fromSC(tw.startTime()), fromSC(tw.endTime()));
}

} // namespace HDDSCAdapter

#endif

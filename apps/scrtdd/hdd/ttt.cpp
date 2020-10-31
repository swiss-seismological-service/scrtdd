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

#include "ttt.h"
#include "utils.h"

#include <seiscomp3/core/strings.h>
#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/math.h>
#include <sstream>
#include <stdexcept>

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp3/logging/log.h>

using namespace std;
using Seiscomp::Core::stringify;

namespace Seiscomp {
namespace HDD {

TravelTimeTable::TravelTimeTable(std::string type, std::string model)
{
  _ttt = TravelTimeTableInterface::Create(type.c_str());
  _ttt->setModel(model.c_str());
}

void TravelTimeTable::compute(double eventLat,
                              double eventLon,
                              double eventDepth,
                              double stationLat,
                              double stationLon,
                              double stationElevation,
                              const std::string &phaseType,
                              double &travelTime,
                              double &takeOffAngle,
                              double &velocityAtSrc)
{
  // when takeOffAngle/velocityAtSrc are not provided by the ttt then
  // the solver will use straight ray path approximation
  // Note: deg2rad(tt.takeoff) doesn't seem to be correct
  travelTime = takeOffAngle = velocityAtSrc = 0;

  double depth  = eventDepth > 0 ? eventDepth : 0;
  TravelTime tt = _ttt->compute(phaseType.c_str(), eventLat, eventLon, depth,
                                stationLat, stationLon, stationElevation);
  travelTime    = tt.time;
}

} // namespace HDD
} // namespace Seiscomp

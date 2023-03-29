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

#include "cvttt.h"
#include "utils.h"

namespace HDD {

ConstantVelocity::ConstantVelocity(double pVel, double sVel)
    : _pVel(pVel), _sVel(sVel)
{}

void ConstantVelocity::freeResources() {}

double ConstantVelocity::compute(double eventLat,
                                 double eventLon,
                                 double eventDepth,
                                 const Catalog::Station &station,
                                 const std::string &phaseType)
{
  double velocity; // [km/s]
  if (phaseType == "P")
    velocity = _pVel;
  else if (phaseType == "S")
    velocity = _sVel;
  else
    throw Exception("Unknown phase type: " + phaseType);

  // straight ray path since we are in a homogeneous media
  double distance = computeDistance(eventLat, eventLon, eventDepth,
                                    station.latitude, station.longitude,
                                    -(station.elevation / 1000.)); // [km]
  return distance / velocity;                                      // [sec]
}

void ConstantVelocity::compute(double eventLat,
                               double eventLon,
                               double eventDepth,
                               const Catalog::Station &station,
                               const std::string &phaseType,
                               double &travelTime,
                               double &azimuth,
                               double &takeOffAngle,
                               double &velocityAtSrc)
{
  travelTime = compute(eventLat, eventLon, eventDepth, station, phaseType);
  computeApproximatedTakeOffAngles(eventLat, eventLon, eventDepth, station,
                                   phaseType, &azimuth, &takeOffAngle);
  if (phaseType == "P")
    velocityAtSrc = _pVel;
  else if (phaseType == "S")
    velocityAtSrc = _sVel;
  else
    throw Exception("Unknown phase type: " + phaseType);
}

} // namespace HDD

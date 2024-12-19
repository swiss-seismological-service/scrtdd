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

#include "homogeneous.h"
#include "../utils.h"

namespace HDD {
namespace TTT {

Homogeneous::Homogeneous(double pVel, double sVel) : _pVel(pVel), _sVel(sVel) {}

double Homogeneous::compute(double eventLat,
                            double eventLon,
                            double eventDepth,
                            double stationLat,
                            double stationLon,
                            double stationElevation,
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
  double distance =
      computeDistance(eventLat, eventLon, eventDepth, stationLat, stationLon,
                      -(stationElevation / 1000.)); // [km]
  return distance / velocity;                       // [sec]
}

void Homogeneous::compute(double eventLat,
                          double eventLon,
                          double eventDepth,
                          double stationLat,
                          double stationLon,
                          double stationElevation,
                          const std::string &phaseType,
                          double &travelTime,
                          double &azimuth,
                          double &takeOffAngle,
                          double &velocityAtSrc)
{
  travelTime = compute(eventLat, eventLon, eventDepth, stationLat, stationLon,
                       stationElevation, phaseType);

  double hDist =
      computeDistance(eventLat, eventLon, eventDepth, stationLat, stationLon,
                      -(stationElevation / 1000.), &azimuth);

  double vDist = eventDepth + stationElevation / 1000.;
  takeOffAngle = std::asin(vDist / hDist);
  takeOffAngle += degToRad(90); // -90(down):+90(up) -> 0(down):180(up)

  if (phaseType == "P")
    velocityAtSrc = _pVel;
  else if (phaseType == "S")
    velocityAtSrc = _sVel;
  else
    throw Exception("Unknown phase type: " + phaseType);
}
} // namespace TTT

} // namespace HDD

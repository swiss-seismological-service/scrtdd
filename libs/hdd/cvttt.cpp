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

#include "cvttt.h"
#include "utils.h"

namespace HDD {

ConstantVelocity::ConstantVelocity(double pVel, double sVel)
    : _pVel(pVel), _sVel(sVel)
{}

void ConstantVelocity::freeResources() {}

void ConstantVelocity::compute(double eventLat,
                               double eventLon,
                               double eventDepth,
                               const Catalog::Station &station,
                               const std::string &phaseType,
                               double &travelTime)
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
  travelTime      = distance / velocity;                           // [sec]
}

void ConstantVelocity::compute(double eventLat,
                               double eventLon,
                               double eventDepth,
                               const Catalog::Station &station,
                               const std::string &phaseType,
                               double &travelTime,
                               double &takeOffAngleAzim,
                               double &takeOffAngleDip,
                               double &velocityAtSrc)
{
  compute(eventLat, eventLon, eventDepth, station, phaseType, travelTime);
  computeApproximatedTakeOfAngles(eventLat, eventLon, eventDepth, station,
                                  phaseType, &takeOffAngleAzim,
                                  &takeOffAngleDip);
  if (phaseType == "P")
    velocityAtSrc = _pVel;
  else if (phaseType == "S")
    velocityAtSrc = _sVel;
  else
    throw Exception("Unknown phase type: " + phaseType);
}

} // namespace HDD

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

#include "ttt.h"
#include "log.h"
#include "utils.h"

#include <sstream>

using namespace std;

namespace HDD {

void TravelTimeTable::computeApproximatedTakeOfAngles(
    double eventLat,
    double eventLon,
    double eventDepth,
    const Catalog::Station &station,
    const std::string &phaseType,
    double *takeOffAngleAzim,
    double *takeOffAngleDip)
{

  if (takeOffAngleAzim || takeOffAngleDip)
  {
    double azimuth;
    double distance = computeDistance(
        eventLat, eventLon, eventDepth, station.latitude, station.longitude,
        -(station.elevation / 1000.), &azimuth, nullptr);

    if (takeOffAngleDip)
    {
      double VertDist  = eventDepth + station.elevation / 1000.;
      *takeOffAngleDip = std::asin(VertDist / distance);
      *takeOffAngleDip += degToRad(90); // -90(down):+90(up) -> 0(down):180(up)
    }

    if (takeOffAngleAzim)
    {
      *takeOffAngleAzim = degToRad(azimuth);
    }
  }
}

} // namespace HDD

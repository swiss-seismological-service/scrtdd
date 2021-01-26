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
#include "nllttt.h"
#include "scttt.h"
#include "utils.h"
#include <seiscomp3/math/math.h>

#include <seiscomp3/core/strings.h>
#include <sstream>
#include <stdexcept>

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp3/logging/log.h>

using namespace std;
using Seiscomp::Core::stringify;

namespace Seiscomp {
namespace HDD {

TravelTimeTable *TravelTimeTable::create(const std::string &type,
                                         const std::string &model)
{
  TravelTimeTable *ttt = nullptr;

  if (type == "LOCSAT" || type == "libtau")
  {
    ttt = new ScTravelTimeTable(type, model);
  }
  else if (type == "NonLinLoc")
  {
    ttt = new NLL::NllTravelTimeTable(type, model);
  }

  if (ttt == nullptr)
  {
    string msg = stringify("Cannot load travel time table: unknown type %s",
                           type.c_str());
    throw runtime_error(msg.c_str());
  }
  return ttt;
}

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
    double backAzimuth;
    double distance = computeDistance(
        eventLat, eventLon, eventDepth, station.latitude, station.longitude,
        -(station.elevation / 1000.), nullptr, &backAzimuth);

    if (takeOffAngleDip)
    {
      double VertDist  = eventDepth + station.elevation / 1000.;
      *takeOffAngleDip = std::asin(VertDist / distance);
    }

    if (takeOffAngleAzim)
    {
      *takeOffAngleAzim = deg2rad(backAzimuth);
    }
  }
}

} // namespace HDD
} // namespace Seiscomp

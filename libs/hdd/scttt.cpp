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

#include "scttt.h"
#include "utils.h"

#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/math.h>
#include <sstream>

using namespace std;

namespace HDD {

ScTravelTimeTable::ScTravelTimeTable(const std::string &type,
                                     const std::string &model,
                                     double depthVelResolution)
    : TravelTimeTable(type, model), _depthVelResolution(depthVelResolution)
{
  _ttt = TravelTimeTableInterface::Create(type.c_str());
  _ttt->setModel(model.c_str());

  /*
  for (int i = 0; i < 500; i++)
  {
    double vel = velocityAtSource(i * _depthVelResolution, "P");
    logInfo("Velocity model phase P depth %.2f [km] vel %.2f [m/sec]",
                   (i * _depthVelResolution), vel);
  }
  for (int i = 0; i < 500; i++)
  {
    double vel = velocityAtSource(i * _depthVelResolution, "S");
    logInfo("Velocity model phase S depth %.2f [km] vel %.2f [m/sec]",
                   (i * _depthVelResolution), vel);
  }
  */
}

void ScTravelTimeTable::compute(double eventLat,
                                double eventLon,
                                double eventDepth,
                                const Catalog::Station &station,
                                const std::string &phaseType,
                                double &travelTime)
{
  double depth = eventDepth > 0 ? eventDepth : 0;
  TravelTime tt =
      _ttt->compute(phaseType.c_str(), eventLat, eventLon, depth,
                    station.latitude, station.longitude, station.elevation);
  travelTime = tt.time;
}

void ScTravelTimeTable::compute(double eventLat,
                                double eventLon,
                                double eventDepth,
                                const Catalog::Station &station,
                                const std::string &phaseType,
                                double &travelTime,
                                double &takeOffAngleAzim,
                                double &takeOffAngleDip,
                                double &velocityAtSrc)
{
  double depth = eventDepth > 0 ? eventDepth : 0;
  TravelTime tt =
      _ttt->compute(phaseType.c_str(), eventLat, eventLon, depth,
                    station.latitude, station.longitude, station.elevation);
  travelTime = tt.time;
  computeApproximatedTakeOfAngles(eventLat, eventLon, eventDepth, station,
                                  phaseType, &takeOffAngleAzim,
                                  &takeOffAngleDip);
  // tt.takeoff is not computed for LOCSAT
  if (type == "libtau")
  {
    takeOffAngleDip = deg2rad(tt.takeoff);
  }
  velocityAtSrc = velocityAtSource(depth, phaseType);
}

// Since the seiscomp travel time API doesn't offer the velocity at source we
// need to reverse-engineer that information.
double ScTravelTimeTable::velocityAtSource(double eventDepth,
                                           const std::string &phaseType)
{
  if (eventDepth < 0) eventDepth = 0;

  const unsigned bin         = std::floor(eventDepth / _depthVelResolution);
  const double binStartDepth = bin * _depthVelResolution;
  const double binEndDepth   = (bin + 1) * _depthVelResolution;

  //
  // first check if we have alrady computed the velocity for this phase/depth
  //
  auto it1 = _depthVel.find(phaseType);
  if (it1 != _depthVel.end())
  {
    const auto &phaseDepthVel = it1->second;
    auto it2                  = phaseDepthVel.find(bin);
    if (it2 != phaseDepthVel.end()) return it2->second;
  }

  //
  // this is a new phase/depth pair
  //
  double tt1 =
      (binStartDepth == 0)
          ? 0
          : _ttt->compute(phaseType.c_str(), 0, 0, binStartDepth, 0, 0, 0).time;
  double tt2 =
      _ttt->compute(phaseType.c_str(), 0, 0, binEndDepth, 0, 0, 0).time;

  double binVelocity = _depthVelResolution / (tt2 - tt1); // [km/sec]

  if (!std::isfinite(binVelocity))
  {
    throw Exception("Unable to compute velocity at source");
  }

  // store the value
  _depthVel[phaseType][bin] = binVelocity;

  return binVelocity;
}

} // namespace HDD

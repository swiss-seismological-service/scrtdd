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
#include "hdd/log.h"
#include "hdd/utils.h"

using namespace std;
using namespace HDD;

namespace HDD {
namespace SCAdapter {

TravelTimeTable::TravelTimeTable(const std::string &type,
                                 const std::string &model,
                                 double depthVelResolution)
    : _type(type), _model(model), _depthVelResolution(depthVelResolution)
{
  load();
  /*
  for (int i = -10; i < 50; i++)
  {
    double vel = velocityAtSource(0, 0, i * _depthVelResolution, "P");
    logInfo("Velocity model phase P depth %.2f [km] vel %.2f [m/sec]",
            (i * _depthVelResolution), vel);
  }
  for (int i = -10; i < 50; i++)
  {
    double vel = velocityAtSource(0, 0, i * _depthVelResolution, "S");
    logInfo("Velocity model phase S depth %.2f [km] vel %.2f [m/sec]",
            (i * _depthVelResolution), vel);
  }
  */
}

void TravelTimeTable::load()
{
  _ttt = Seiscomp::TravelTimeTableInterface::Create(_type.c_str());
  if (!_ttt || !_ttt->setModel(_model.c_str()))
  {
    throw Exception(strf("Unable to set travel time format %s and model %s",
                         _type.c_str(), _model.c_str()));
  }
}

void TravelTimeTable::freeResources() { _ttt = nullptr; }

void TravelTimeTable::compute(double eventLat,
                              double eventLon,
                              double eventDepth,
                              const Catalog::Station &station,
                              const std::string &phaseType,
                              double &travelTime)
{
  if (!_ttt) load();

  double depth = eventDepth > 0 ? eventDepth : 0;
  Seiscomp::TravelTime tt =
      _ttt->compute(phaseType.c_str(), eventLat, eventLon, depth,
                    station.latitude, station.longitude, station.elevation);
  travelTime = tt.time;
}

void TravelTimeTable::compute(double eventLat,
                              double eventLon,
                              double eventDepth,
                              const Catalog::Station &station,
                              const std::string &phaseType,
                              double &travelTime,
                              double &takeOffAngleAzim,
                              double &takeOffAngleDip,
                              double &velocityAtSrc)
{
  if (!_ttt) load();

  Seiscomp::TravelTime tt =
      _ttt->compute(phaseType.c_str(), eventLat, eventLon, eventDepth,
                    station.latitude, station.longitude, station.elevation);
  travelTime = tt.time;
  computeApproximatedTakeOfAngles(eventLat, eventLon, eventDepth, station,
                                  phaseType, &takeOffAngleAzim,
                                  &takeOffAngleDip);
  // tt.takeoff is not computed for LOCSAT
  if (_type != "LOCSAT")
  {
    takeOffAngleDip = degToRad(tt.takeoff);
  }
  velocityAtSrc = velocityAtSource(eventLat, eventLon, eventDepth, phaseType);
}

// Since the seiscomp travel time API doesn't offer the velocity at source we
// need to reverse-engineer that information.
double TravelTimeTable::velocityAtSource(double eventLat,
                                         double eventLon,
                                         double eventDepth,
                                         const std::string &phaseType)
{
  const int bin              = std::floor(eventDepth / _depthVelResolution);
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
  double tt1 = (binStartDepth == 0)
                   ? 0
                   : _ttt->compute(phaseType.c_str(), eventLat, eventLon,
                                   binStartDepth, eventLat, eventLon, 0)
                         .time;
  double tt2 = (binEndDepth == 0)
                   ? 0
                   : _ttt->compute(phaseType.c_str(), eventLat, eventLon,
                                   binEndDepth, eventLat, eventLon, 0)
                         .time;

  double binVelocity = _depthVelResolution / std::abs(tt2 - tt1); // [km/sec]

  if (!std::isfinite(binVelocity))
  {
    throw Exception("Unable to compute velocity at source");
  }

  // store the value
  _depthVel[phaseType][bin] = binVelocity;

  return binVelocity;
}
} // namespace SCAdapter
} // namespace HDD

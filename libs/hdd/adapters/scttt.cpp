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
#include <seiscomp/math/geo.h>

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

  //
  // Test for computeVelocityAtSource()
  //
  // logInfo("TravelTimeTable %s %s", type.c_str(), model.c_str());
  //
  // for (int i = -10; i < 50; i++)
  //{
  //  try
  //  {
  //    double vel = computeVelocityAtSource(0, 0, i * _depthVelResolution,
  //    "P"); logInfo("Velocity model phase P depth %.2f [km] vel %.2f [m/sec]",
  //            (i * _depthVelResolution), vel);
  //  }
  //  catch (...)
  //  {}
  //}
  // for (int i = -10; i < 50; i++)
  //{
  //  try
  //  {
  //    double vel = computeVelocityAtSource(0, 0, i * _depthVelResolution,
  //    "S"); logInfo("Velocity model phase S depth %.2f [km] vel %.2f [m/sec]",
  //            (i * _depthVelResolution), vel);
  //  }
  //  catch (...)
  //  {}
  //}
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

double TravelTimeTable::compute(double eventLat,
                                double eventLon,
                                double eventDepth,
                                const Catalog::Station &station,
                                const std::string &phaseType)
{
  if (!_ttt) load();

  Seiscomp::TravelTime tt =
      _ttt->compute(phaseType.c_str(), eventLat, eventLon, eventDepth,
                    station.latitude, station.longitude, station.elevation);
  if (tt.time < 0)
  {
    throw Exception("No travel time data available");
  }

  return tt.time;
}

void TravelTimeTable::compute(double eventLat,
                              double eventLon,
                              double eventDepth,
                              const Catalog::Station &station,
                              const std::string &phaseType,
                              double &travelTime,
                              double &azimuth,
                              double &takeOffAngle,
                              double &velocityAtSrc)
{
  if (!_ttt) load();

  Seiscomp::TravelTime tt =
      _ttt->compute(phaseType.c_str(), eventLat, eventLon, eventDepth,
                    station.latitude, station.longitude, station.elevation);
  if (tt.time < 0)
  {
    throw Exception("No travel time data available");
  }

  travelTime = tt.time;

  computeApproximatedTakeOffAngles(eventLat, eventLon, eventDepth, station,
                                   phaseType, &azimuth, &takeOffAngle);
  //
  // Check for not fully implemented Travel Time Tables
  // i.e. LOCSAT doesn't compute tt.takeoff, tt.dtdh
  // and tt.dtdd is very different from other models
  //
  if (_type == "LOCSAT")
  {
    velocityAtSrc =
        computeVelocityAtSource(eventLat, eventLon, eventDepth, phaseType);
  }
  else
  {
    velocityAtSrc =
        1.0 / std::sqrt(square(tt.dtdh) +
                        square(Seiscomp::Math::Geo::km2deg(tt.dtdd)));
  }

  if (_type != "LOCSAT")
  {
    takeOffAngle = degToRad(tt.takeoff);
  }
}

// reverse-engineer the velocity at source
double TravelTimeTable::computeVelocityAtSource(double eventLat,
                                                double eventLon,
                                                double eventDepth,
                                                const std::string &phaseType)
{
  const int bin              = std::floor(eventDepth / _depthVelResolution);
  const double binStartDepth = bin * _depthVelResolution;
  const double binEndDepth   = (bin + 1) * _depthVelResolution;

  //
  // first check if we have already computed the velocity for this
  // phase/depth
  //
  auto it1 = _depthVel.find(phaseType);
  if (it1 != _depthVel.end())
  {
    const auto &phaseDepthVel = it1->second;
    auto it2                  = phaseDepthVel.find(bin);
    if (it2 != phaseDepthVel.end())
    {
      return it2->second;
    }
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
  if (tt1 < 0 || tt2 < 0)
  {
    throw Exception("Unable to compute velocity at source");
  }

  double binVelocity = _depthVelResolution / std::abs(tt2 - tt1); // [km/sec]

  if (!std::isfinite(binVelocity))
  {
    throw Exception("Unable to compute velocity at source");
  }

  // cache the value
  _depthVel[phaseType][bin] = binVelocity;

  return binVelocity;
}

} // namespace SCAdapter
} // namespace HDD

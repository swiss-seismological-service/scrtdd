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

TravelTimeTable::TravelTimeTable(std::string type,
                                 std::string model,
                                 double depthVelResolution)
    : _depthVelResolution(depthVelResolution)
{
  _ttt = TravelTimeTableInterface::Create(type.c_str());
  _ttt->setModel(model.c_str());

  /*
  for (int i = 0; i < 500; i++)
  {
    double vel = velocityAtSource(i * _depthVelResolution, "P");
    SEISCOMP_INFO("Velocity model phase P depth %.2f [km] vel %.2f [m/sec]",
                   (i * _depthVelResolution), vel);
  }
  for (int i = 0; i < 500; i++)
  {
    double vel = velocityAtSource(i * _depthVelResolution, "S");
    SEISCOMP_INFO("Velocity model phase S depth %.2f [km] vel %.2f [m/sec]",
                   (i * _depthVelResolution), vel);
  }
  */
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
  double depth  = eventDepth > 0 ? eventDepth : 0;
  TravelTime tt = _ttt->compute(phaseType.c_str(), eventLat, eventLon, depth,
                                stationLat, stationLon, stationElevation);
  travelTime    = tt.time;
  velocityAtSrc = velocityAtSource(depth, phaseType);
  // tt.takeOff is not computed for LOCSAT and for libTau it seems wrong
  takeOffAngle = 0; 
}

// Since the seiscomp travel time api doesn't offer the velocity at
// source information we need to reverse-engineer that information
double TravelTimeTable::velocityAtSource(double eventDepth,
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
  // This is a new phase/depth pair
  //
  double binVelocity = 0;
  try
  {
    double tt1 = (binStartDepth == 0) ? 0 :
        _ttt->compute(phaseType.c_str(), 0, 0, binStartDepth, 0, 0, 0).time;
    double tt2 = 
        _ttt->compute(phaseType.c_str(), 0, 0, binEndDepth, 0, 0, 0).time;

    binVelocity = _depthVelResolution / (tt2 - tt1); // [km/sec]
  }
  catch (exception &e)
  {
    SEISCOMP_WARNING("TTT: %s", e.what());
  }

  // store the value
  _depthVel[phaseType][bin] = binVelocity;

  //SEISCOMP_DEBUG("velocityAtSource %s bin %u depth %.3f[km] vel %.3f [m/sec]",
  //               phaseType.c_str(), bin, bin * _depthVelResolution,
  //               binVelocity);

  return binVelocity;
}

} // namespace HDD
} // namespace Seiscomp

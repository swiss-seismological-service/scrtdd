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
#include "hdd/log.h"
#include "hdd/utils.h"
#include <seiscomp/math/geo.h>

using namespace std;
using namespace HDD;

namespace HDD {
namespace SCAdapter {

TravelTimeTable::TravelTimeTable(const std::string &type,
                                 const std::string &model)
    : _type(type), _model(model)
{
  _ttt = Seiscomp::TravelTimeTableInterface::Create(_type.c_str());
  if (!_ttt || !_ttt->setModel(_model.c_str()))
  {
    throw Exception(strf("Unable to set travel time format %s and model %s",
                         _type.c_str(), _model.c_str()));
  }
}

double TravelTimeTable::compute(double eventLat,
                                double eventLon,
                                double eventDepth,
                                double stationLat,
                                double stationLon,
                                double stationElevation,
                                const std::string &phaseType)
{
#if SC_API_VERSION >= SC_API_VERSION_CHECK(16, 0, 0)
  double ttime =
      _ttt->computeTime(phaseType.c_str(), eventLat, eventLon, eventDepth,
                        stationLat, stationLon, stationElevation);
#else
  double ttime =
      _ttt->compute(phaseType.c_str(), eventLat, eventLon, eventDepth,
                    stationLat, stationLon, stationElevation)
          .time;
#endif

  if (ttime < 0)
  {
    throw Exception("No travel time data available");
  }

  return ttime;
}

void TravelTimeTable::compute(double eventLat,
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
  Seiscomp::TravelTime tt =
      _ttt->compute(phaseType.c_str(), eventLat, eventLon, eventDepth,
                    stationLat, stationLon, stationElevation);
  if (tt.time < 0)
  {
    throw Exception("No travel time data available");
  }

  travelTime = tt.time;

  // We want dtdd and dtdh to use the same units:
  // - transform dtdd [sec/rad] -> [sec/km]
  double dtdd2  = tt.dtdd / kmOfDegree(eventDepth);
  velocityAtSrc = 1.0 / std::sqrt(square(tt.dtdh) + square(dtdd2));

#if SC_API_VERSION >= SC_API_VERSION_CHECK(16, 0, 0)
  if (tt.azi) // 3D model
  {
    azimuth = *tt.azi;
  }
  else
#endif
  {
    azimuth = computeAzimuth(eventLat, eventLon, stationLat, stationLon);
    azimuth = radToDeg(azimuth);
  }

  takeOffAngle = tt.takeoff;
}

} // namespace SCAdapter
} // namespace HDD

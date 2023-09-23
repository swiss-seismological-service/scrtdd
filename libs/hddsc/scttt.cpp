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
                                 const std::string &model)
    : _type(type), _model(model)
{
  load();
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

  double ttime = _ttt->computeTravelTime(phaseType.c_str(), eventLat, eventLon,
                                         eventDepth, station.latitude,
                                         station.longitude, station.elevation);
  if (ttime < 0)
  {
    throw Exception("No travel time data available");
  }

  return ttime;
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

  // We want dtdd and dtdh to use the same units:
  // - transform dtdd [sec/rad] -> [sec/km]
  double dtdd2  = tt.dtdd / kmOfDegree(eventDepth);
  velocityAtSrc = 1.0 / std::sqrt(square(tt.dtdh) + square(dtdd2));

  if (tt.azi) // 3D model
  {
    azimuth = *tt.azi;
  }
  else
  {
    azimuth =
        computeAzimuth(eventLat, eventLon, station.latitude, station.longitude);
  }

  takeOffAngle = degToRad(tt.takeoff);
}

} // namespace SCAdapter
} // namespace HDD

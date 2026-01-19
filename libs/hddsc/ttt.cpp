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
  double ttime;
  try
  {
#if SC_API_VERSION >= SC_API_VERSION_CHECK(16, 0, 0)
    ttime = _ttt->computeTime(phaseType.c_str(), eventLat, eventLon, eventDepth,
                              stationLat, stationLon, stationElevation);
#else
    ttime = _ttt->compute(phaseType.c_str(), eventLat, eventLon, eventDepth,
                          stationLat, stationLon, stationElevation)
                .time;
#endif
  }
  catch (exception &e)
  {
    throw HDD::Exception(e.what());
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
                              double &takeOffAzi,
                              double &takeOffDip,
                              double &dtdd,
                              double &dtdh)
{
  Seiscomp::TravelTime tt;

  try
  {
    tt = _ttt->compute(phaseType.c_str(), eventLat, eventLon, eventDepth,
                       stationLat, stationLon, stationElevation);
  }
  catch (exception &e)
  {
    throw HDD::Exception(e.what());
  }

  travelTime = tt.time;
  dtdd       = tt.dtdd;
  dtdh       = tt.dtdh;
  takeOffDip = tt.takeoff;

#if SC_API_VERSION >= SC_API_VERSION_CHECK(16, 0, 0)
  if (tt.azi) // 3D model
  {
    takeOffAzi = *tt.azi;
  }
  else
#endif
  {
    takeOffAzi = computeAzimuth(eventLat, eventLon, stationLat, stationLon);
    takeOffAzi = radToDeg(takeOffAzi);
  }
}

} // namespace SCAdapter
} // namespace HDD

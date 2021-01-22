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

#ifndef __HDD_TTT_H__
#define __HDD_TTT_H__

#include "catalog.h"

#include <seiscomp3/core/baseobject.h>
#include <seiscomp3/seismology/ttt.h>

#include <unordered_map>

namespace Seiscomp {
namespace HDD {

DEFINE_SMARTPOINTER(TravelTimeTable);

class TravelTimeTable : public Core::BaseObject
{
public:
  TravelTimeTable(std::string type,
                  std::string model,
                  double depthVelResolution = 0.1);
  virtual ~TravelTimeTable() {}

  void compute(double eventLat,
               double eventLon,
               double eventDepth,
               double stationLat,
               double stationLon,
               double stationElevation,
               const std::string &phaseType,
               double &travelTime,
               double &takeOffAngle,
               double &velocityAtSrc);

  void compute(const Catalog::Event &event,
               const Catalog::Station &station,
               const std::string &phaseType,
               double &travelTime,
               double &takeOffAngle,
               double &velocityAtSrc)
  {
    return compute(event.latitude, event.longitude, event.depth,
                   station.latitude, station.longitude, station.elevation,
                   phaseType, travelTime, takeOffAngle, velocityAtSrc);
  }

private:
  double velocityAtSource(double eventDepth, const std::string &phaseType);

  TravelTimeTableInterfacePtr _ttt;
  const double _depthVelResolution; // km
  // key 1 = phase type. key 2 = depth bin
  std::unordered_map<std::string, std::unordered_map<int, double>> _depthVel;
};

DEFINE_SMARTPOINTER(TravelTimeTable);

} // namespace HDD
} // namespace Seiscomp

#endif

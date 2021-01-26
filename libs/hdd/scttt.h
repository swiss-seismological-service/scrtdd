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

#ifndef __HDD_SCTTT_H__
#define __HDD_SCTTT_H__

#include "ttt.h"

#include <seiscomp3/core/baseobject.h>
#include <seiscomp3/seismology/ttt.h>

#include <unordered_map>

namespace Seiscomp {
namespace HDD {

class ScTravelTimeTable : public TravelTimeTable
{
public:
  ScTravelTimeTable(const std::string &type,
                    const std::string &model,
                    double depthVelResolution = 0.1);
  virtual ~ScTravelTimeTable() {}

  virtual void compute(double eventLat,
                       double eventLon,
                       double eventDepth,
                       const Catalog::Station &station,
                       const std::string &phaseType,
                       double &travelTime,
                       double &takeOfAngleAzim,
                       double &takeOfAngleDip,
                       double &velocityAtSrc);

  virtual void compute(double eventLat,
                       double eventLon,
                       double eventDepth,
                       const Catalog::Station &station,
                       const std::string &phaseType,
                       double &travelTime);

private:
  double velocityAtSource(double eventDepth, const std::string &phaseType);

  TravelTimeTableInterfacePtr _ttt;
  const double _depthVelResolution; // km
  // key 1 = phase type. key 2 = depth bin
  std::unordered_map<std::string, std::unordered_map<int, double>> _depthVel;
};

} // namespace HDD
} // namespace Seiscomp

#endif

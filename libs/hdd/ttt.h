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

#ifndef __HDD_TTT_H__
#define __HDD_TTT_H__

#include "catalog.h"

#include <seiscomp3/core/baseobject.h>
#include <seiscomp3/seismology/ttt.h>

namespace Seiscomp {
namespace HDD {

DEFINE_SMARTPOINTER(TravelTimeTable);

class TravelTimeTable : public Core::BaseObject
{
public:
  static TravelTimeTable *create(const std::string &type,
                                 const std::string &model);

  /*
   * The implementation of this interface MUST compute:
   * - travelTime [seconds]
   * - takeOffAngleAzim [degree]
   * - takeOffAngleDip [degree]
   * - velocityAtSrc [km/sec]
   */
  virtual void compute(double eventLat,
                       double eventLon,
                       double eventDepth,
                       const Catalog::Station &station,
                       const std::string &phaseType,
                       double &travelTime,
                       double &takeOfAngleAzim,
                       double &takeOfAngleDip,
                       double &velocityAtSrc) = 0;

  void compute(const Catalog::Event &event,
               const Catalog::Station &station,
               const std::string &phaseType,
               double &travelTime,
               double &takeOfAngleAzim,
               double &takeOfAngleDip,
               double &velocityAtSrc)
  {
    return compute(event.latitude, event.longitude, event.depth, station,
                   phaseType, travelTime, takeOfAngleAzim, takeOfAngleDip,
                   velocityAtSrc);
  }

  /*
   * Sometimes only travel time is required (save computation)
   */
  virtual void compute(double eventLat,
                       double eventLon,
                       double eventDepth,
                       const Catalog::Station &station,
                       const std::string &phaseType,
                       double &travelTime) = 0;

  void compute(const Catalog::Event &event,
               const Catalog::Station &station,
               const std::string &phaseType,
               double &travelTime)
  {
    return compute(event.latitude, event.longitude, event.depth, station,
                   phaseType, travelTime);
  }

  const std::string type;
  const std::string model;

protected:
  /*
   * Utility function to compute straight ray path approximation
   * for takeOfAngles angles when those are not available in the
   * travel time tables
   */
  static void computeApproximatedTakeOfAngles(double eventLat,
                                              double eventLon,
                                              double eventDepth,
                                              const Catalog::Station &station,
                                              const std::string &phaseType,
                                              double *takeOfAngleAzim = nullptr,
                                              double *takeOfAngleDip = nullptr);

  TravelTimeTable(const std::string &t, const std::string &m)
      : type(t), model(m)
  {}
  virtual ~TravelTimeTable() {}
};

} // namespace HDD
} // namespace Seiscomp

#endif

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

#ifndef __HDD_CVTTT_H__
#define __HDD_CVTTT_H__

#include "ttt.h"

namespace HDD {

class ConstantVelocity : public HDD::TravelTimeTable
{
public:
  ConstantVelocity(double pVel, double sVel);

  virtual ~ConstantVelocity() = default;

  void freeResources() override;

  void compute(double eventLat,
               double eventLon,
               double eventDepth,
               const Catalog::Station &station,
               const std::string &phaseType,
               double &travelTime) override;

  void compute(double eventLat,
               double eventLon,
               double eventDepth,
               const Catalog::Station &station,
               const std::string &phaseType,
               double &travelTime,
               double &takeOffAngleAzim,
               double &takeOffAngleDip,
               double &velocityAtSrc) override;

private:
  const double _pVel, _sVel; // [km/sec]
};

} // namespace HDD

#endif

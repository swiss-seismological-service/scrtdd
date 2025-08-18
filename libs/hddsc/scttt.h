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

#ifndef __HDD_SCTTT_H__
#define __HDD_SCTTT_H__

#include "hdd/ttt.h"

#include <seiscomp/seismology/ttt.h>

#include <unordered_map>

namespace HDD {
namespace SCAdapter {

class TravelTimeTable : public HDD::TravelTimeTable
{
public:
  TravelTimeTable(const std::string &type,
                  const std::string &model);
  virtual ~TravelTimeTable() = default;

  void freeResources() override;

  void compute(double eventLat,
               double eventLon,
               double eventDepth,
               const Catalog::Station &station,
               const std::string &phaseType,
               double &travelTime,
               double &azimuth,
               double &takeOffAngle,
               double &velocityAtSrc) override;

  double compute(double eventLat,
                 double eventLon,
                 double eventDepth,
                 const Catalog::Station &station,
                 const std::string &phaseType) override;

private:
  void load();

  const std::string _type;
  const std::string _model;
  Seiscomp::TravelTimeTableInterfacePtr _ttt;
};

} // namespace SCAdapter
} // namespace HDD
#endif

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

#ifndef __HDD_NLLTTT_H__
#define __HDD_NLLTTT_H__

#include "3rd-party/lrucache.h"
#include "nllgrid.h"
#include "ttt.h"

#include <unordered_set>

namespace HDD {
namespace NLL {

class TravelTimeTable : public HDD::TravelTimeTable
{
public:
  TravelTimeTable(const std::string &velGridPath,
                  const std::string &timeGridPath,
                  const std::string &angleGridPath,
                  bool swapBytes,
                  unsigned cacheSize = 255);

  virtual ~TravelTimeTable() = default;

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
  std::string _velGridPath;
  std::string _timeGridPath;
  std::string _angleGridPath;
  bool _swapBytes;
  lru_cache<std::string, std::shared_ptr<VelGrid>> _velGrids;
  lru_cache<std::string, std::shared_ptr<TimeGrid>> _timeGrids;
  lru_cache<std::string, std::shared_ptr<AngleGrid>> _angleGrids;
  std::unordered_set<std::string> _unloadableGrids;
};

} // namespace NLL
} // namespace HDD

#endif

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

#ifndef __HDD_CLUSTERING_H__
#define __HDD_CLUSTERING_H__

#include "catalog.h"
#include <list>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace HDD {

// DD background catalog
struct Neighbours
{
  unsigned refEvId;
  std::unordered_set<unsigned> ids; // neighbouring event id
  std::unordered_map<
      unsigned,                       // indexed by event id
      std::unordered_map<std::string, // indexed by station id
                         std::unordered_set<Catalog::Phase::Type>>>
      phases;

  void add(unsigned neighbourId,
           const std::string &stationId,
           const Catalog::Phase::Type &phase)
  {
    ids.insert(neighbourId);
    phases[neighbourId][stationId].insert(phase);
  }

  bool has(unsigned neighbourId) const
  {
    return ids.find(neighbourId) != ids.end();
  }

  bool has(unsigned neighbourId, const std::string stationId) const
  {
    const auto &neighPhases = phases.find(neighbourId);
    if (neighPhases != phases.end())
      return neighPhases->second.find(stationId) != neighPhases->second.end();
    return false;
  }

  bool has(unsigned neighbourId,
           const std::string stationId,
           Catalog::Phase::Type type) const
  {
    const auto &neighPhases = phases.find(neighbourId);
    if (neighPhases != phases.end())
    {
      const auto &neighPhaseTypes = neighPhases->second.find(stationId);
      if (neighPhaseTypes != neighPhases->second.end())
        return neighPhaseTypes->second.find(type) !=
               neighPhaseTypes->second.end();
    }
    return false;
  }

  std::unordered_map<std::string, std::unordered_set<Catalog::Phase::Type>>
  allPhases() const;

  std::unordered_map<std::string, std::unordered_set<Catalog::Phase::Type>>
  allPhases(unsigned neighbourId) const;

  std::unique_ptr<Catalog> toCatalog(const Catalog &catalog,
                                     bool includeRefEv = false) const;
};

std::unique_ptr<Neighbours>
selectNeighbouringEvents(const Catalog &catalog,
                         const Catalog::Event &refEv,
                         const Catalog &refEvCatalog,
                         double minPhaseWeight   = 0,
                         double minESdis         = 0,
                         double maxESdis         = -1, // -1 = no limits
                         double minEStoIEratio   = 0,
                         unsigned minDTperEvt    = 1,
                         unsigned maxDTperEvt    = 0, // 0 = no limits
                         unsigned minNumNeigh    = 1,
                         unsigned maxNumNeigh    = 0, // 0 = no limits
                         unsigned numEllipsoids  = 5,
                         double maxEllipsoidSize = 10,
                         bool keepUnmatched      = false);

// find Neighbours for each event in the catalog
std::unordered_map<unsigned, std::unique_ptr<Neighbours>>
selectNeighbouringEventsCatalog(const Catalog &catalog,
                                double minPhaseWeight,
                                double minESdis,
                                double maxESdis,
                                double minEStoIEratio,
                                unsigned minDTperEvt,
                                unsigned maxDTperEvt,
                                unsigned minNumNeigh,
                                unsigned maxNumNeigh,
                                unsigned numEllipsoids,
                                double maxEllipsoidSize,
                                bool keepUnmatched);

// Organize the neighbours by not connected clusters. In addition,
// don't report the same pair multiple times (e.g. ev1-ev2 and ev2-ev1)
// since we only need one observation per pair in the DD solver.
// The input will be moved to the return value
std::list<std::unordered_map<unsigned, std::unique_ptr<Neighbours>>>
clusterizeNeighbouringEvents(
    std::unordered_map<unsigned, std::unique_ptr<Neighbours>>
        &neighboursByEvent);

} // namespace HDD

#endif

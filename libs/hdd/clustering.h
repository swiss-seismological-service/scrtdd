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
#include <deque>
#include <list>
#include <seiscomp3/core/baseobject.h>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace Seiscomp {
namespace HDD {

DEFINE_SMARTPOINTER(Neighbours);

// DD background catalog
struct Neighbours : public Core::BaseObject
{
  unsigned refEvId;

  std::unordered_set<unsigned> ids; // neighbouring event id

  std::unordered_map<unsigned,                       // indexed by event id
                     std::unordered_map<std::string, // indexed by station id
                                        std::set<Catalog::Phase::Type>>>
      phases;

  void add(unsigned neighbourId,
           const std::string &stationId,
           const Catalog::Phase::Type &phase)
  {
    ids.insert(neighbourId);
    phases[neighbourId][stationId].insert(phase);
  }

  std::unordered_set<unsigned>::size_type numNeighbours() const
  {
    return ids.size();
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

  std::unordered_map<std::string, std::set<Catalog::Phase::Type>>
  allPhases() const;

  CatalogPtr toCatalog(const CatalogCPtr &catalog,
                       bool includeRefEv = false) const;
};

NeighboursPtr
selectNeighbouringEvents(const CatalogCPtr &catalog,
                         const Catalog::Event &refEv,
                         const CatalogCPtr &refEvCatalog,
                         double minPhaseWeight   = 0,
                         double minESdis         = 0,
                         double maxESdis         = -1,
                         double minEStoIEratio   = 0,
                         unsigned minDTperEvt    = 1,
                         unsigned maxDTperEvt    = 0, // 0 = no limits
                         unsigned minNumNeigh    = 1,
                         unsigned maxNumNeigh    = 0, // 0 = no limits
                         unsigned numEllipsoids  = 5,
                         double maxEllipsoidSize = 10,
                         bool keepUnmatched      = false);

std::list<NeighboursPtr>
selectNeighbouringEventsCatalog(const CatalogCPtr &catalog,
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

std::deque<std::list<NeighboursPtr>>
clusterizeNeighbouringEvents(const std::list<NeighboursPtr> &neighboursList);

} // namespace HDD
} // namespace Seiscomp

#endif

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

#ifndef __RTDD_APPLICATIONS_CLUSTERING_H__
#define __RTDD_APPLICATIONS_CLUSTERING_H__

#include "catalog.h"
#include <seiscomp3/core/baseobject.h>
#include <unordered_map>
#include <unordered_set>
#include <set>


namespace Seiscomp {
namespace HDD { 

DEFINE_SMARTPOINTER(Neighbours);

// DD background catalog
struct Neighbours : public Core::BaseObject
{
    unsigned refEvId;

    unsigned numNeighbours;

    std::unordered_set<unsigned> ids; // neighbouring event id 

    std::unordered_map<unsigned, // indexed by event id
         std::unordered_map<std::string, // indexed by station id
                            std::set<Catalog::Phase::Type> > > phases;

    std::unordered_map<std::string, std::set<Catalog::Phase::Type> > allPhases() const
    {
         std::unordered_map<std::string, std::set<Catalog::Phase::Type> > allPhases;
         for ( const auto& kw1 : phases )
            for ( const auto& kw2 : kw1.second )
                 allPhases[kw2.first].insert(kw2.second.begin(), kw2.second.end());
         return allPhases;
    }

    bool has(unsigned neighbourId) const
    {
        return ids.find(neighbourId) != ids.end();
    }

    bool has(unsigned neighbourId, const std::string stationId) const
    {
        const auto& neighPhases = phases.find(neighbourId);
        if ( neighPhases != phases.end() )
            return neighPhases->second.find(stationId) != neighPhases->second.end();
        return false;
    }

    bool has(unsigned neighbourId, const std::string stationId, Catalog::Phase::Type type) const
    {
        const auto& neighPhases = phases.find(neighbourId);
        if ( neighPhases != phases.end() ) {
            const auto& neighPhaseTypes = neighPhases->second.find(stationId);
            if ( neighPhaseTypes != neighPhases->second.end() )
                return  neighPhaseTypes->second.find(type) != neighPhaseTypes->second.end();
        }
        return false;
    } 

    CatalogPtr fromNeighbours(const CatalogCPtr& catalog, bool includeRefEv=false) const
    {
        CatalogPtr returnCat( new Catalog() );
        for (unsigned neighbourId : ids)
            returnCat->add(neighbourId, *catalog, true);
        if ( includeRefEv )
            returnCat->add(refEvId, *catalog, true);
        return returnCat;
    }
};


NeighboursPtr
selectNeighbouringEvents(const CatalogCPtr& catalog,
                         const Catalog::Event& refEv,
                         const CatalogCPtr& refEvCatalog,
                         double minPhaseWeight = 0,
                         double minESdis=0,
                         double maxESdis=-1,
                         double minEStoIEratio=0,
                         int minDTperEvt=1,
                         int maxDTperEvt=-1,
                         int minNumNeigh=1,
                         int maxNumNeigh=-1,
                         int numEllipsoids=5,
                         double maxEllipsoidSize=10,
                         bool keepUnmatched=false);

std::list<NeighboursPtr>
selectNeighbouringEventsCatalog(const CatalogCPtr& catalog,
                                double minPhaseWeight,
                                double minESdis,
                                double maxESdis,
                                double minEStoIEratio,
                                int minDTperEvt,
                                int maxDTperEvt,
                                int minNumNeigh,
                                int maxNumNeigh,
                                int numEllipsoids,
                                double maxEllipsoidSize,
                                bool keepUnmatched);

}
}

#endif 


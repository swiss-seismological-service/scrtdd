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

#include "clustering.h"
#include "utils.h"
#include "ellipsoid.ipp"

#include <seiscomp3/core/strings.h>

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp3/logging/log.h>

using namespace std;
using namespace Seiscomp;
using Seiscomp::Core::stringify;
using Event = HDD::Catalog::Event;
using Phase = HDD::Catalog::Phase;
using Station = HDD::Catalog::Station;

namespace Seiscomp {
namespace HDD {


NeighboursPtr
selectNeighbouringEvents(const CatalogCPtr& catalog,
                         const Event& refEv,
                         const CatalogCPtr& refEvCatalog,
                         double minPhaseWeight,
                         double minESdist,
                         double maxESdist,
                         double minEStoIEratio,
                         unsigned minDTperEvt,
                         unsigned maxDTperEvt,
                         unsigned minNumNeigh, 
                         unsigned maxNumNeigh,
                         unsigned numEllipsoids,
                         double maxEllipsoidSize,
                         bool keepUnmatched)
{
    SEISCOMP_INFO("Selecting Neighbouring Events for event %s lat %.6f lon %.6f depth %.4f mag %.2f time %s",
                   string(refEv).c_str(), refEv.latitude, refEv.longitude, refEv.depth,
                   refEv.magnitude, refEv.time.iso().c_str());

    // Optimization: make code faster but the result will be the same
    if ( maxNumNeigh <= 0 )
    {
        SEISCOMP_INFO("Disabling ellipsoid algorithm since maxNumNeigh is not set");
        numEllipsoids = 0;
    }

    /*
     * Build ellipsoids
     *
     * From Waldhauser 2009: to assure a spatially homogeneous subsampling, reference
     * events are selected within each of five concentric, vertically elongated
     * ellipsoidal layers of increasing thickness. Each layer has 8 quadrants.
     */
    vector<HddEllipsoidPtr> ellipsoids;
    double verticalSize = maxEllipsoidSize * 2; // horizontal to vertical axis length 
    for ( unsigned i = 0; i < (numEllipsoids-1); i++ )
    {
        ellipsoids.push_back( new HddEllipsoid(verticalSize, refEv.latitude, refEv.longitude, refEv.depth) );
        verticalSize /= 2;
    }
    ellipsoids.push_back( new HddEllipsoid(0, verticalSize, refEv.latitude, refEv.longitude, refEv.depth) );

    //
    // sort catalog events by distance and drop the ones further than the outmost ellipsoid
    //
    unordered_map<unsigned,double> distanceByEvent; // eventid, distance
    unordered_map<unsigned,double> azimuthByEvent;  // eventid, azimuth

    for (const auto& kv : catalog->getEvents() )
    {
        const Event& event = kv.second;

        if (event == refEv)
            continue;

        // drop event if outside the outmost ellipsod boundaries
        HddEllipsoidPtr outmostEllip = ellipsoids[0];
        if ( ! outmostEllip->getOuterEllipsoid().isInside(event.latitude, event.longitude, event.depth) )
            continue;

        // compute distance between current event and reference origin
        double azimuth;
        double distance = computeDistance(refEv, event, &azimuth);

        // keep a list of added events
        distanceByEvent[event.id] = distance;
        azimuthByEvent[event.id]  = azimuth;
    }

    //
    // Select stations within configured distance
    //
    unordered_map<string,double> validatedStationDistance;
    for (const auto& kv : catalog->getStations()  )
    {
        const string& staId = kv.first;
        const Station& station = kv.second;

        // compute distance between reference event and station
        double staRefEvDistance = computeDistance(refEv, station);

        // check this station distance is ok
        if ( ( maxESdist <= 0 || staRefEvDistance <= maxESdist ) ||  // too far away ?
             ( staRefEvDistance >= minESdist )                  )     // too close ?
        {
            validatedStationDistance[staId] = staRefEvDistance;
        }
    } 

    //
    // Select from the events within distance the ones who respect the constraints
    //
    struct SelectedEventEntry {
        Event event;
        unordered_map<string, set<Phase::Type> > phases;
    };
    multimap<double,SelectedEventEntry> selectedEvents; // distance, struct
    unordered_map<unsigned,int> dtCountByEvent; // eventid, dtCount
    unordered_set<string> includedStations;
    unordered_set<string> excludedStations;

    for (const auto& kv : distanceByEvent)
    {
        const Event& event = catalog->getEvents().at(kv.first);
        const double eventDistance = kv.second;

        //
        // Loop through this even phases and keep track of
        // the valid phases and their station distance
        //
        multimap<double, pair<string,Phase::Type> > stationByDistance; // distance, <stationid,phaseType>
        multimap<double, pair<string,Phase::Type> > unmatchedPhases; // distance, <stationid,phaseType>

        auto eqlrng = catalog->getPhases().equal_range(event.id);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
            const Phase& phase = it->second;
            const Station& station = catalog->getStations().at(phase.stationId);

            // check pick weight
            if (phase.procInfo.weight < minPhaseWeight)
                continue;

            // check this station distance to reference event is ok
            const auto& staRefEvDistanceIt = validatedStationDistance.find(phase.stationId);
            if ( staRefEvDistanceIt == validatedStationDistance.end() )
                continue;

            double staRefEvDistance = staRefEvDistanceIt->second;

            if ( (staRefEvDistance / eventDistance) < minEStoIEratio )  // ratio too small ?
            {
                continue;
            } 

            // check this station distance to current event is ok
            if ( maxESdist > 0 )
            {
                // compute distance between current event and station
                double stationDistance = computeDistance(event, station);

                if ( ( stationDistance > maxESdist )                 ||      // too far away ?
                     ( stationDistance < minESdist )                 ||       // too close ?
                     ( (stationDistance / eventDistance) < minEStoIEratio ) ) // ratio too small ?
                {
                    continue;
                }
            }

            // now find corresponding phase in reference event phases
            bool peer_found = false;
            auto itRef = refEvCatalog->searchPhase(refEv.id, phase.stationId, phase.procInfo.type);
            if ( itRef != refEvCatalog->getPhases().end() )
            {
                const Phase& refPhase = itRef->second;
                if (refPhase.procInfo.weight >= minPhaseWeight)
                    peer_found = true;
            }

            if ( ! peer_found )
            {
                unmatchedPhases.emplace(staRefEvDistance, pair<string,Phase::Type>(phase.stationId, phase.procInfo.type));
                continue;
            }

            stationByDistance.emplace(staRefEvDistance, pair<string,Phase::Type>(phase.stationId, phase.procInfo.type));
        }

        // Check enough phases (> minDTperEvt) ? if not skip event
        if ( stationByDistance.size() < minDTperEvt )
        {
            continue;
        }

        unsigned numObservations = 0;

        // Since the constraints are met evSelEntry will be added to selectedEvents
        SelectedEventEntry evSelEntry;
        evSelEntry.event = event;

        // Copy the phases
        for ( const auto& kw : stationByDistance )
        {
            // if maxDTperEvt is set then make sure to stay within limits
            if ( maxDTperEvt > 0 && numObservations >= maxDTperEvt )
                break;
            const pair<string,Phase::Type>& staPh = kw.second;
            evSelEntry.phases[ staPh.first ].insert( staPh.second );
            numObservations++;
        }

        if ( keepUnmatched )
        {
            // Those neighbour's phases respects the constraints and if we add
            // theoretical picks to refEv those phases might become useful for xcorr
            for ( const auto& kw : unmatchedPhases )
            {
                // if maxDTperEvt is set then make sure to stay within limits
                if ( maxDTperEvt > 0 && numObservations >= maxDTperEvt )
                    break;
                const pair<string,Phase::Type>& staPh = kw.second;
                evSelEntry.phases[ staPh.first ].insert( staPh.second );
                numObservations++;
            }
        }

        // add this event to the selected ones
        selectedEvents.emplace(eventDistance, evSelEntry);
        dtCountByEvent.emplace(event.id, numObservations);
    }

    //
    // Finally build the catalog of neighboring events using the elipsoids
    // algorithm or simply the nearest neighbour method
    //
    NeighboursPtr neighboringEventCat( new Neighbours() );
    neighboringEventCat->refEvId = refEv.id;
    neighboringEventCat->numNeighbours = 0;

    if ( numEllipsoids <= 0 )
    {
        //
        // if numEllipsoids is 0 then disable the ellipsoid algorithm and simply select events
        // on the nearest neighbor basis
        // Since selectedEvents is sorted by distance we get closer events first
        //
        for (auto kv : selectedEvents)
        {
            const SelectedEventEntry& evSelEntry = kv.second;
            const Event& ev = evSelEntry.event;

            // add this event to the catalog
            neighboringEventCat->ids.insert( ev.id );
            neighboringEventCat->phases[ ev.id ] = evSelEntry.phases;
            neighboringEventCat->numNeighbours++;

            SEISCOMP_INFO("Neighbour: #obsers %2d distance %5.2f azimuth %3.f "
                          "depth-diff %6.3f depth %5.3f event %s",
                          dtCountByEvent[ev.id], distanceByEvent[ev.id], azimuthByEvent[ev.id], 
                          refEv.depth-ev.depth, ev.depth, string(ev).c_str() );

            if ( maxNumNeigh > 0 && neighboringEventCat->numNeighbours >= maxNumNeigh) break;
        }
    }
    else
    {
        //
        // Apply Waldauser's concentric ellipsoids algorithm
        //
        vector<int> quadrants = {1,2,3,4,5,6,7,8};
        bool workToDo = true;

        while ( workToDo )
        {
            for(int elpsNum = ellipsoids.size() -1; elpsNum >= 0;  elpsNum--)
            {
                for (int quadrant : quadrants)
                {
                    // if we don't have events or we have already selected maxNumNeigh neighbors exit
                    if ( selectedEvents.empty() ||
                        (maxNumNeigh > 0 && neighboringEventCat->numNeighbours >= maxNumNeigh) )
                    {
                        workToDo = false;
                        break;
                    }

                    // since selectedEvents is sorted by distance we get closer events first
                    for (auto it = selectedEvents.begin(); it != selectedEvents.end(); it++)
                    {
                        const SelectedEventEntry& evSelEntry = it->second;
                        const Event& ev = evSelEntry.event;

                        bool found = ellipsoids[elpsNum]->isInside(ev.latitude, ev.longitude,
                                                                   ev.depth, quadrant);
                        if ( found )
                        {
                            // add this event to the catalog
                            neighboringEventCat->ids.insert( ev.id );
                            neighboringEventCat->phases[ ev.id ] = evSelEntry.phases;

                            neighboringEventCat->numNeighbours++;
                            selectedEvents.erase(it);

                            SEISCOMP_INFO("Neighbour: ellipsoid %2d quadrant %d #observs %2d "
                              "distance %5.2f azimuth %3.f depth-diff %6.3f depth %5.3f event %s ",
                              elpsNum, quadrant, dtCountByEvent[ev.id], distanceByEvent[ev.id],
                              azimuthByEvent[ev.id], refEv.depth-ev.depth, ev.depth,
                              string(ev).c_str() );

                            break;
                        }
                    }
                }
            }
        }
    }

    // Check if enough neighbors were found
    if ( neighboringEventCat->numNeighbours < minNumNeigh )
    {
        string msg = stringify("Skipping event %s, insufficient number of neighbors (%d)",
                               string(refEv).c_str(), neighboringEventCat->numNeighbours);
        SEISCOMP_DEBUG("%s", msg.c_str());
        throw runtime_error(msg);
    }

    return neighboringEventCat;
}



std::list<NeighboursPtr>
selectNeighbouringEventsCatalog(const CatalogCPtr& catalog,
                                double minPhaseWeight,
                                double minESdist,
                                double maxESdist,
                                double minEStoIEratio,
                                unsigned minDTperEvt,
                                unsigned maxDTperEvt,
                                unsigned minNumNeigh,
                                unsigned maxNumNeigh,
                                unsigned numEllipsoids,
                                double maxEllipsoidSize,
                                bool keepUnmatched)
{
    SEISCOMP_INFO("Selecting Catalog Neighbouring Events ");

    // neighbours for each event
    list<NeighboursPtr> neighboursByEvent;

    // for each event find the neighbours
    CatalogPtr validCatalog = new Catalog(*catalog);
    list<unsigned> todoEvents;
    for (const auto& kv : validCatalog->getEvents() )
        todoEvents.push_back( kv.first );

    while ( ! todoEvents.empty() )
    {
        list<NeighboursPtr> newNeighbourCats;
        vector<unsigned> removedEvents;

        // for each event find the neighbours
        for (unsigned evIdTodo : todoEvents )
        {
            Catalog::Event event = validCatalog->getEvents().find(evIdTodo)->second;

            NeighboursPtr neighbours;
            try {
                neighbours = selectNeighbouringEvents(
                    validCatalog, event, validCatalog,  minPhaseWeight, minESdist, maxESdist,
                    minEStoIEratio, minDTperEvt, maxDTperEvt,  minNumNeigh, maxNumNeigh,
                    numEllipsoids, maxEllipsoidSize, keepUnmatched
                );
            } catch ( ... ) { }

            if ( ! neighbours )
            {
                // event discarded because it doesn't satisfies requirements
                removedEvents.push_back( event.id );
                // next loop we don't want other events to pick this as neighbour
                validCatalog->removeEvent( event.id );
                todoEvents.remove( event.id ); // this invalidates the loop !
                // stop here because we dont' want to keep building potentially wrong neighbours
                break;
            }

            newNeighbourCats.push_back(neighbours);
        }

        // add newly computed neighbors catalogs to previous ones
        for ( NeighboursPtr& neighbours : newNeighbourCats )
        {
            neighboursByEvent.push_back( neighbours );
            // make sure we won't recompute what has been already done
            todoEvents.remove( neighbours->refEvId );
        }

        // check if the removed events were used as neighbour of any other event
        // if so rebuild neighbours for those events
        bool redo;
        do {
            redo = false;
            list<NeighboursPtr> validNeighbourCats;

            for ( NeighboursPtr& neighbours : neighboursByEvent )
            {
                bool currCatInvalid = false;
                for (unsigned removedEventId : removedEvents)
                {
                    if( neighbours->ids.count( removedEventId ) != 0 )
                    {
                        currCatInvalid = true;
                        break;
                    }
                }

                if ( currCatInvalid )
                {
                    removedEvents.push_back( neighbours->refEvId );
                    todoEvents.push_back( neighbours->refEvId );
                    redo = true;
                    continue;
                }
                validNeighbourCats.push_back( neighbours );
            }

            neighboursByEvent.clear();
            neighboursByEvent = validNeighbourCats;

        } while( redo );
    }

    // We don't want to report the same pairs multiple times
    // when creating double-difference observations. So we'll
    // remove the pairs that appeared in previous catalogs from
    // the following catalogs
    std::unordered_multimap<unsigned,unsigned> existingPairs;

    for ( NeighboursPtr& neighbours : neighboursByEvent )
    {
        unsigned currEventId = neighbours->refEvId;

        // remove from currrent catalog the existing pairs
        auto eqlrng = existingPairs.equal_range(currEventId );
        for (auto existingPair = eqlrng.first; existingPair != eqlrng.second; existingPair++)
        {
            neighbours->ids.erase( existingPair->second );
            neighbours->phases.erase( existingPair->second );
        }

        // remove current pairs from following catalogs
        for ( unsigned neighEvId : neighbours->ids )
        {
            existingPairs.emplace(neighEvId, currEventId);
        }
    }

    return neighboursByEvent;
}


}
}

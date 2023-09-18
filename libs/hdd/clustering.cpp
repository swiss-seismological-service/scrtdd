/***************************************************************************
 * MIT License                                                             *
 *                                                                         *
 * Copyright (C) by ETHZ/SED                                               *
 *                                                                         *
 * Permission is hereby granted, free of charge, to any person obtaining a *
 * copy of this software and associated documentation files (the           *
 * “Software”), to deal in the Software without restriction, including     *
 * without limitation the rights to use, copy, modify, merge, publish,     *
 * distribute, sublicense, and/or sell copies of the Software, and to      *
 * permit persons to whom the Software is furnished to do so, subject to   *
 * the following conditions:                                               *
 *                                                                         *
 * The above copyright notice and this permission notice shall be          *
 * included in all copies or substantial portions of the Software.         *
 *                                                                         *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,         *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      *
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  *
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY    *
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,    *
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE       *
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                  *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#include "clustering.h"
#include "ellipsoid.h"
#include "log.h"
#include "utils.h"

using namespace std;
using Event   = HDD::Catalog::Event;
using Phase   = HDD::Catalog::Phase;
using Station = HDD::Catalog::Station;

namespace HDD {

unordered_map<string, unordered_set<Catalog::Phase::Type>>
Neighbours::allPhases() const
{
  unordered_map<string, unordered_set<Catalog::Phase::Type>> allPhases;
  for (const auto &kv1 : phases)
    for (const auto &kv2 : kv1.second)
      allPhases[kv2.first].insert(kv2.second.begin(), kv2.second.end());
  return allPhases;
}

unordered_map<string, unordered_set<Catalog::Phase::Type>>
Neighbours::allPhases(unsigned neighbourId) const
{
  unordered_map<string, unordered_set<Catalog::Phase::Type>> allPhases;
  try
  {
    for (const auto &kv : phases.at(neighbourId))
      allPhases[kv.first].insert(kv.second.begin(), kv.second.end());
  }
  catch (...)
  {}
  return allPhases;
}

unique_ptr<Catalog> Neighbours::toCatalog(const Catalog &catalog,
                                          bool includeRefEv) const
{
  unique_ptr<Catalog> returnCat(new Catalog());
  for (unsigned neighbourId : ids) returnCat->add(neighbourId, catalog, true);
  if (includeRefEv) returnCat->add(refEvId, catalog, true);
  return returnCat;
}

unique_ptr<Neighbours> selectNeighbouringEvents(const Catalog &catalog,
                                                const Event &refEv,
                                                const Catalog &refEvCatalog,
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
  logDebug("Selecting Neighbouring Events for event %s lat %g lon %g depth %g",
           string(refEv).c_str(), refEv.latitude, refEv.longitude, refEv.depth);

  // Optimization: make code faster but the result will be the same.
  if (maxNumNeigh <= 0)
  {
    logDebug("Disabling ellipsoid algorithm since maxNumNeigh is not set");
    numEllipsoids = 0;
  }

  /*
   * Build ellipsoids
   *
   * From Waldhauser 2009: to assure a spatially homogeneous subsampling,
   * reference events are selected within each of five concentric, vertically
   * elongated ellipsoidal layers of increasing thickness. Each layer has 8
   * quadrants.
   */
  vector<HddEllipsoid> ellipsoids;
  double verticalSize =
      maxEllipsoidSize * 2; // horizontal to vertical axis length
  for (unsigned i = 1; i < numEllipsoids; i++)
  {
    ellipsoids.push_back(HddEllipsoid(verticalSize, refEv.latitude,
                                      refEv.longitude, refEv.depth));
    verticalSize /= 2;
  }
  ellipsoids.push_back(HddEllipsoid(0, verticalSize, refEv.latitude,
                                    refEv.longitude, refEv.depth));

  //
  // Sort catalog events by distance and drop the ones further than the outmost
  // ellipsoid.
  //
  multimap<double, unsigned> eventByDistance;      // distance, eventid
  unordered_map<unsigned, double> distanceByEvent; // eventid, distance

  for (const auto &kv : catalog.getEvents())
  {
    const Event &event = kv.second;

    if (event == refEv) continue;

    // drop event if outside the outmost ellipsod boundaries
    const HddEllipsoid &outmostEllip = ellipsoids.at(0);
    if (!outmostEllip.getOuterEllipsoid().isInside(
            event.latitude, event.longitude, event.depth))
      continue;

    // compute distance between current event and reference origin
    double distance = computeDistance(refEv, event);

    // keep a list of events in range sorted by distance
    eventByDistance.emplace(distance, event.id);
    distanceByEvent.emplace(event.id, distance);
  }

  //
  // select stations within configured distance
  //
  unordered_map<string, double> validatedStationDistance;
  for (const auto &kv : catalog.getStations())
  {
    const string &staId    = kv.first;
    const Station &station = kv.second;

    // compute distance between reference event and station
    double staRefEvDistance = computeDistance(refEv, station);

    // check if station distance is ok
    if ((maxESdist <= 0 || staRefEvDistance <= maxESdist) || // too far away ?
        (staRefEvDistance >= minESdist))                     // too close ?
    {
      validatedStationDistance.emplace(staId, staRefEvDistance);
    }
  }

  //
  // Select from the events within distance the ones which respect the
  // constraints.
  //
  struct SelectedEventEntry
  {
    Event event;
    unordered_map<string, unordered_set<Phase::Type>> phases;
  };
  multimap<double, SelectedEventEntry> selectedEvents; // distance, struct
  unordered_map<unsigned, int> dtCountByEvent;         // eventid, dtCount

  for (const auto &kv : eventByDistance)
  {
    const double eventDistance = kv.first;
    const Event &event         = catalog.getEvents().at(kv.second);

    // Loop through event phases and keep track of the valid phases and their
    // station distance.
    multimap<double, pair<string, Phase::Type>>
        stationByDistance; // distance, <stationid,phaseType>
    multimap<double, pair<string, Phase::Type>>
        unmatchedPhases; // distance, <stationid,phaseType>

    auto eqlrng = catalog.getPhases().equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      const Phase &phase     = it->second;
      const Station &station = catalog.getStations().at(phase.stationId);

      // check pick weight
      if (phase.procInfo.weight < minPhaseWeight) continue;

      // check this station distance to reference event is ok
      const auto &staRefEvDistanceIt =
          validatedStationDistance.find(phase.stationId);
      if (staRefEvDistanceIt == validatedStationDistance.end()) continue;

      double staRefEvDistance = staRefEvDistanceIt->second;

      if ((staRefEvDistance / eventDistance) <
          minEStoIEratio) // ratio too small ?
      {
        continue;
      } // check if the station distance to the current event is valid
      if (maxESdist > 0 || minESdist > 0 || minEStoIEratio > 0)
      {
        // compute distance between current event and station
        double stationDistance = computeDistance(event, station);

        if ((maxESdist > 0 && stationDistance > maxESdist) || // too far away ?
            (stationDistance < minESdist) ||                  // too close ?
            ((stationDistance / eventDistance) <
             minEStoIEratio)) // ratio too small ?
        {
          continue;
        }
      }

      // now find corresponding phase in reference event phases
      bool peer_found = false;
      auto itRef      = refEvCatalog.searchPhase(refEv.id, phase.stationId,
                                                 phase.procInfo.type);
      if (itRef != refEvCatalog.getPhases().end())
      {
        const Phase &refPhase = itRef->second;
        if (refPhase.procInfo.weight >= minPhaseWeight) peer_found = true;
      }

      if (!peer_found)
      {
        unmatchedPhases.emplace(
            std::piecewise_construct, std::forward_as_tuple(staRefEvDistance),
            std::forward_as_tuple(phase.stationId, phase.procInfo.type));
        continue;
      }

      stationByDistance.emplace(
          std::piecewise_construct, std::forward_as_tuple(staRefEvDistance),
          std::forward_as_tuple(phase.stationId, phase.procInfo.type));
    }

    // check if enough phases (> `minDTperEvt`); if not skip event
    if (stationByDistance.size() < minDTperEvt)
    {
      continue;
    }

    unsigned numObservations = 0;

    // Since the constraints are met `evSelEntry` will be added to
    // `selectedEvents`
    SelectedEventEntry evSelEntry;
    evSelEntry.event = event;

    // copy the phases
    for (const auto &kw : stationByDistance)
    {
      // If `maxDTperEvt` is set, make sure to stay within limits.
      if (maxDTperEvt > 0 && numObservations >= maxDTperEvt) break;
      const pair<string, Phase::Type> &staPh = kw.second;
      evSelEntry.phases[staPh.first].insert(staPh.second);
      numObservations++;
    }

    if (keepUnmatched)
    {
      // Add theoretical picks to `refEv` of those neighbors' phases respecting
      // the constraints. Those phases might be used for cross-correlation.
      for (const auto &kw : unmatchedPhases)
      {
        // if `maxDTperEvt` is set, make sure to stay within limits
        if (maxDTperEvt > 0 && numObservations >= maxDTperEvt) break;
        const pair<string, Phase::Type> &staPh = kw.second;
        evSelEntry.phases[staPh.first].insert(staPh.second);
        numObservations++;
      }
    }

    // add this event to the selected ones
    selectedEvents.emplace(eventDistance, evSelEntry);
    dtCountByEvent.emplace(event.id, numObservations);

    // If ellipsoid algorithm is disabled and maxNumNeigh is set, then
    // we can stop when we have found maxNumNeigh events
    if (numEllipsoids <= 0 && maxNumNeigh > 0 &&
        selectedEvents.size() >= maxNumNeigh)
    {
      break;
    }
  }

  // Finally, build the catalog of neighboring events using either the
  // ellipsoids algorithm or the nearest neighbour method.
  unique_ptr<Neighbours> neighbours(new Neighbours());
  neighbours->refEvId = refEv.id;

  if (numEllipsoids <= 0)
  {
    // If `numEllipsoids` is 0, disable the ellipsoid algorithm and simply
    // select events by means of the nearest neighbor algorithm. Since
    // `selectedEvents` is sorted by distance, we obtain closer events first.
    for (const auto &kv : selectedEvents)
    {
      const SelectedEventEntry &evSelEntry = kv.second;
      const Event &ev                      = evSelEntry.event;

      // add this event to the catalog
      neighbours->ids.insert(ev.id);
      neighbours->phases.emplace(ev.id, evSelEntry.phases);

      logDebug("Neighbour: #phases %2d distance %g depth-diff %g event %s",
               dtCountByEvent[ev.id], distanceByEvent[ev.id],
               refEv.depth - ev.depth, string(ev).c_str());
    }
  }
  else
  {
    //
    // apply Waldauser's concentric ellipsoids algorithm
    //
    bool workToDo = true;

    while (workToDo)
    {
      for (int elpsNum = ellipsoids.size() - 1; elpsNum >= 0; elpsNum--)
      {
        for (int quadrant : {1, 2, 3, 4, 5, 6, 7, 8})
        {
          // if we either don't have events or we have already selected
          // `maxNumNeigh` neighbors, exit
          if (selectedEvents.empty() ||
              (maxNumNeigh > 0 && neighbours->ids.size() >= maxNumNeigh))
          {
            workToDo = false;
            break;
          }

          // since `selectedEvents` is sorted by distance, we get closer events
          // first
          for (auto it = selectedEvents.begin(); it != selectedEvents.end();
               it++)
          {
            const SelectedEventEntry &evSelEntry = it->second;
            const Event &ev                      = evSelEntry.event;

            bool found = ellipsoids[elpsNum].isInside(ev.latitude, ev.longitude,
                                                      ev.depth, quadrant);
            if (found)
            {
              // add this event to the catalog
              neighbours->ids.insert(ev.id);
              neighbours->phases.emplace(ev.id, evSelEntry.phases);

              logDebug("Neighbour: ellipsoid %2d quadrant %d #phases %2d "
                       "distance %5.2f depth-diff %6.3f event %s",
                       elpsNum, quadrant, dtCountByEvent[ev.id],
                       distanceByEvent[ev.id], refEv.depth - ev.depth,
                       string(ev).c_str());

              selectedEvents.erase(it);
              break;
            }
          }
        }
      }
    }
  }

  // check if enough neighbors were found
  if (neighbours->ids.size() < minNumNeigh)
  {
    string msg =
        strf("Skipping event %s, insufficient number of neighbors (%zu)",
             string(refEv).c_str(), neighbours->ids.size());
    logDebug("%s", msg.c_str());
    throw Exception(msg);
  }

  return neighbours;
}

unordered_map<unsigned, unique_ptr<Neighbours>>
selectNeighbouringEventsCatalog(const Catalog &catalog,
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
  logInfo("Selecting Catalog Neighbouring Events");

  // output: neighbours for each event in the catalog
  unordered_map<unsigned, unique_ptr<Neighbours>> neighboursList;

  // validCatalog contains events not discarded by user criteria
  Catalog validCatalog(catalog);

  // events discarded by user criteria
  unordered_set<unsigned> removedEvents;

  // for each event find its neighbours
  for (const auto &kv : catalog.getEvents())
  {
    const Catalog::Event& event = kv.second;

    unique_ptr<Neighbours> neighbours;
    try
    {
      neighbours = selectNeighbouringEvents(
          validCatalog, event, validCatalog, minPhaseWeight, minESdist,
          maxESdist, minEStoIEratio, minDTperEvt, maxDTperEvt, minNumNeigh,
          maxNumNeigh, numEllipsoids, maxEllipsoidSize, keepUnmatched);
    }
    catch (...)
    {}

    if (!neighbours)
    {
      // event discarded because it doesn't satisfy requirements
      removedEvents.insert(event.id);
      // we don't want other events to pick this as neighbour
      validCatalog.removeEvent(event.id);
      continue;
    }

    neighboursList.emplace(neighbours->refEvId, std::move(neighbours));
  }

  // if the removed events were used as neighbours of any other valid event
  // then try to rebuild their neighbours
  bool redo;
  do
  {
    logDebug("Found neighbours for %zu events (%zu events don't satisfy the constraints)",
            neighboursList.size(), removedEvents.size());

    logDebug("Fix events whose neighbours are the events not satisfying the constraints");

    redo = false;
    unordered_map<unsigned, unique_ptr<Neighbours>> validNeighbours;

    for (auto &kv : neighboursList)
    {
      unique_ptr<Neighbours> &neighbours = kv.second;
      bool invalid                       = false;

      // check if the neighbours are in the removed event list
      for (unsigned nbId : neighbours->ids)
      {
        if (removedEvents.count(nbId) != 0)
        {
          invalid = true;
          break;
        }
      }

      // if the current event uses at least a removed event, then rebuild
      // its neighbours
      if (invalid)
      {
        const Catalog::Event& event = validCatalog.getEvents().find(neighbours->refEvId)->second;

        try
        {
          neighbours = selectNeighbouringEvents(
              validCatalog, event, validCatalog, minPhaseWeight, minESdist,
              maxESdist, minEStoIEratio, minDTperEvt, maxDTperEvt, minNumNeigh,
              maxNumNeigh, numEllipsoids, maxEllipsoidSize, keepUnmatched);
        }
        catch (...)
        {}

        if (neighbours)
        {
          invalid = false;
        }
      }

      // failed to rebuild its neighbours: remove this event too
      if (invalid)
      {
        removedEvents.insert(neighbours->refEvId);
        validCatalog.removeEvent(neighbours->refEvId);
        redo = true;
        continue;
      }

      validNeighbours.emplace(neighbours->refEvId, std::move(neighbours));
    }

    neighboursList = std::move(validNeighbours);

  } while (redo);

  return neighboursList;
}

/*
 * Arrange neighbours in not connected clusters.
 *
 * Also, we don't want to report the same pair multiple times
 * (e.g. ev1-ev2 and ev2-ev1) since we only want one observation
 * per pair
 */
list<unordered_map<unsigned, unique_ptr<Neighbours>>>
clusterizeNeighbouringEvents(
    unordered_map<unsigned, unique_ptr<Neighbours>> &neighboursByEvent)
{
  map<unsigned, vector<unique_ptr<Neighbours>>> clusters; // cluster id, cluster
  unordered_map<unsigned, unsigned> clusterIdByEvent; // event id, cluster id

  logInfo("Searching for not connected clusters...");

  while (!neighboursByEvent.empty())
  {
    // keep track of event pairs found (useful to drop identical pairs)
    unordered_multimap<unsigned, unsigned> discoveredPairs;

    // start traversal with first unseen event neighbours
    unordered_set<unsigned> clusterEvs{neighboursByEvent.begin()->first};

    vector<unique_ptr<Neighbours>> currentCluster;
    unordered_set<unsigned> connectedClusters;

    while (!clusterEvs.empty()) // when empty, the cluster is fully built
    {
      unsigned currentEv = *clusterEvs.begin();
      clusterEvs.erase(currentEv);

      // keep track of clusters connected to the current one
      const auto &clusterIdByEventIt = clusterIdByEvent.find(currentEv);
      if (clusterIdByEventIt != clusterIdByEvent.end())
        connectedClusters.insert(clusterIdByEventIt->second);

      const auto &neighboursByEventIt = neighboursByEvent.find(currentEv);

      // skip already processed events
      if (neighboursByEventIt == neighboursByEvent.end()) continue;

      unique_ptr<Neighbours> neighbours =
          std::move(neighboursByEventIt->second);
      neighboursByEvent.erase(neighbours->refEvId);

      // update the set for the traversal of this cluster
      for (unsigned neighEvId : neighbours->ids) clusterEvs.insert(neighEvId);

      // remove from current neighbours the pairs that appeared previously
      // in this cluster
      auto eqlrng = discoveredPairs.equal_range(neighbours->refEvId);
      for (auto existingPair = eqlrng.first; existingPair != eqlrng.second;
           existingPair++)
      {
        neighbours->ids.erase(existingPair->second);
        neighbours->phases.erase(existingPair->second);
      }

      // keep track of new event pairs for following iterations and
      // update the queue for the breadth-first traversal
      for (unsigned neighEvId : neighbours->ids)
        discoveredPairs.emplace(neighEvId, neighbours->refEvId);

      // populate current cluster
      currentCluster.push_back(std::move(neighbours));
    }

    // merge all connected clusters to the current one
    for (unsigned clusterId : connectedClusters)
    {
      vector<unique_ptr<Neighbours>> &clstr = clusters.at(clusterId);
      currentCluster.insert(currentCluster.end(),
                            std::make_move_iterator(clstr.begin()),
                            std::make_move_iterator(clstr.end()));
      clusters.erase(clusterId);
    }

    // save current cluster
    unsigned maxKey       = clusters.empty() ? 0 : clusters.rbegin()->first;
    unsigned newClusterId = maxKey + 1;

    clusters.emplace(newClusterId, std::move(currentCluster));

    for (const unique_ptr<Neighbours> &n : clusters.at(newClusterId))
      clusterIdByEvent[n->refEvId] = newClusterId;
  }

  list<unordered_map<unsigned, unique_ptr<Neighbours>>> returnClusters;
  for (auto &kv : clusters)
  {
    vector<unique_ptr<Neighbours>> &currentCluster = kv.second;
    unordered_map<unsigned, unique_ptr<Neighbours>> convertedCluster;
    for (unique_ptr<Neighbours> &n : currentCluster)
      convertedCluster.emplace(n->refEvId, std::move(n));
    returnClusters.push_back(std::move(convertedCluster));
  }
  return returnClusters;
}

} // namespace HDD

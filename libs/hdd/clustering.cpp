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
#include "csvreader.h"
#include "ellipsoid.h"
#include "kdtree.h"
#include "log.h"
#include "utils.h"

using namespace std;
using namespace HDD::Logger;
using Event   = HDD::Catalog::Event;
using Phase   = HDD::Catalog::Phase;
using Station = HDD::Catalog::Station;

namespace HDD {

size_t Neighbours::amount() const { return _phases.size(); }

std::unordered_set<unsigned> Neighbours::ids() const
{
  std::unordered_set<unsigned> keys;
  keys.reserve(_phases.size());
  for (const auto &kv : _phases)
  {
    keys.insert(kv.first);
  }
  return keys;
}

void Neighbours::add(unsigned neighbourId,
                     const std::string &stationId,
                     const std::string &phase)
{
  _phases[neighbourId][phase].insert(stationId);
}

void Neighbours::remove(unsigned neighbourId) { _phases.erase(neighbourId); }

bool Neighbours::has(unsigned neighbourId) const
{
  return _phases.find(neighbourId) != _phases.end();
}

bool Neighbours::has(unsigned neighbourId, const std::string &stationId) const
{
  try
  {
    const auto &map = _phases.at(neighbourId);
    for (const auto &kv2 : map)
    {
      const string phase = kv2.first;
      const auto &set    = kv2.second;
      if (set.count(stationId) > 0)
      {
        return true;
      }
    }
  }
  catch (const std::out_of_range &e)
  {}
  return false;
}

bool Neighbours::has(unsigned neighbourId,
                     const std::string &stationId,
                     const std::string &phase) const
{
  try
  {
    return _phases.at(neighbourId).at(phase).count(stationId) > 0;
  }
  catch (const std::out_of_range &e)
  {
    return false;
  }
}

std::unordered_set<std::string> Neighbours::stations() const
{
  unordered_set<string> results;
  for (const auto &kv1 : _phases)
  {
    // unsigned neighbourId = kv1.first;
    const auto &map = kv1.second;
    for (const auto &kv2 : map)
    {
      const string phase = kv2.first;
      const auto &set    = kv2.second;
      for (const string &station : set)
      {
        results.insert(station);
      }
    }
  }
  return results;
}

std::vector<std::tuple<std::string, std::string, unsigned>>
Neighbours::phases() const
{
  vector<tuple<string, string, unsigned>> result;
  for (const auto &kv1 : _phases)
  {
    unsigned neighbourId = kv1.first;
    const auto &map      = kv1.second;
    for (const auto &kv2 : map)
    {
      const string phase = kv2.first;
      const auto &set    = kv2.second;
      for (const string &station : set)
      {
        result.emplace_back(station, phase, neighbourId);
      }
    }
  }
  return result;
}

std::vector<std::tuple<std::string, std::string>>
Neighbours::phases(unsigned neighbourId) const
{
  vector<tuple<string, string>> result;
  try
  {
    const auto &map = _phases.at(neighbourId);
    for (const auto &kv2 : map)
    {
      const string phase = kv2.first;
      const auto &set    = kv2.second;
      for (const string &station : set)
      {
        result.emplace_back(station, phase);
      }
    }
  }
  catch (const std::out_of_range &e)
  {}
  return result;
}

Catalog Neighbours::toCatalog(const Catalog &catalog, bool includeRefEv) const
{
  Catalog returnCat{};
  for (const auto &kv : _phases)
  {
    unsigned neighbourId = kv.first;
    returnCat.add(neighbourId, catalog, true);
  }
  if (includeRefEv)
  {
    returnCat.add(_refEvId, catalog, true);
  }
  return returnCat;
}

std::ofstream Neighbours::writeToFile(const Catalog &cat,
                                      const std::string &file) const
{
  ofstream os(file);
  os << "eventId1,eventId2,networkCode,stationCode,locationCode,phaseType"
     << endl;
  appendToStream(cat, os);
  return os;
}

void Neighbours::appendToStream(const Catalog &cat, std::ostream &os) const
{
  for (const auto &kv1 : _phases)
  {
    unsigned neighbourId = kv1.first;
    const auto &staPhs   = kv1.second;
    for (const auto &kv2 : staPhs)
    {
      const std::string phase                    = kv2.first;
      const std::unordered_set<string> &stations = kv2.second;
      for (const string &stationId : stations)
      {
        const Catalog::Station &sta = cat.getStations().at(stationId);
        os << strf("%u,%u,%s,%s,%s,%s", _refEvId, neighbourId,
                   sta.networkCode.c_str(), sta.stationCode.c_str(),
                   sta.locationCode.c_str(), phase.c_str())
           << endl;
      }
    }
  }
}

void Neighbours::writeToFile(
    const std::unordered_map<unsigned, Neighbours> &neighboursByEvent,
    const Catalog &cat,
    const std::string &file)
{
  ofstream os;
  bool first = true;

  for (const auto &kv1 : neighboursByEvent)
  {
    const Neighbours &n = kv1.second;
    if (first)
    {
      os    = n.writeToFile(cat, file);
      first = false;
    }
    else
    {
      n.appendToStream(cat, os);
    }
  }
}

unordered_map<unsigned, Neighbours>
Neighbours::readFromFile(const Catalog &cat, const std::string &file)
{
  unordered_map<unsigned, Neighbours> neighbours;

  vector<unordered_map<string, string>> lines = CSV::readWithHeader(file);
  int row_count                               = 0;
  try
  {
    for (const auto &row : lines)
    {
      row_count++;
      unsigned ev1        = std::stoul(row.at("eventId1"));
      unsigned ev2        = std::stoul(row.at("eventId2"));
      string networkCode  = row.at("networkCode");
      string stationCode  = row.at("stationCode");
      string locationCode = row.at("locationCode");
      string phase        = row.at("phaseType");

      std::string stationId =
          cat.searchStation(networkCode, stationCode, locationCode)->second.id;

      if (neighbours.count(ev1) <= 0)
      {
        neighbours.emplace(ev1, Neighbours(ev1));
      }

      if (neighbours.count(ev2) <= 0)
      {
        neighbours.emplace(ev2, Neighbours(ev2));
      }

      Neighbours &current = neighbours.find(ev1)->second;
      current.add(ev2, stationId, phase);
    }
  }
  catch (const std::exception &e)
  {
    string msg = strf("Error while parsing file '%s' at row %d: %s",
                      file.c_str(), row_count, e.what());
    throw Exception(msg);
  }
  return neighbours;
}

EventTree createEventTree(const Catalog &catalog)
{
  vector<EventTree::Point> points;
  points.reserve(catalog.getEvents().size());
  for (const auto &kv : catalog.getEvents())
  {
    const Event &event = kv.second;
    EventTree::Point point{event.latitude, event.longitude, event.depth,
                           event.id};
    points.push_back(std::move(point));
  }
  return EventTree(std::move(points));
}

Neighbours selectNeighbouringEvents(const EventTree &evTree,
                                    const Catalog &catalog,
                                    const Event &refEv,
                                    const Catalog &refEvCatalog,
                                    double minESdist,
                                    double maxESdist,
                                    double minEStoIEratio,
                                    unsigned minPhase,
                                    unsigned maxPhase,
                                    unsigned maxNumNeigh,
                                    unsigned numEllipsoids,
                                    double maxNeighbourDist)
{
  // Optimization: make code faster but the result will be the same.
  if (maxNumNeigh <= 0)
  {
    numEllipsoids = 0;
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
  // Helper that tries to add an event as neighbour to refEv if all the
  // constraints are met and returns true. Returns false otherwise
  //
  auto addNeighbour = [&](Neighbours &neighbours, unsigned eventId,
                          double eventDistance) -> bool {
    if (eventId == refEv.id)
    {
      return false;
    }

    auto search = catalog.getEvents().find(eventId);
    if (search == catalog.getEvents().end())
    {
      return false;
    }

    const Event &event = search->second;

    multimap<double, pair<string, string>>
        staPhByDistance; // distance, <stationid,phase>

    // Loop through event phases and keep track of the valid phases and their
    // station distance.
    auto eqlrng = catalog.getPhases().equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      const Phase &phase     = it->second;
      const Station &station = catalog.getStations().at(phase.stationId);

      // skip unwanted phase types
      if (phase.procInfo.type == Phase::Type::NO)
      {
        continue;
      }

      // check this station distance to reference event is ok
      const auto &staRefEvDistanceIt =
          validatedStationDistance.find(phase.stationId);
      if (staRefEvDistanceIt == validatedStationDistance.end())
      {
        continue;
      }

      double staRefEvDistance = staRefEvDistanceIt->second;

      // station distance to inter-event distance ratio too small ?
      if ((staRefEvDistance / eventDistance) < minEStoIEratio)
      {
        continue;
      }

      // check if the station distance to the current event is valid
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
      auto itRef =
          refEvCatalog.searchPhase(refEv.id, phase.stationId, phase.type);
      if (itRef == refEvCatalog.getPhases().end() ||
          itRef->second.procInfo.type == Phase::Type::NO)
      {
        // phase not found
        continue;
      }

      staPhByDistance.emplace(
          std::piecewise_construct, std::forward_as_tuple(staRefEvDistance),
          std::forward_as_tuple(phase.stationId, phase.type));
    }

    // check if enough phases (> `minPhase`); if not skip event
    if (staPhByDistance.size() < minPhase)
    {
      return false;
    }

    // The constraints are met so add the neighbour with the matched phases
    unsigned addedPhases = 0;
    for (const auto &kw : staPhByDistance)
    {
      // If `maxPhase` is set, make sure to stay within limits.
      if (maxPhase > 0 && addedPhases >= maxPhase)
      {
        break;
      }
      const pair<string, string> &staPh = kw.second; // station, phase
      neighbours.add(event.id, staPh.first, staPh.second);
      addedPhases++;
    }

    return true;
  };

  //
  // Find the neighboring events using either the ellipsoids algorithm or the
  // nearest neighbour method.
  //
  Neighbours neighbours(refEv.id);

  if (numEllipsoids <= 0) // nearest neighbour
  {

    auto callback = [&](const EventTree::Point &point, double eventDistance) {
      // we reached the search distance limit
      if (eventDistance > maxNeighbourDist)
      {
        return true; // stops bestFirstSearch
      }

      // add eventId as a valid neighbour
      const unsigned eventId = point.data;
      bool added             = addNeighbour(neighbours, eventId, eventDistance);

      // we have found maxNumNeigh events
      if (added && maxNumNeigh > 0 && neighbours.amount() >= maxNumNeigh)
      {
        return true; // stops bestFirstSearch
      }

      return false; // continue bestFirstSearch
    };

    // bestFirstSearch calls 'callback' by nearest neighbour until callback
    // retuns true
    evTree.bestFirstSearch(refEv.latitude, refEv.longitude, refEv.depth,
                           callback);
  }
  else // apply Waldauser's concentric ellipsoids algorithm
  {
    //
    // Build ellipsoids
    //
    vector<HddEllipsoid> ellipsoids(numEllipsoids, HddEllipsoid(0, 0, 0, 0));
    double verticalSemiAxisLen = maxNeighbourDist * 2;
    for (int i = numEllipsoids - 1; i > 0; --i)
    {
      ellipsoids[i] = HddEllipsoid(verticalSemiAxisLen, refEv.latitude,
                                   refEv.longitude, refEv.depth);
      verticalSemiAxisLen /= 2;
    }
    ellipsoids[0] = HddEllipsoid(0, verticalSemiAxisLen, refEv.latitude,
                                 refEv.longitude, refEv.depth);

    const HddEllipsoid &outmostEllip = ellipsoids.at(numEllipsoids - 1);

    //
    // Find all neighbours within radius equal to the max VerticalSemiAxisLen
    // of the largest ellipsoid
    //
    double radius = outmostEllip.getOuterEllipsoidVerticalSemiAxisLen();

    // distance km, point idx
    multimap<double, size_t> pointsByDistance = evTree.radiusSearch(
        refEv.latitude, refEv.longitude, refEv.depth, radius);

    // drop events outside the outmost ellipsod boundaries
    for (auto it = pointsByDistance.begin(); it != pointsByDistance.end();)
    {
      const size_t idx              = it->second;
      const EventTree::Point &point = evTree.points().at(idx);

      if (outmostEllip.getOuterEllipsoid().isInside(
              point.latitude, point.longitude, point.depth))
      {
        ++it;
        continue;
      }
      it = pointsByDistance.erase(it);
    }

    //
    // Select neighbours for each ellipsoid/quadrant combination in a round
    // robin fashion until 'maxNumNeigh' is reached
    //
    bool workToDo = !pointsByDistance.empty();

    while (workToDo)
    {
      bool match = false;

      // loop through ellipsoids
      for (const HddEllipsoid &ellpsd : ellipsoids)
      {
        auto minSearchDist = pointsByDistance.lower_bound(
            ellpsd.getInnerEllipsoidHorizontalSemiAxesLen());
        auto maxSearchDist = pointsByDistance.upper_bound(
            ellpsd.getOuterEllipsoidVerticalSemiAxisLen());

        // loop through quadrants
        for (const int quadrant : {1, 2, 3, 4, 5, 6, 7, 8})
        {

          // Search an event within the current ellipdoid/quadrant
          for (auto it = minSearchDist; it != maxSearchDist;)
          {
            const double eventDistance    = it->first;
            const size_t idx              = it->second;
            const EventTree::Point &point = evTree.points().at(idx);
            const unsigned eventId        = point.data;

            // check this event (point) is in the current ellipdoid/quadrant
            if (!ellpsd.isInside(point.latitude, point.longitude, point.depth,
                                 quadrant))
            {
              ++it;
              continue; // try next event
            }

            // event in current ellipdoid/quadrant -> try addinh it to the
            // neighbours
            bool added = addNeighbour(neighbours, eventId, eventDistance);

            //
            // remove the event from future selection
            //
            bool updateMinSearchDist = (minSearchDist == it);
            it                       = pointsByDistance.erase(it);
            if (updateMinSearchDist) minSearchDist = it;

            match = true;

            if (added) break; // next ellipsoid/quadrant
          }

          // if we either don't have events or we have already selected
          // `maxNumNeigh` neighbors, exit
          if (pointsByDistance.empty() ||
              (maxNumNeigh > 0 && neighbours.amount() >= maxNumNeigh))
          {
            workToDo = false;
            break;
          }
        }
        if (!workToDo) break;
      }

      if (!match && workToDo) // this should never happen, just safety belt
      {
        workToDo = false;
        logWarningF("Internal logic error in ellipsoid algorithm (remaining "
                    "events %zu)",
                    pointsByDistance.size());
      }
    }
  }

  return neighbours;
}

unordered_map<unsigned, Neighbours>
selectNeighbouringEventsCatalog(const EventTree &evTree,
                                const Catalog &catalog,
                                double minESdist,
                                double maxESdist,
                                double minEStoIEratio,
                                unsigned minPhase,
                                unsigned maxPhase,
                                unsigned minNumNeigh,
                                unsigned maxNumNeigh,
                                unsigned numEllipsoids,
                                double maxNeighbourDist)
{
  logInfo("Searching for event neighbours in the catalog");

  // output: neighbours for each event in the catalog
  unordered_map<unsigned, Neighbours> neighboursList;

  // validCatalog contains events not discarded by user criteria
  Catalog validCatalog(catalog);

  // events discarded by user criteria
  unordered_set<unsigned> removedEvents;

  // for each event find its neighbours
  for (const auto &kv : catalog.getEvents())
  {
    const Catalog::Event &event = kv.second;
    Neighbours neighbours       = selectNeighbouringEvents(
        evTree, validCatalog, event, validCatalog, minESdist, maxESdist,
        minEStoIEratio, minPhase, maxPhase, maxNumNeigh, numEllipsoids,
        maxNeighbourDist);

    if (neighbours.amount() >= minNumNeigh)
    {
      neighboursList.emplace(neighbours.referenceId(), std::move(neighbours));
    }
    else
    {
      // event discarded because it doesn't satisfy requirements
      removedEvents.insert(event.id);
      // we don't want other events to pick this as neighbour
      validCatalog.removeEvent(event.id);
    }
  }

  // if the removed events were used as neighbours of any other valid event
  // then try to rebuild their neighbours
  bool redo;
  do
  {
    redo = false;
    unordered_map<unsigned, Neighbours> validNeighbours;

    for (auto &kv : neighboursList)
    {
      Neighbours &neighbours = kv.second;
      bool invalid           = false;

      // check if the neighbours are in the removed event list
      for (unsigned nbId : neighbours.ids())
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
        const Catalog::Event &event =
            validCatalog.getEvents().find(neighbours.referenceId())->second;

        neighbours = selectNeighbouringEvents(
            evTree, validCatalog, event, validCatalog, minESdist, maxESdist,
            minEStoIEratio, minPhase, maxPhase, maxNumNeigh, numEllipsoids,
            maxNeighbourDist);

        invalid = (neighbours.amount() < minNumNeigh);
      }

      // failed to rebuild its neighbours: remove this event too
      if (invalid)
      {
        removedEvents.insert(neighbours.referenceId());
        validCatalog.removeEvent(neighbours.referenceId());
        redo = true;
        continue;
      }

      validNeighbours.emplace(neighbours.referenceId(), std::move(neighbours));
    }

    neighboursList = std::move(validNeighbours);

  } while (redo);

  return neighboursList;
}

/*
 * Organize the neighbours by not connected clusters.
 *
 * In addition, don't report the same pair multiple times (e.g. ev1-ev2 and
 * ev2-ev1) since we only need one observation per pair in the DD solver.
 *
 */
list<unordered_map<unsigned, Neighbours>> clusterizeNeighbouringEvents(
    unordered_map<unsigned, Neighbours> neighboursByEvent)
{
  map<unsigned, vector<Neighbours>> clusters;         // cluster id, cluster
  unordered_map<unsigned, unsigned> clusterIdByEvent; // event id, cluster id

  logInfoF("Forming clusters from the selected %zu events",
           neighboursByEvent.size());

  while (!neighboursByEvent.empty())
  {
    // keep track of event pairs found (useful to drop identical pairs)
    unordered_multimap<unsigned, unsigned> discoveredPairs;

    // start traversal with first unseen event neighbours
    unordered_set<unsigned> clusterEvs{neighboursByEvent.begin()->first};

    vector<Neighbours> currentCluster;
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

      Neighbours neighbours = std::move(neighboursByEventIt->second);
      neighboursByEvent.erase(neighbours.referenceId());

      // update the set for the traversal of this cluster
      for (unsigned neighEvId : neighbours.ids()) clusterEvs.insert(neighEvId);

      // remove from current neighbours the pairs that appeared previously
      // in this cluster
      auto eqlrng = discoveredPairs.equal_range(neighbours.referenceId());
      for (auto existingPair = eqlrng.first; existingPair != eqlrng.second;
           existingPair++)
      {
        neighbours.remove(existingPair->second);
      }

      // keep track of new event pairs for following iterations and
      // update the queue for the breadth-first traversal
      for (unsigned neighEvId : neighbours.ids())
        discoveredPairs.emplace(neighEvId, neighbours.referenceId());

      // populate current cluster
      currentCluster.push_back(std::move(neighbours));
    }

    // merge all connected clusters to the current one
    for (unsigned clusterId : connectedClusters)
    {
      vector<Neighbours> &clstr = clusters.at(clusterId);
      currentCluster.insert(currentCluster.end(),
                            std::make_move_iterator(clstr.begin()),
                            std::make_move_iterator(clstr.end()));
      clusters.erase(clusterId);
    }

    // save current cluster
    unsigned maxKey       = clusters.empty() ? 0 : clusters.rbegin()->first;
    unsigned newClusterId = maxKey + 1;

    clusters.emplace(newClusterId, std::move(currentCluster));

    for (const Neighbours &n : clusters.at(newClusterId))
    {
      clusterIdByEvent[n.referenceId()] = newClusterId;
    }
  }

  list<unordered_map<unsigned, Neighbours>> returnClusters;
  for (auto &kv : clusters)
  {
    vector<Neighbours> &currentCluster = kv.second;
    unordered_map<unsigned, Neighbours> convertedCluster;
    for (Neighbours &n : currentCluster)
    {
      convertedCluster.emplace(n.referenceId(), std::move(n));
    }
    returnClusters.push_back(std::move(convertedCluster));
  }
  return returnClusters;
}

} // namespace HDD

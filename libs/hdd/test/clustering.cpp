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

#define SEISCOMP_TEST_MODULE hdd
#include <seiscomp/unittest/unittests.h>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>

#include "hdd/catalog.h"
#include "hdd/clustering.h"
#include "hdd/utils.h"

using namespace std;
using namespace HDD;
using Event     = HDD::Catalog::Event;
using Phase     = HDD::Catalog::Phase;
using Station   = HDD::Catalog::Station;
namespace bdata = boost::unit_test::data;

namespace {

// 4 stations at `distance` km N, E, S and W of lat/lon. They have slightly
// different elevations so that their distance from an event is never the same:
// this makes the order in which the stations are selected deterministic (the
// closest first, ST01 then ST02...)
void addStationsToCatalog(HDD::Catalog &cat,
                          double lat,
                          double lon,
                          double distance)
{
  Station sta;
  double staLat, staLon;

  computeCoordinates(distance, degToRad(0), lat, lon, staLat, staLon);
  sta = {"NET.ST01", staLat, staLon, 250, "NET", "ST01", ""};
  cat.addStation(sta);
  computeCoordinates(distance, degToRad(90), lat, lon, staLat, staLon);
  sta = {"NET.ST02", staLat, staLon, 295, "NET", "ST02", ""};
  cat.addStation(sta);
  computeCoordinates(distance, degToRad(180), lat, lon, staLat, staLon);
  sta = {"NET.ST03", staLat, staLon, 301, "NET", "ST03", ""};
  cat.addStation(sta);
  computeCoordinates(distance, degToRad(270), lat, lon, staLat, staLon);
  sta = {"NET.ST04", staLat, staLon, 395, "NET", "ST04", ""};
  cat.addStation(sta);
}

void addEventToCatalog(HDD::Catalog &cat, double lat, double lon, double depth)
{
  Event ev{};
  ev.latitude            = lat;
  ev.longitude           = lon;
  ev.depth               = depth;
  const unsigned eventId = cat.addEvent(ev, false);

  for (const auto &kv : cat.getStations())
  {
    const Station &sta = kv.second;

    Phase ph;
    ph.eventId          = eventId;
    ph.stationId        = sta.id;
    ph.lowerUncertainty = 0;
    ph.upperUncertainty = 0;
    ph.type             = "P";
    ph.networkCode      = sta.networkCode;
    ph.stationCode      = sta.stationCode;
    ph.locationCode     = sta.locationCode;
    ph.channelCode      = "";
    cat.addPhase(ph);

    ph.type = "S";
    cat.addPhase(ph);
  }
}

void addNeighboursToCatalog(HDD::Catalog &cat,
                            const Event &event,
                            const vector<double> &neighDist)
{
  for (double distance : neighDist)
  {
    double depth = event.depth;

    double neighbourLat, neighbourLon;
    computeCoordinates(distance, degToRad(45), event.latitude, event.longitude,
                       neighbourLat, neighbourLon, depth);
    addEventToCatalog(cat, neighbourLat, neighbourLon, depth);
    computeCoordinates(distance, degToRad(135), event.latitude, event.longitude,
                       neighbourLat, neighbourLon, depth);
    addEventToCatalog(cat, neighbourLat, neighbourLon, depth);
    computeCoordinates(distance, degToRad(225), event.latitude, event.longitude,
                       neighbourLat, neighbourLon, depth);
    addEventToCatalog(cat, neighbourLat, neighbourLon, depth);
    computeCoordinates(distance, degToRad(315), event.latitude, event.longitude,
                       neighbourLat, neighbourLon, depth);
    addEventToCatalog(cat, neighbourLat, neighbourLon, depth);
  }
}

HDD::Catalog buildCatalog(double lat,
                          double lon,
                          double depth,
                          double stationsDistance,
                          const vector<double> &neighDist,
                          const std::vector<std::string> &PphaseToKeep,
                          const std::vector<std::string> &SphaseToKeep)
{
  HDD::Catalog cat;
  addStationsToCatalog(cat, lat, lon, stationsDistance);
  addEventToCatalog(cat, lat, lon, depth);
  const Event &event = cat.getEvents().begin()->second;
  addNeighboursToCatalog(cat, event, neighDist);
  return HDD::Catalog::fillProcessingInfo(
      cat, Phase::Source::CATALOG, PphaseToKeep, SphaseToKeep, {});
}

void checkNeighbours(const HDD::Neighbours &neighbours,
                     const unordered_set<unsigned> &mustBePresentIds,
                     const unordered_set<string> &mustBePresentStations,
                     const unordered_set<string> &mustBePresentPhases)
{
  BOOST_CHECK(mustBePresentStations == neighbours.stations());
  for (const unsigned id : mustBePresentIds)
  {
    for (const string &sta : mustBePresentStations)
    {
      BOOST_CHECK(neighbours.has(id, sta));
      for (const string &pht : mustBePresentPhases)
      {
        BOOST_CHECK(neighbours.has(id, sta, pht));
      }
    }
  }
}

// Check that `n` contains the neighbours of event `refId` with exactly the
// expected (station id, phase, neighbour id) entries
void checkNeighbourPhases(
    const HDD::Neighbours &n,
    unsigned refId,
    vector<tuple<string, string, unsigned>> expectedPhases)
{
  unordered_set<unsigned> expectedIds;
  unordered_set<string> expectedStations;
  for (const auto &t : expectedPhases)
  {
    expectedStations.insert(std::get<0>(t));
    expectedIds.insert(std::get<2>(t));
  }

  BOOST_CHECK(n.referenceId() == refId);
  BOOST_CHECK(n.ids() == expectedIds);
  BOOST_CHECK(n.stations() == expectedStations);

  vector<tuple<string, string, unsigned>> phases = n.phases();
  std::sort(expectedPhases.begin(), expectedPhases.end());
  std::sort(phases.begin(), phases.end());
  BOOST_CHECK(phases == expectedPhases);

  for (const auto &t : expectedPhases)
  {
    const string &stationId  = std::get<0>(t);
    const string &phase      = std::get<1>(t);
    const unsigned neighEvId = std::get<2>(t);
    BOOST_CHECK(n.has(neighEvId, stationId));
    BOOST_CHECK(n.has(neighEvId, stationId, phase));
  }
}

Neighbours selectNeighbouringEventsThrow(const EventTree &evTree,
                                         const Catalog &catalog,
                                         const Event &refEv,
                                         const Catalog &refEvCatalog,
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
  Neighbours neighbours = HDD::selectNeighbouringEvents(
      evTree, catalog, refEv, refEvCatalog, minESdist, maxESdist,
      minEStoIEratio, minPhase, maxPhase, maxNumNeigh, numEllipsoids,
      maxNeighbourDist);

  // check if enough neighbors were found
  if (neighbours.amount() < minNumNeigh)
  {
    throw HDD::Exception("Insufficient number of neighbors");
  }

  return neighbours;
}

struct Origin
{
  double lat;
  double lon;
  double depth;
};

const vector<Origin> orgList{
    {0, 0, 1},       {45, 90, 3},     {45, -90, 3},     {-45, 90, 3},
    {-45, -90, 3},   {20, 180, 8},    {-20, 180, 8},    {20, -180, 8},
    {-20, -180, 8},  {30, 179.99, 5}, {-30, 179.99, 5}, {-30, -179.99, 5},
    {30, -179.99, 5}};

} // namespace

// Neighbour selection for a single event with 4 neighbours at 1 km and 4
// stations at 50 km: verify that each selection constraint (station distance,
// station to inter-event distance ratio, number of phases and neighbours,
// neighbour distance) rejects the neighbours when it cannot be satisfied and
// limits them as expected otherwise. All tests are repeated at the locations
// in `orgList`, which includes the poles and the antimeridian.
BOOST_DATA_TEST_CASE(test_clustering1, bdata::xrange(orgList.size()), orgIdx)
{
  const Origin &org = orgList.at(orgIdx);
  const HDD::Catalog cat =
      buildCatalog(org.lat, org.lon, org.depth, 50, {1}, {"P"}, {});
  const EventTree tree = createEventTree(cat);


  const Event &event = cat.getEvents().at(1);
  const unordered_set<unsigned> neighbourIds{2, 3, 4, 5};

  // neighbours out of range: maxNeighbourDist
  BOOST_CHECK_THROW(
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh   0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    0.5)   //  double maxNeighbourDist
      ,
      HDD::Exception);

  // stations out of range: minESdis
  BOOST_CHECK_THROW(
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    52, //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2)  //  double maxNeighbourDist
      ,
      HDD::Exception);

  // stations out of range: maxESdis
  BOOST_CHECK_THROW(
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    48, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2)  //  double maxNeighbourDist
      ,
      HDD::Exception);

  // minEStoIEratio not enough
  BOOST_CHECK_THROW(
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    60, //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2)  //  double maxNeighbourDist
      ,
      HDD::Exception);

  // too many DTs per event requested
  BOOST_CHECK_THROW(
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    5,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2)  //  double maxNeighbourDist
      ,
      HDD::Exception);

  // too many neighbours requested
  BOOST_CHECK_THROW(
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    5,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2)  //  double maxNeighbourDist
      ,
      HDD::Exception);

  // maxNeighbourDist/numEllipsoids
  HDD::Neighbours neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    0,  //  unsigned numEllipsoids
                                    2);  //  double maxNeighbourDist

  BOOST_CHECK_EQUAL(neighbours.referenceId(), event.id);
  BOOST_CHECK(neighbours.ids() == neighbourIds);
  checkNeighbours(neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."}, {"P"});

  // maxNeighbourDist/numEllipsoids
  neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2); //  double maxNeighbourDist

  BOOST_CHECK_EQUAL(neighbours.referenceId(), event.id);
  BOOST_CHECK(neighbours.ids() == neighbourIds);
  checkNeighbours(neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."}, {"P"});

  // maxNeighbourDist/numEllipsoids
  neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    5,  //  unsigned numEllipsoids
                                    2); //  double maxNeighbourDist

  BOOST_CHECK_EQUAL(neighbours.referenceId(), event.id);
  BOOST_CHECK(neighbours.ids() == neighbourIds);
  checkNeighbours(neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."}, {"P"});

  // minEStoIEratio
  neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    48, //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh   0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2); //  double maxNeighbourDist

  BOOST_CHECK_EQUAL(neighbours.referenceId(), event.id);
  BOOST_CHECK(neighbours.ids() == neighbourIds);
  checkNeighbours(neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."}, {"P"});

  // maxNumNeigh
  neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    2,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2); //  double maxNeighbourDist

  BOOST_CHECK_EQUAL(neighbours.referenceId(), event.id);
  BOOST_CHECK(neighbours.ids().size() == 2);
  BOOST_CHECK(neighbours.amount() == 2);

  // maxNumNeigh
  neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    3,  //  unsigned maxNumNeigh   0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2); //  double maxNeighbourDist

  BOOST_CHECK_EQUAL(neighbours.referenceId(), event.id);
  BOOST_CHECK(neighbours.ids().size() == 3);
  BOOST_CHECK(neighbours.amount() == 3);

  // minDTperEvt
  neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    4,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    5,  //  unsigned numEllipsoids
                                    2); //  double maxNeighbourDist

  BOOST_CHECK_EQUAL(neighbours.referenceId(), event.id);
  BOOST_CHECK(neighbours.ids() == neighbourIds);
  checkNeighbours(neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."}, {"P"});

  // maxDTperEvt
  neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    1,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2); //  double maxNeighbourDist

  BOOST_CHECK_EQUAL(neighbours.referenceId(), event.id);
  BOOST_CHECK(neighbours.ids() == neighbourIds);
  checkNeighbours(neighbours, neighbourIds, {"NET.ST01."}, {"P"});

  // maxDTperEvt
  neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    2,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2); //  double maxNeighbourDist

  BOOST_CHECK_EQUAL(neighbours.referenceId(), event.id);
  BOOST_CHECK(neighbours.ids() == neighbourIds);
  checkNeighbours(neighbours, neighbourIds, {"NET.ST01.", "NET.ST02."}, {"P"});
}

// Neighbour selection with the ellipsoid algorithm, for a single event with
// neighbours at 0.5, 0.6, 1.5, 1.6, 3.0 and 3.1 km in 4 directions (two events
// at 0.5, 1.5 and 3.0 km): verify which neighbours are selected for different
// numbers of neighbours and constraints
BOOST_DATA_TEST_CASE(test_clustering2, bdata::xrange(orgList.size()), orgIdx)
{
  const Origin &org = orgList.at(orgIdx);
  const HDD::Catalog cat =
      buildCatalog(org.lat, org.lon, org.depth, 70,
                   {0.5, 0.5, 0.6, 1.5, 1.5, 1.6, 3.0, 3.0, 3.1}, {}, {"S"});
  const EventTree tree = createEventTree(cat);


  const Event &event = cat.getEvents().at(1);
  unordered_set<unsigned> neighbourIds;

  // 3 ellipsoidal layers, up to 24 neighbours: in a round-robin over layers
  // and quadrants the closest events are selected, which are the events at
  // 0.5, 1.5 and 3.0 km, while those at 0.6, 1.6 and 3.1 km are left out
  HDD::Neighbours neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    24, //  unsigned maxNumNeigh  0 = no limits
                                    3,  //  unsigned numEllipsoids
                                    4); //  double maxNeighbourDist

  neighbourIds = {2,  3,  4,  5,  6,  7,  8,  9,  14, 15, 16, 17,
                  18, 19, 20, 21, 26, 27, 28, 29, 30, 31, 32, 33};

  BOOST_CHECK_EQUAL(neighbours.referenceId(), event.id);
  BOOST_CHECK(neighbours.ids() == neighbourIds);
  checkNeighbours(neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},{"S"});

  // As above, with station distance and phase constraints and exactly 16
  // neighbours: only the events at 0.5 and 1.5 km are selected
  neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    50, //  double minESdis
                                    91, //  double maxESdis      -1 = no limits
                                    10, //  double minEStoIEratio
                                    4,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    16, //  unsigned minNumNeigh
                                    16, //  unsigned maxNumNeigh  0 = no limits
                                    3,  //  unsigned numEllipsoids
                                    4); //  double maxNeighbourDist

  neighbourIds = {2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21};

  BOOST_CHECK_EQUAL(neighbours.referenceId(), event.id);
  BOOST_CHECK(neighbours.ids() == neighbourIds);
  checkNeighbours(neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."}, {"S"});

  // As above, but at most 2 phases per event pair: only the 2 closest
  // stations are used
  neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    50, //  double minESdis
                                    91, //  double maxESdis      -1 = no limits
                                    10, //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    2,  //  unsigned maxDTperEvt  0 = no limits
                                    16, //  unsigned minNumNeigh
                                    16, //  unsigned maxNumNeigh  0 = no limits
                                    3,  //  unsigned numEllipsoids
                                    4); //  double maxNeighbourDist

  neighbourIds = {2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21};

  BOOST_CHECK_EQUAL(neighbours.referenceId(), event.id);
  BOOST_CHECK(neighbours.ids() == neighbourIds);
  checkNeighbours(neighbours, neighbourIds, {"NET.ST01.", "NET.ST02."}, {"S"});

  // As above, with exactly 8 neighbours: only the events at 0.5 km are
  // selected
  neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    50, //  double minESdis
                                    91, //  double maxESdis      -1 = no limits
                                    10, //  double minEStoIEratio
                                    4,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    8,  //  unsigned minNumNeigh
                                    8,  //  unsigned maxNumNeigh  0 = no limits
                                    3,  //  unsigned numEllipsoids
                                    4);  //  double maxNeighbourDist

  neighbourIds = {2, 3, 4, 5, 6, 7, 8, 9};

  BOOST_CHECK_EQUAL(neighbours.referenceId(), event.id);
  BOOST_CHECK(neighbours.ids() == neighbourIds);
  checkNeighbours(neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."}, {"S"});
}

// Neighbour selection for all events of a catalog and clustering: a central
// event with 4 events at 7.5 km and 4 at 15 km (max neighbour distance 8 km).
// Verify the neighbours of each event and that all events form a single
// cluster, where each pair is kept only once
BOOST_DATA_TEST_CASE(test_clustering3, bdata::xrange(orgList.size()), orgIdx)
{
  const Origin &org = orgList.at(orgIdx);
  const HDD::Catalog cat =
      buildCatalog(org.lat, org.lon, org.depth, 100, {7.5, 15}, {"P"}, {"S"});
  const EventTree tree = createEventTree(cat);


  // find Neighbours for each event in the catalog
  unordered_map<unsigned, HDD::Neighbours> neighboursByEvent =
      HDD::selectNeighbouringEventsCatalog(
          tree, cat,
          0,      //  double minESdis
          -1,     //  double maxESdis      -1 = no limits
          0,      //  double minEStoIEratio
          1,      //  unsigned minDTperEvt
          0,      //  unsigned maxDTperEvt  0 = no limits
          1,      //  unsigned minNumNeigh
          0,      //  unsigned maxNumNeigh   0 = no limits
          0,      //  unsigned numEllipsoids
          8);     //  double maxNeighbourDist

  BOOST_REQUIRE(neighboursByEvent.size() == 9);

  const unordered_map<unsigned, unordered_set<unsigned>> expectedNeighbours = {
      {1, {2, 3, 4, 5}}, {2, {1, 6}}, {3, {1, 7}}, {4, {1, 8}}, {5, {1, 9}},
      {6, {2}},          {7, {3}},    {8, {4}},    {9, {5}}};

  for (const auto &kv : neighboursByEvent)
  {
    const HDD::Neighbours &neighbours = kv.second;
    BOOST_REQUIRE(neighbours.referenceId() == kv.first);
    BOOST_REQUIRE_NO_THROW(expectedNeighbours.at(neighbours.referenceId()));
    const unordered_set<unsigned> &neighbourIds =
        expectedNeighbours.at(neighbours.referenceId());
    BOOST_REQUIRE(neighbours.ids() == neighbourIds);
    checkNeighbours(neighbours, neighbourIds,
                    {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                    {"P", "S"});
  }

  list<unordered_map<unsigned, HDD::Neighbours>> clusters =
      HDD::clusterizeNeighbouringEvents(std::move(neighboursByEvent));

  BOOST_REQUIRE(clusters.size() == 1);

  for (const auto &cluster : clusters)
    for (const auto &kv : cluster)
    {
      const HDD::Neighbours &neighbours = kv.second;
      BOOST_REQUIRE(neighbours.referenceId() == kv.first);
      BOOST_REQUIRE_NO_THROW(expectedNeighbours.at(neighbours.referenceId()));
      const unordered_set<unsigned> &neighbourIds =
          expectedNeighbours.at(neighbours.referenceId());
      for (const unsigned neighId : neighbourIds)
      {
        if (neighbours.ids().count(neighId) == 1) // is in our Neighbours
        {
          checkNeighbours(neighbours, {neighId},
                          {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                          {"P", "S"});
        }
        else // must be in the other event's Neighbours
        {
          const HDD::Neighbours &otherNeighbours = cluster.at(neighId);
          checkNeighbours(otherNeighbours, {neighbours.referenceId()},
                          {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                          {"P", "S"});
        }
      }
    }
}

// As test_clustering3, with events at 3 and 4 km and max neighbour distance
// 1.5 km: the central event has no neighbours and is discarded, the other
// events form 4 independent clusters of 2 events
BOOST_DATA_TEST_CASE(test_clustering4, bdata::xrange(orgList.size()), orgIdx)
{
  const Origin &org = orgList.at(orgIdx);
  const HDD::Catalog cat =
      buildCatalog(org.lat, org.lon, org.depth, 10, {3, 4}, {}, {"S"});
  const EventTree tree = createEventTree(cat);


  // find Neighbours for each event in the catalog
  unordered_map<unsigned, HDD::Neighbours> neighboursByEvent =
      HDD::selectNeighbouringEventsCatalog(
          tree, cat,
          0,      //  double minESdis
          -1,     //  double maxESdis      -1 = no limits
          0,      //  double minEStoIEratio
          1,      //  unsigned minDTperEvt
          0,      //  unsigned maxDTperEvt  0 = no limits
          1,      //  unsigned minNumNeigh
          0,      //  unsigned maxNumNeigh   0 = no limits
          0,      //  unsigned numEllipsoids
          1.5);   //  double maxNeighbourDist

  BOOST_REQUIRE(neighboursByEvent.size() == 8);

  const unordered_map<unsigned, unordered_set<unsigned>> expectedNeighbours = {
      {2, {6}}, {3, {7}}, {4, {8}}, {5, {9}},
      {6, {2}}, {7, {3}}, {8, {4}}, {9, {5}}};

  for (const auto &kv : neighboursByEvent)
  {
    const HDD::Neighbours &neighbours = kv.second;
    BOOST_REQUIRE(neighbours.referenceId() == kv.first);
    BOOST_REQUIRE_NO_THROW(expectedNeighbours.at(neighbours.referenceId()));
    const unordered_set<unsigned> &neighbourIds =
        expectedNeighbours.at(neighbours.referenceId());
    BOOST_REQUIRE(neighbours.ids() == neighbourIds);
    checkNeighbours(neighbours, neighbourIds,
                    {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                    {"S"});
  }

  list<unordered_map<unsigned, HDD::Neighbours>> clusters =
      HDD::clusterizeNeighbouringEvents(std::move(neighboursByEvent));

  BOOST_REQUIRE(clusters.size() == 4);

  for (const auto &cluster : clusters)
    for (const auto &kv : cluster)
    {
      const HDD::Neighbours &neighbours = kv.second;
      BOOST_REQUIRE(neighbours.referenceId() == kv.first);
      BOOST_REQUIRE_NO_THROW(expectedNeighbours.at(neighbours.referenceId()));
      const unordered_set<unsigned> &neighbourIds =
          expectedNeighbours.at(neighbours.referenceId());
      for (const unsigned neighId : neighbourIds)
      {
        if (neighbours.ids().count(neighId) == 1) // is in our Neighbours
        {
          checkNeighbours(neighbours, {neighId},
                          {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                          {"S"});
        }
        else // must be in the other event's Neighbours
        {
          const HDD::Neighbours &otherNeighbours = cluster.at(neighId);
          checkNeighbours(otherNeighbours, {neighbours.referenceId()},
                          {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                          {"S"});
        }
      }
    }
}

// Phase matching: event pairs use only the phases with the same name at the
// same station (e.g. Pg with Pg, not with P or Pn). Four events 10 km apart on
// a line, with different phase names: verify the phases of each pair, before
// and after clustering
BOOST_DATA_TEST_CASE(test_clustering5, bdata::xrange(orgList.size()), orgIdx)
{
  const Origin &org = orgList.at(orgIdx);

  UTCTime time{}; // time is not relevant in clustering

  std::vector<HDD::Catalog::Station> stations = {
      {"N1.STA01.", 0.0, 0.0, 0.0, "N1", "STA01", ""},
      {"N1.STA02.", 0.0, 0.0, 0.0, "N1", "STA02", ""},
  };

  std::vector<HDD::Catalog::Event> events = {
      {0, time, 0.0, 0.0, 0.0, 0.0},
      {1, time, 0.0, 0.0, 0.0, 0.0},
      {2, time, 0.0, 0.0, 0.0, 0.0},
      {3, time, 0.0, 0.0, 0.0, 0.0},
  };

  const std::vector<HDD::Catalog::Phase> phases = {
      {0, "N1.STA01.", time, 0, 0, "P",  "N1", "STA01", "", "XY1"},
      {0, "N1.STA01.", time, 0, 0, "S",  "N1", "STA01", "", "XY1"},
      {0, "N1.STA01.", time, 0, 0, "Pg", "N1", "STA01", "", "XY1"},
      {0, "N1.STA01.", time, 0, 0, "Sg", "N1", "STA01", "", "XY1"},
      {0, "N1.STA02.", time, 0, 0, "P",  "N1", "STA02", "", "XY1"},
      {0, "N1.STA02.", time, 0, 0, "S",  "N1", "STA02", "", "XY1"},
      {0, "N1.STA02.", time, 0, 0, "Pg", "N1", "STA02", "", "XY1"},
      {0, "N1.STA02.", time, 0, 0, "Sg", "N1", "STA02", "", "XY1"},
      {1, "N1.STA01.", time, 0, 0, "Pn", "N1", "STA01", "", "XY1"},
      {1, "N1.STA01.", time, 0, 0, "Sn", "N1", "STA01", "", "XY1"},
      {1, "N1.STA01.", time, 0, 0, "Pg", "N1", "STA01", "", "XY1"},
      {1, "N1.STA01.", time, 0, 0, "Sg", "N1", "STA01", "", "XY1"},
      {1, "N1.STA02.", time, 0, 0, "Pn", "N1", "STA02", "", "XY1"},
      {1, "N1.STA02.", time, 0, 0, "Sn", "N1", "STA02", "", "XY1"},
      {1, "N1.STA02.", time, 0, 0, "Pg", "N1", "STA02", "", "XY1"},
      {1, "N1.STA02.", time, 0, 0, "Sg", "N1", "STA02", "", "XY1"},
      {2, "N1.STA01.", time, 0, 0, "Pn", "N1", "STA01", "", "XY1"},
      {2, "N1.STA01.", time, 0, 0, "Sn", "N1", "STA01", "", "XY1"},
      {2, "N1.STA01.", time, 0, 0, "P",  "N1", "STA01", "", "XY1"},
      {2, "N1.STA01.", time, 0, 0, "S",  "N1", "STA01", "", "XY1"},
      {2, "N1.STA02.", time, 0, 0, "Pn", "N1", "STA02", "", "XY1"},
      {2, "N1.STA02.", time, 0, 0, "Sn", "N1", "STA02", "", "XY1"},
      {2, "N1.STA02.", time, 0, 0, "P",  "N1", "STA02", "", "XY1"},
      {2, "N1.STA02.", time, 0, 0, "S",  "N1", "STA02", "", "XY1"},
      {3, "N1.STA01.", time, 0, 0, "P",  "N1", "STA01", "", "XY1"},
      {3, "N1.STA01.", time, 0, 0, "S",  "N1", "STA01", "", "XY1"},
      {3, "N1.STA01.", time, 0, 0, "Pg", "N1", "STA01", "", "XY1"},
      {3, "N1.STA01.", time, 0, 0, "Sg", "N1", "STA01", "", "XY1"},
      {3, "N1.STA02.", time, 0, 0, "P",  "N1", "STA02", "", "XY1"},
      {3, "N1.STA02.", time, 0, 0, "S",  "N1", "STA02", "", "XY1"},
      {3, "N1.STA02.", time, 0, 0, "Pg", "N1", "STA02", "", "XY1"},
      {3, "N1.STA02.", time, 0, 0, "Sg", "N1", "STA02", "", "XY1"},
  };

  double refLat = org.lat;
  double refLon = org.lon;

  computeCoordinates(10., degToRad(0), refLat, refLon,
      stations[0].latitude, stations[0].longitude);
  computeCoordinates(10., degToRad(180), refLat, refLon,
      stations[1].latitude, stations[1].longitude);

  computeCoordinates(15., degToRad(90),  refLat, refLon,
      events[0].latitude, events[0].longitude);
  computeCoordinates(5.,  degToRad(90),  refLat, refLon,
      events[1].latitude, events[1].longitude);
  computeCoordinates(5.,  degToRad(270), refLat, refLon,
      events[2].latitude, events[2].longitude);
  computeCoordinates(15., degToRad(270), refLat, refLon,
      events[3].latitude, events[3].longitude);

  const HDD::Catalog cat = HDD::Catalog::fillProcessingInfo(
       HDD::Catalog(stations, events, phases), Phase::Source::CATALOG,
       {"P", "Pg", "Pn"}, {"S", "Sg", "Sn"}, {});
  const EventTree tree = createEventTree(cat);

  // find Neighbours for each event in the catalog
  const unordered_map<unsigned, HDD::Neighbours> neighboursByEvent =
      HDD::selectNeighbouringEventsCatalog(
          tree, cat,
          0,      //  double minESdis
          -1,     //  double maxESdis      -1 = no limits
          0,      //  double minEStoIEratio
          1,      //  unsigned minDTperEvt
          0,      //  unsigned maxDTperEvt  0 = no limits
          1,      //  unsigned minNumNeigh
          0,      //  unsigned maxNumNeigh   0 = no limits
          0,      //  unsigned numEllipsoids
       10.1);     //  double maxNeighbourDist

  BOOST_REQUIRE(neighboursByEvent.size() == 4);

  // each event is paired with the events 10 km away, on the common phases
  checkNeighbourPhases(neighboursByEvent.at(0), 0,
                       {{"N1.STA01.", "Pg", 1},
                        {"N1.STA01.", "Sg", 1},
                        {"N1.STA02.", "Pg", 1},
                        {"N1.STA02.", "Sg", 1}});
  checkNeighbourPhases(neighboursByEvent.at(1), 1,
                       {{"N1.STA01.", "Pg", 0},
                        {"N1.STA01.", "Sg", 0},
                        {"N1.STA02.", "Pg", 0},
                        {"N1.STA02.", "Sg", 0},
                        {"N1.STA01.", "Pn", 2},
                        {"N1.STA01.", "Sn", 2},
                        {"N1.STA02.", "Pn", 2},
                        {"N1.STA02.", "Sn", 2}});
  checkNeighbourPhases(neighboursByEvent.at(2), 2,
                       {{"N1.STA01.", "Pn", 1},
                        {"N1.STA01.", "Sn", 1},
                        {"N1.STA02.", "Pn", 1},
                        {"N1.STA02.", "Sn", 1},
                        {"N1.STA01.", "P", 3},
                        {"N1.STA01.", "S", 3},
                        {"N1.STA02.", "P", 3},
                        {"N1.STA02.", "S", 3}});
  checkNeighbourPhases(neighboursByEvent.at(3), 3,
                       {{"N1.STA01.", "P", 2},
                        {"N1.STA01.", "S", 2},
                        {"N1.STA02.", "P", 2},
                        {"N1.STA02.", "S", 2}});

  const list<unordered_map<unsigned, HDD::Neighbours>> clusters =
      HDD::clusterizeNeighbouringEvents(neighboursByEvent);

  BOOST_REQUIRE(clusters.size() == 1);

  const unordered_map<unsigned, HDD::Neighbours> &cluster = clusters.front();
  BOOST_CHECK(cluster.size() == neighboursByEvent.size());

  // each pair is kept only once, by one of its two events
  checkNeighbourPhases(cluster.at(0), 0,
                       {{"N1.STA01.", "Pg", 1},
                        {"N1.STA01.", "Sg", 1},
                        {"N1.STA02.", "Pg", 1},
                        {"N1.STA02.", "Sg", 1}});
  checkNeighbourPhases(cluster.at(1), 1,
                       {{"N1.STA01.", "Pn", 2},
                        {"N1.STA01.", "Sn", 2},
                        {"N1.STA02.", "Pn", 2},
                        {"N1.STA02.", "Sn", 2}});
  checkNeighbourPhases(cluster.at(2), 2,
                       {{"N1.STA01.", "P", 3},
                        {"N1.STA01.", "S", 3},
                        {"N1.STA02.", "P", 3},
                        {"N1.STA02.", "S", 3}});
  checkNeighbourPhases(cluster.at(3), 3, {});
}

// Verify that the station distance limits apply to the distance between the
// stations and the reference event too, not only to the distance between the
// stations and its neighbours. The limits are set so that all the stations are
// too close (too far) from the reference event, but not from some of its
// neighbours (those 1 km further away from (closer to) the station).
BOOST_DATA_TEST_CASE(test_clustering_station_distance,
                     bdata::xrange(orgList.size()),
                     orgIdx)
{
  const Origin &org = orgList.at(orgIdx);
  const HDD::Catalog cat =
      buildCatalog(org.lat, org.lon, org.depth, 50, {1}, {"P"}, {});
  const EventTree tree = createEventTree(cat);
  const Event &event = cat.getEvents().at(1);

  double minDist = std::numeric_limits<double>::max(), maxDist = 0;
  for (const auto &kv : cat.getStations())
  {
    const double distance = computeDistance(event, kv.second);
    minDist               = std::min(minDist, distance);
    maxDist               = std::max(maxDist, distance);
  }

  // all stations too close to the reference event
  BOOST_CHECK_THROW(
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    maxDist + 0.1, //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2)  //  double maxNeighbourDist
      ,
      HDD::Exception);

  // all stations too far from the reference event
  BOOST_CHECK_THROW(
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    0,             //  double minESdis
                                    minDist - 0.1, //  double maxESdis
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2)  //  double maxNeighbourDist
      ,
      HDD::Exception);

  // all stations within the limits for the reference event and its neighbours
  const HDD::Neighbours neighbours =
      selectNeighbouringEventsThrow(tree, cat, event, cat,
                                    minDist - 0.9, //  double minESdis
                                    maxDist + 0.9, //  double maxESdis
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2); //  double maxNeighbourDist
  const unordered_set<unsigned> neighbourIds{2, 3, 4, 5};
  BOOST_CHECK(neighbours.ids() == neighbourIds);
  checkNeighbours(neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."}, {"P"});
}

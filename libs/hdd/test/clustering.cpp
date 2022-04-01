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

void addStationsToCatalog(HDD::Catalog &cat,
                          double lat,
                          double lon,
                          double distance)
{
  Station sta;
  double staLat, staLon;

  computeCoordinates(distance, 0, lat, lon, staLat, staLon);
  sta = {"NET.ST01", staLat, staLon, 250, "NET", "ST01", ""};
  cat.addStation(sta);
  computeCoordinates(distance, 90, lat, lon, staLat, staLon);
  sta = {"NET.ST02", staLat, staLon, 295, "NET", "ST02", ""};
  cat.addStation(sta);
  computeCoordinates(distance, 180, lat, lon, staLat, staLon);
  sta = {"NET.ST03", staLat, staLon, 301, "NET", "ST03", ""};
  cat.addStation(sta);
  computeCoordinates(distance, 270, lat, lon, staLat, staLon);
  sta = {"NET.ST04", staLat, staLon, 395, "NET", "ST04", ""};
  cat.addStation(sta);
}

void addEventToCatalog(HDD::Catalog &cat, double lat, double lon, double depth)
{
  Event ev{0};
  ev.latitude            = lat;
  ev.longitude           = lon;
  ev.depth               = depth;
  const unsigned eventId = cat.addEvent(ev);

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
    ph.isManual         = true;
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
    computeCoordinates(distance, 45, event.latitude, event.longitude,
                       neighbourLat, neighbourLon);
    addEventToCatalog(cat, neighbourLat, neighbourLon, depth);
    computeCoordinates(distance, 135, event.latitude, event.longitude,
                       neighbourLat, neighbourLon);
    addEventToCatalog(cat, neighbourLat, neighbourLon, depth);
    computeCoordinates(distance, 225, event.latitude, event.longitude,
                       neighbourLat, neighbourLon);
    addEventToCatalog(cat, neighbourLat, neighbourLon, depth);
    computeCoordinates(distance, 315, event.latitude, event.longitude,
                       neighbourLat, neighbourLon);
    addEventToCatalog(cat, neighbourLat, neighbourLon, depth);
  }
}

HDD::Catalog
buildCatalog(double lat,
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
  return HDD::Catalog::filterPhasesAndSetWeights(cat, Phase::Source::CATALOG,
                                                 PphaseToKeep, SphaseToKeep);
}

void checkNeighbours(const HDD::Neighbours &neighbours,
                     const unordered_set<unsigned> &mustBePresentIds,
                     const unordered_set<string> &mustBePresentStations,
                     const unordered_set<Phase::Type> &mustBePresentPhases)
{
  for (const unsigned id : mustBePresentIds)
    for (const string &sta : mustBePresentStations)
      for (const Phase::Type pht : mustBePresentPhases)
        BOOST_CHECK(neighbours.has(id, sta, pht));
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

BOOST_DATA_TEST_CASE(test_clustering1, bdata::xrange(orgList.size()), orgIdx)
{
  const Origin &org = orgList.at(orgIdx);
  const HDD::Catalog cat =
      buildCatalog(org.lat, org.lon, org.depth, 50, {1}, {"P"}, {});

  // cat->writeToFile(strf("test_clustering1-cat%d-event.csv",orgIdx),
  //                 strf("test_clustering1-cat%d-phase.csv",orgIdx),
  //                 strf("test_clustering1-cat%d-station.csv",orgIdx));

  const Event &event = cat.getEvents().at(1);
  const unordered_set<unsigned> neighbourIds{2, 3, 4, 5};

  unique_ptr<HDD::Neighbours> neighbours;

  // neighbours out of rangei: maxEllipsoidSize
  BOOST_CHECK_THROW(
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh   0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    0.5,   //  double maxEllipsoidSize
                                    false) //  bool keepUnmatched
      ,
      HDD::Exception);

  // stations out of range: minESdis
  BOOST_CHECK_THROW(
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    52, //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false) //  bool keepUnmatched
      ,
      HDD::Exception);

  // stations out of range: maxESdis
  BOOST_CHECK_THROW(
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    48, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false) //  bool keepUnmatched
      ,
      HDD::Exception);

  // minEStoIEratio not enough
  BOOST_CHECK_THROW(
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    60, //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false) //  bool keepUnmatched
      ,
      HDD::Exception);

  // too many DTs per event requested
  BOOST_CHECK_THROW(
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    5,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false) //  bool keepUnmatched
      ,
      HDD::Exception);

  // too many neighbours requested
  BOOST_CHECK_THROW(
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    5,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false) //  bool keepUnmatched
      ,
      HDD::Exception);

  // maxEllipsoidSize/numEllipsoids
  neighbours =
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    0,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false); //  bool keepUnmatched

  BOOST_CHECK_EQUAL(neighbours->refEvId, event.id);
  BOOST_CHECK(neighbours->ids == neighbourIds);
  checkNeighbours(*neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                  {Phase::Type::P});

  // maxEllipsoidSize/numEllipsoids
  neighbours =
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false); //  bool keepUnmatched

  BOOST_CHECK_EQUAL(neighbours->refEvId, event.id);
  BOOST_CHECK(neighbours->ids == neighbourIds);
  checkNeighbours(*neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                  {Phase::Type::P});

  // maxEllipsoidSize/numEllipsoids
  neighbours =
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    5,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false); //  bool keepUnmatched

  BOOST_CHECK_EQUAL(neighbours->refEvId, event.id);
  BOOST_CHECK(neighbours->ids == neighbourIds);
  checkNeighbours(*neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                  {Phase::Type::P});

  // minEStoIEratio
  neighbours =
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    48, //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh   0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false); //  bool keepUnmatched

  BOOST_CHECK_EQUAL(neighbours->refEvId, event.id);
  BOOST_CHECK(neighbours->ids == neighbourIds);
  checkNeighbours(*neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                  {Phase::Type::P});

  // maxNumNeigh
  neighbours =
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    2,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false); //  bool keepUnmatched

  BOOST_CHECK_EQUAL(neighbours->refEvId, event.id);
  BOOST_CHECK(neighbours->ids.size() == 2);

  // maxNumNeigh
  neighbours =
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    3,  //  unsigned maxNumNeigh   0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false); //  bool keepUnmatched

  BOOST_CHECK_EQUAL(neighbours->refEvId, event.id);
  BOOST_CHECK(neighbours->ids.size() == 3);

  // minDTperEvt
  neighbours =
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    4,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    5,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false); //  bool keepUnmatched

  BOOST_CHECK_EQUAL(neighbours->refEvId, event.id);
  BOOST_CHECK(neighbours->ids == neighbourIds);
  checkNeighbours(*neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                  {Phase::Type::P});

  // maxDTperEvt
  neighbours =
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    1,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false); //  bool keepUnmatched

  BOOST_CHECK_EQUAL(neighbours->refEvId, event.id);
  BOOST_CHECK(neighbours->ids == neighbourIds);
  checkNeighbours(*neighbours, neighbourIds, {"NET.ST01."}, {Phase::Type::P});

  // maxDTperEvt
  neighbours =
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    2,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    0,  //  unsigned maxNumNeigh  0 = no limits
                                    1,  //  unsigned numEllipsoids
                                    2,  //  double maxEllipsoidSize
                                    false); //  bool keepUnmatched

  BOOST_CHECK_EQUAL(neighbours->refEvId, event.id);
  BOOST_CHECK(neighbours->ids == neighbourIds);
  checkNeighbours(*neighbours, neighbourIds, {"NET.ST01.", "NET.ST02."},
                  {Phase::Type::P});
}

BOOST_DATA_TEST_CASE(test_clustering2, bdata::xrange(orgList.size()), orgIdx)
{
  const Origin &org = orgList.at(orgIdx);
  const HDD::Catalog cat =
      buildCatalog(org.lat, org.lon, org.depth, 70,
                   {0.5, 0.5, 0.6, 1.5, 1.5, 1.6, 3.0, 3.0, 3.1}, {}, {"S"});

  // cat->writeToFile(strf("test_clustering2-cat%d-event.csv",orgIdx),
  //                 strf("test_clustering2-cat%d-phase.csv",orgIdx),
  //                 strf("test_clustering2-cat%d-station.csv",orgIdx));

  const Event &event = cat.getEvents().at(1);
  unordered_set<unsigned> neighbourIds;

  // test multiple ellipsoids
  unique_ptr<HDD::Neighbours> neighbours;
  neighbours =
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    0,  //  double minESdis
                                    -1, //  double maxESdis      -1 = no limits
                                    0,  //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    1,  //  unsigned minNumNeigh
                                    24, //  unsigned maxNumNeigh  0 = no limits
                                    3,  //  unsigned numEllipsoids
                                    4,  //  double maxEllipsoidSize
                                    false); //  bool keepUnmatched

  neighbourIds = {2,  3,  4,  5,  6,  7,  8,  9,  14, 15, 16, 17,
                  18, 19, 20, 21, 26, 27, 28, 29, 30, 31, 32, 33};

  BOOST_CHECK_EQUAL(neighbours->refEvId, event.id);
  BOOST_CHECK(neighbours->ids == neighbourIds);
  checkNeighbours(*neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                  {Phase::Type::S});

  // test multiple ellipsoids
  neighbours =
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    50, //  double minESdis
                                    91, //  double maxESdis      -1 = no limits
                                    10, //  double minEStoIEratio
                                    4,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    16, //  unsigned minNumNeigh
                                    16, //  unsigned maxNumNeigh  0 = no limits
                                    3,  //  unsigned numEllipsoids
                                    4,  //  double maxEllipsoidSize
                                    false); //  bool keepUnmatched

  neighbourIds = {2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21};

  BOOST_CHECK_EQUAL(neighbours->refEvId, event.id);
  BOOST_CHECK(neighbours->ids == neighbourIds);
  checkNeighbours(*neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                  {Phase::Type::S});

  // test multiple ellipsoids
  neighbours =
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    50, //  double minESdis
                                    91, //  double maxESdis      -1 = no limits
                                    10, //  double minEStoIEratio
                                    1,  //  unsigned minDTperEvt
                                    2,  //  unsigned maxDTperEvt  0 = no limits
                                    16, //  unsigned minNumNeigh
                                    16, //  unsigned maxNumNeigh  0 = no limits
                                    3,  //  unsigned numEllipsoids
                                    4,  //  double maxEllipsoidSize
                                    false); //  bool keepUnmatched

  neighbourIds = {2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21};

  BOOST_CHECK_EQUAL(neighbours->refEvId, event.id);
  BOOST_CHECK(neighbours->ids == neighbourIds);
  checkNeighbours(*neighbours, neighbourIds, {"NET.ST01.", "NET.ST02."},
                  {Phase::Type::S});

  // test multiple ellipsoids
  neighbours =
      HDD::selectNeighbouringEvents(cat, event, cat,
                                    0,  //  double minPhaseWeight
                                    50, //  double minESdis
                                    91, //  double maxESdis      -1 = no limits
                                    10, //  double minEStoIEratio
                                    4,  //  unsigned minDTperEvt
                                    0,  //  unsigned maxDTperEvt  0 = no limits
                                    8,  //  unsigned minNumNeigh
                                    8,  //  unsigned maxNumNeigh  0 = no limits
                                    3,  //  unsigned numEllipsoids
                                    4,  //  double maxEllipsoidSize
                                    false); //  bool keepUnmatched

  neighbourIds = {2, 3, 4, 5, 6, 7, 8, 9};

  BOOST_CHECK_EQUAL(neighbours->refEvId, event.id);
  BOOST_CHECK(neighbours->ids == neighbourIds);
  checkNeighbours(*neighbours, neighbourIds,
                  {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                  {Phase::Type::S});
}

BOOST_DATA_TEST_CASE(test_clustering3, bdata::xrange(orgList.size()), orgIdx)
{
  const Origin &org = orgList.at(orgIdx);
  const HDD::Catalog cat =
      buildCatalog(org.lat, org.lon, org.depth, 100, {7.5, 15}, {"P"}, {"S"});

  // cat->writeToFile(strf("test_clustering3-cat%d-event.csv",orgIdx),
  //                 strf("test_clustering3-cat%d-phase.csv",orgIdx),
  //                 strf("test_clustering3-cat%d-station.csv",orgIdx));

  // find Neighbours for each event in the catalog
  unordered_map<unsigned, unique_ptr<HDD::Neighbours>> neighboursByEvent =
      HDD::selectNeighbouringEventsCatalog(
          cat,
          0,      //  double minPhaseWeight
          0,      //  double minESdis
          -1,     //  double maxESdis      -1 = no limits
          0,      //  double minEStoIEratio
          1,      //  unsigned minDTperEvt
          0,      //  unsigned maxDTperEvt  0 = no limits
          1,      //  unsigned minNumNeigh
          0,      //  unsigned maxNumNeigh   0 = no limits
          0,      //  unsigned numEllipsoids
          8,      //  double maxEllipsoidSize
          false); //  bool keepUnmatched

  BOOST_REQUIRE(neighboursByEvent.size() == 9);

  const unordered_map<unsigned, unordered_set<unsigned>> expectedNeighbours = {
      {1, {2, 3, 4, 5}}, {2, {1, 6}}, {3, {1, 7}}, {4, {1, 8}}, {5, {1, 9}},
      {6, {2}},          {7, {3}},    {8, {4}},    {9, {5}}};

  for (const auto &kv : neighboursByEvent)
  {
    const unique_ptr<HDD::Neighbours> &neighbours = kv.second;
    BOOST_REQUIRE(neighbours->refEvId == kv.first);
    BOOST_REQUIRE_NO_THROW(expectedNeighbours.at(neighbours->refEvId));
    const unordered_set<unsigned> &neighbourIds =
        expectedNeighbours.at(neighbours->refEvId);
    BOOST_REQUIRE(neighbours->ids == neighbourIds);
    checkNeighbours(*neighbours, neighbourIds,
                    {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                    {Phase::Type::P, Phase::Type::S});
  }

  list<unordered_map<unsigned, unique_ptr<HDD::Neighbours>>> clusters =
      HDD::clusterizeNeighbouringEvents(neighboursByEvent);

  BOOST_REQUIRE(clusters.size() == 1);

  for (const auto &cluster : clusters)
    for (const auto &kv : cluster)
    {
      const unique_ptr<HDD::Neighbours> &neighbours = kv.second;
      BOOST_REQUIRE(neighbours->refEvId == kv.first);
      BOOST_REQUIRE_NO_THROW(expectedNeighbours.at(neighbours->refEvId));
      const unordered_set<unsigned> &neighbourIds =
          expectedNeighbours.at(neighbours->refEvId);
      for (const unsigned neighId : neighbourIds)
      {
        if (neighbours->ids.count(neighId) == 1) // is in our Neighbours
        {
          checkNeighbours(*neighbours, {neighId},
                          {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                          {Phase::Type::P, Phase::Type::S});
        }
        else // must be in the other event's Neighbours
        {
          const unique_ptr<HDD::Neighbours> &otherNeighbours =
              cluster.at(neighId);
          checkNeighbours(*otherNeighbours, {neighbours->refEvId},
                          {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                          {Phase::Type::P, Phase::Type::S});
        }
      }
    }
}

BOOST_DATA_TEST_CASE(test_clustering4, bdata::xrange(orgList.size()), orgIdx)
{
  const Origin &org = orgList.at(orgIdx);
  const HDD::Catalog cat =
      buildCatalog(org.lat, org.lon, org.depth, 10, {3, 4}, {}, {"S"});

  // cat->writeToFile(strf("test_clustering4-cat%d-event.csv",orgIdx),
  //                 strf("test_clustering4-cat%d-phase.csv",orgIdx),
  //                 strf("test_clustering4-cat%d-station.csv",orgIdx));

  // find Neighbours for each event in the catalog
  unordered_map<unsigned, unique_ptr<HDD::Neighbours>> neighboursByEvent =
      HDD::selectNeighbouringEventsCatalog(
          cat,
          0,      //  double minPhaseWeight
          0,      //  double minESdis
          -1,     //  double maxESdis      -1 = no limits
          0,      //  double minEStoIEratio
          1,      //  unsigned minDTperEvt
          0,      //  unsigned maxDTperEvt  0 = no limits
          1,      //  unsigned minNumNeigh
          0,      //  unsigned maxNumNeigh   0 = no limits
          0,      //  unsigned numEllipsoids
          1.5,    //  double maxEllipsoidSize
          false); //  bool keepUnmatched

  BOOST_REQUIRE(neighboursByEvent.size() == 8);

  const unordered_map<unsigned, unordered_set<unsigned>> expectedNeighbours = {
      {2, {6}}, {3, {7}}, {4, {8}}, {5, {9}},
      {6, {2}}, {7, {3}}, {8, {4}}, {9, {5}}};

  for (const auto &kv : neighboursByEvent)
  {
    const unique_ptr<HDD::Neighbours> &neighbours = kv.second;
    BOOST_REQUIRE(neighbours->refEvId == kv.first);
    BOOST_REQUIRE_NO_THROW(expectedNeighbours.at(neighbours->refEvId));
    const unordered_set<unsigned> &neighbourIds =
        expectedNeighbours.at(neighbours->refEvId);
    BOOST_REQUIRE(neighbours->ids == neighbourIds);
    checkNeighbours(*neighbours, neighbourIds,
                    {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                    {Phase::Type::S});
  }

  list<unordered_map<unsigned, unique_ptr<HDD::Neighbours>>> clusters =
      HDD::clusterizeNeighbouringEvents(neighboursByEvent);

  BOOST_REQUIRE(clusters.size() == 4);

  for (const auto &cluster : clusters)
    for (const auto &kv : cluster)
    {
      const unique_ptr<HDD::Neighbours> &neighbours = kv.second;
      BOOST_REQUIRE(neighbours->refEvId == kv.first);
      BOOST_REQUIRE_NO_THROW(expectedNeighbours.at(neighbours->refEvId));
      const unordered_set<unsigned> &neighbourIds =
          expectedNeighbours.at(neighbours->refEvId);
      for (const unsigned neighId : neighbourIds)
      {
        if (neighbours->ids.count(neighId) == 1) // is in our Neighbours
        {
          checkNeighbours(*neighbours, {neighId},
                          {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                          {Phase::Type::S});
        }
        else // must be in the other event's Neighbours
        {
          const unique_ptr<HDD::Neighbours> &otherNeighbours =
              cluster.at(neighId);
          checkNeighbours(*otherNeighbours, {neighbours->refEvId},
                          {"NET.ST01.", "NET.ST02.", "NET.ST03.", "NET.ST04."},
                          {Phase::Type::S});
        }
      }
    }
}

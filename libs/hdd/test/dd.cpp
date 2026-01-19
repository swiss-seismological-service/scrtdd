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

#define BOOST_TEST_MODULE libhdd
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>

#include "common.ipp"
#include "hdd/catalog.h"
#include "hdd/dd.h"
#include "hdd/random.h"
#include "hdd/ttt.h"

using namespace std;
using namespace HDD;
using Event     = HDD::Catalog::Event;
using Phase     = HDD::Catalog::Phase;
using Station   = HDD::Catalog::Station;
namespace bdata = boost::unit_test::data;

namespace {

void addStationsToCatalog(HDD::Catalog &cat, double lat, double lon)
{
  Station sta;
  double staLat, staLon;

  double distance = 25;

  computeCoordinates(distance, degToRad(0), lat, lon, staLat, staLon);
  sta = {"NT.ST01", staLat, staLon,  500, "NT", "STA01", ""};
  cat.addStation(sta);
  computeCoordinates(distance, degToRad(90), lat, lon, staLat, staLon);
  sta = {"NT.ST02", staLat, staLon, 1000, "NT", "STA02", ""};
  cat.addStation(sta);
  computeCoordinates(distance, degToRad(180), lat, lon, staLat, staLon);
  sta = {"NT.ST03", staLat, staLon, 1500, "NT", "STA03", ""};
  cat.addStation(sta);
  computeCoordinates(distance, degToRad(270), lat, lon, staLat, staLon);
  sta = {"NT.ST04", staLat, staLon, 2000, "NT", "STA04", ""};
  cat.addStation(sta);

  distance = 15;

  computeCoordinates(distance, degToRad(45), lat, lon, staLat, staLon);
  sta = {"NT.ST05", staLat, staLon, -2000, "NT", "STA05", ""};
  cat.addStation(sta);
  computeCoordinates(distance, degToRad(135), lat, lon, staLat, staLon);
  sta = {"NT.ST06", staLat, staLon, -1500, "NT", "STA06", ""};
  cat.addStation(sta);
  computeCoordinates(distance, degToRad(225), lat, lon, staLat, staLon);
  sta = {"NT.ST07", staLat, staLon, -1000, "NT", "STA07", ""};
  cat.addStation(sta);
  computeCoordinates(distance, degToRad(315), lat, lon, staLat, staLon);
  sta = {"NT.ST08", staLat, staLon,  -500, "NT", "STA08", ""};
  cat.addStation(sta);
}

void addNLLStationsToCatalog(HDD::Catalog &cat)
{
  // with NLL we have no freedom to put the stations where we want
  // because the grid files have to be generated with stations locations
  for (const auto &sta : ::stationList)
  {
    cat.addStation(sta);
  }
}

void addEventToCatalog(HDD::Catalog &cat,
                       HDD::TravelTimeTable &ttt,
                       const UTCTime &time,
                       double lat,
                       double lon,
                       double depth)
{
  Event ev{};
  ev.time                = time;
  ev.latitude            = lat;
  ev.longitude           = normalizeLon(lon);
  ev.depth               = depth;
  ev.magnitude           = 1.0;
  const unsigned eventId = cat.addEvent(ev, false);

  for (const auto &kv : cat.getStations())
  {
    const Station &sta = kv.second;

    Phase ph;
    ph.eventId          = eventId;
    ph.stationId        = sta.id;
    ph.lowerUncertainty = 0;
    ph.upperUncertainty = 0;
    ph.networkCode      = sta.networkCode;
    ph.stationCode      = sta.stationCode;
    ph.locationCode     = sta.locationCode;
    ph.channelCode      = "";

    double travelTime = ttt.compute(ev, sta, "P");
    ph.time           = ev.time + secToDur(travelTime);
    ph.type           = "P";
    cat.addPhase(ph);

    travelTime = ttt.compute(ev, sta, "S");
    ph.time    = ev.time + secToDur(travelTime);
    ph.type    = "S";
    cat.addPhase(ph);
  }
}

void addEvents1ToCatalog(HDD::Catalog &cat,
                         HDD::TravelTimeTable &ttt,
                         const UTCTime &time,
                         double lat,
                         double lon,
                         double depth,
                         int numEvents,
                         double extent)
{
  double distance = extent / 2;
  double startLon, endLon, dummy;
  computeCoordinates(distance, degToRad(-90), lat, lon, dummy, startLon, depth);
  computeCoordinates(distance, degToRad(90), lat, lon, dummy, endLon, depth);
  const double lonStep = (endLon - startLon) / (numEvents - 1);
  for (int evn = -(numEvents / 2); evn < std::ceil(numEvents / 2.0); evn++)
  {
    addEventToCatalog(cat, ttt, time + secToDur(evn * 61), lat,
                      lon + evn * lonStep, depth);
  }
}

void addEvents2ToCatalog(HDD::Catalog &cat,
                         HDD::TravelTimeTable &ttt,
                         const UTCTime &time,
                         double lat,
                         double lon,
                         double depth,
                         int numEvents,
                         double extent)
{
  double distance = extent / 2;
  double startLat, endLat, dummy;
  computeCoordinates(distance, degToRad(180), lat, lon, startLat, dummy, depth);
  computeCoordinates(distance, degToRad(0), lat, lon, endLat, dummy, depth);
  double latStep = (endLat - startLat) / (numEvents - 1);
  for (int evn = -(numEvents / 2); evn < std::ceil(numEvents / 2.0); evn++)
  {
    addEventToCatalog(cat, ttt, time + secToDur(evn * 122), lat + evn * latStep,
                      lon, depth);
  }
}

void addEvents3ToCatalog(HDD::Catalog &cat,
                         HDD::TravelTimeTable &ttt,
                         const UTCTime &time,
                         double lat,
                         double lon,
                         double depth,
                         int numEvents,
                         double extent)
{
  double depthStep = extent / (numEvents - 1);
  for (int evn = -(numEvents / 2); evn < std::ceil(numEvents / 2.0); evn++)
  {
    addEventToCatalog(cat, ttt, time + secToDur(evn * 183), lat, lon,
                      depth + evn * depthStep);
  }
}

HDD::Catalog buildCatalog(HDD::TravelTimeTable &ttt,
                          bool nllStations,
                          const UTCTime &time,
                          double lat,
                          double lon,
                          double depth,
                          int numEvents, // multiple of 6
                          double extent)
{
  HDD::Catalog cat;

  if (nllStations)
    addNLLStationsToCatalog(cat);
  else
    addStationsToCatalog(cat, lat, lon);

  double distance = extent * 2.5;
  double clusterLat, clusterLon;
  // cluster 1
  computeCoordinates(distance, degToRad(135), lat, lon, clusterLat, clusterLon,
                     depth);
  addEvents1ToCatalog(cat, ttt, time + secToDur(1), clusterLat, clusterLon,
                      depth, numEvents / 6, extent);
  // cluster 2
  computeCoordinates(distance, degToRad(315), lat, lon, clusterLat, clusterLon,
                     depth);
  addEvents2ToCatalog(cat, ttt, time + secToDur(2), clusterLat, clusterLon,
                      depth, numEvents / 6, extent);
  // cluster 3
  computeCoordinates(distance, degToRad(45), lat, lon, clusterLat, clusterLon,
                     depth);
  addEvents3ToCatalog(cat, ttt, time + secToDur(3), clusterLat, clusterLon,
                      depth, numEvents / 6, extent);
  // cluster 4
  computeCoordinates(distance, degToRad(225), lat, lon, clusterLat, clusterLon,
                     depth);
  addEvents1ToCatalog(cat, ttt, time + secToDur(4), clusterLat, clusterLon,
                      depth, numEvents / 6, extent);
  addEvents2ToCatalog(cat, ttt, time + secToDur(4), clusterLat, clusterLon,
                      depth, numEvents / 6, extent);
  addEvents3ToCatalog(cat, ttt, time + secToDur(4), clusterLat, clusterLon,
                      depth, numEvents / 6, extent);

  return Catalog::fillProcessingInfo(cat, Phase::Source::CATALOG, {"P"},
                                            {"S"}, {});
}

HDD::Catalog buildBackgroundCatalog(HDD::TravelTimeTable &ttt,
                                    bool nllStations,
                                    const UTCTime &time,
                                    double lat,
                                    double lon,
                                    double depth,
                                    int numEvents, // multiple of 3
                                    double extent)
{
  HDD::Catalog cat;

  if (nllStations)
    addNLLStationsToCatalog(cat);
  else
    addStationsToCatalog(cat, lat, lon);

  addEvents1ToCatalog(cat, ttt, time + secToDur(1), lat, lon, depth,
                      numEvents / 3, extent);
  addEvents2ToCatalog(cat, ttt, time + secToDur(2), lat, lon, depth,
                      numEvents / 3, extent);
  addEvents3ToCatalog(cat, ttt, time + secToDur(3), lat, lon, depth,
                      numEvents / 3, extent);

  return Catalog::fillProcessingInfo(cat, Phase::Source::CATALOG, {"P"},
                                            {"S"}, {});
}

void randomPerturbation(HDD::Catalog &cat)
{
  // random changes to all events (mean of all changes is != 0)
  HDD::NormalRandomer timeDist(0.1, 0.2, 0x1004);  // sec
  HDD::NormalRandomer latDist(1.5, 0.3, 0x1001);   // km
  HDD::NormalRandomer lonDist(-0.7, 0.8, 0x1002);  // km
  HDD::NormalRandomer depthDist(1.5, 0.7, 0x1003); // km
  for (const auto &kv : cat.getEvents())
  {
    Event ev = kv.second;
    ev.time += secToDur(timeDist.next());

    ev.depth += depthDist.next();
    if (ev.depth < 0) ev.depth = 0;

    double distance = latDist.next();
    if (distance >= 0)
      computeCoordinates(distance, degToRad(180), ev.latitude, ev.longitude,
                         ev.latitude, ev.longitude, ev.depth);
    else
      computeCoordinates(-distance, degToRad(0), ev.latitude, ev.longitude,
                         ev.latitude, ev.longitude, ev.depth);

    distance = lonDist.next();
    if (distance >= 0)
      computeCoordinates(distance, degToRad(90), ev.latitude, ev.longitude,
                         ev.latitude, ev.longitude, ev.depth);
    else
      computeCoordinates(-distance, degToRad(270), ev.latitude, ev.longitude,
                         ev.latitude, ev.longitude, ev.depth);

    cat.updateEvent(ev);
  }
}

HDD::Catalog relocateCatalog(const HDD::Catalog &cat,
                             unique_ptr<HDD::TravelTimeTable> ttt,
                             const string &workingDir)
{
  HDD::Config ddCfg;

  HDD::DD dd(cat, ddCfg, std::move(ttt));

  bool saveProcessing = false; // for debugging

  HDD::ClusteringOptions clusterOpt;
  clusterOpt.numEllipsoids    = 0;
  clusterOpt.maxNeighbourDist = 15;
  clusterOpt.minEvStaToInterEvRatio   = 0;

  HDD::XcorrOptions xcorrOpt;
  xcorrOpt.enable = false; // disable cross-correlation

  HDD::SolverOptions solverOpt;
  solverOpt.algoIterations               = 20;
  solverOpt.absLocConstraintStart        = 0.3;
  solverOpt.absLocConstraintEnd          = 0.3;
  solverOpt.dampingFactorStart           = 0.01;
  solverOpt.dampingFactorEnd             = 0.01;
  solverOpt.downWeightingByResidualStart = 0;
  solverOpt.downWeightingByResidualEnd   = 0;

  std::list<std::unordered_map<unsigned, Neighbours>> clusters;
  XCorrCache xcorrData;
  return dd.relocateMultiEvents(clusters, xcorrData, clusterOpt, xcorrOpt,
                                solverOpt, saveProcessing);
}

HDD::Catalog relocateSingleEvent(const HDD::Catalog &bgCat,
                                 unique_ptr<HDD::TravelTimeTable> ttt,
                                 const string &workingDir,
                                 const HDD::Catalog &realTimeCat)
{
  HDD::Config ddCfg;

  HDD::DD dd(bgCat, ddCfg, std::move(ttt));

  bool saveProcessing = false; // for debugging

  HDD::ClusteringOptions clusterOpt;
  clusterOpt.numEllipsoids    = 5;
  clusterOpt.maxNeighbourDist = 15;
  clusterOpt.maxNumNeigh      = 40;
  clusterOpt.minEvStaToInterEvRatio   = 0;

  HDD::XcorrOptions xcorrOpt;
  xcorrOpt.enable = false; // disable cross-correlation

  HDD::SolverOptions solverOpt;
  solverOpt.algoIterations        = 20;
  solverOpt.absLocConstraintStart = 0.3;
  solverOpt.absLocConstraintEnd   = 0.3;
  solverOpt.dampingFactorStart    = 0.01;
  solverOpt.dampingFactorEnd      = 0.01;

  HDD::Catalog relocCat;
  for (const auto &kv : realTimeCat.getEvents())
  {
    const Event &ev            = kv.second;
    HDD::Catalog orgToRelocate = realTimeCat.extractEvent(ev.id, false);
    HDD::Catalog relocatedEvent =
        dd.relocateSingleEvent(orgToRelocate, true, clusterOpt, clusterOpt,
                               xcorrOpt, solverOpt, saveProcessing);
    relocCat.add(relocatedEvent, false);
  }

  return relocCat;
}

void testCatalogEqual(const HDD::Catalog &cat1, const HDD::Catalog &cat2)
{
  BOOST_REQUIRE_EQUAL(cat1.getEvents().size(), cat2.getEvents().size());
  for (const auto &kv : cat1.getEvents())
  {
    const Event &ev1 = kv.second;
    BOOST_REQUIRE_EQUAL(cat2.getEvents().count(ev1.id), 1);
    const Event &ev2 = cat2.getEvents().at(ev1.id);
    BOOST_CHECK_SMALL(durToSec((ev1.time - ev2.time)),
                      0.01); // tolerance 10 millisec
    BOOST_CHECK_SMALL(computeDistance(ev1.latitude, ev1.longitude,
                                      ev2.latitude, ev2.longitude),
                      0.01); // tolerance 10 meter
    BOOST_CHECK_SMALL(ev1.depth-ev2.depth, 0.05); // tolerance 50 m
  }
}

struct Centroid
{
  double lat;
  double lon;
  double depth;
};

const vector<Centroid> centroidList{
  {nllGridCentroidLat, nllGridCentroidLon, 0},
  {nllGridCentroidLat, nllGridCentroidLon, -1.5},
  {nllGridCentroidLat, nllGridCentroidLon, 1.5},
  {nllGridCentroidLat, nllGridCentroidLon, 8},
  {0, 0, -1},
  { 85, -90, 0},
  {-85, 90, 3},
  { 75, 170, 8},
  { -60, -170, 2},
  { 0, 179.99, 4},
};

} // namespace

BOOST_DATA_TEST_CASE(test_dd_multi_event1,
                     bdata::xrange(centroidList.size()) *
                         bdata::xrange(tttList.size()),
                     cIdx,
                     tttIdx)
{
  const Centroid &centroid = centroidList.at(cIdx);
  bool nllStations         = (tttList.at(tttIdx).type == "Native:NonLinLoc");

  if (nllStations && 
   (centroid.lat != nllGridCentroidLat || centroid.lon != nllGridCentroidLon))
  {
    BOOST_TEST_MESSAGE(
        "Skipping NonLinLoc test in a location without grid files");
    BOOST_CHECK(true); // Needed to suppress "Test case [...] did not check any
                       // assertions"
    return;
  }

  unique_ptr<HDD::TravelTimeTable> ttt = createTTT(tttList.at(tttIdx));

  const HDD::Catalog baseCat =
      buildCatalog(*ttt, nllStations, UTCClock::fromDate(2001, 1, 2, 0, 0, 0, 0),
                   centroid.lat, centroid.lon, centroid.depth, 42, 3.0);

  // Test no event changes, the relocation should not move the events
  string workingDir =
      strf("./data/test_dd_multi_event1_%lu_%lu_reloc", cIdx, tttIdx);
  const HDD::Catalog relocCat =
      relocateCatalog(baseCat, std::move(ttt), workingDir);
  testCatalogEqual(baseCat, relocCat);
}

BOOST_DATA_TEST_CASE(test_dd_multi_event2,
                     bdata::xrange(centroidList.size()) *
                         bdata::xrange(tttList.size()),
                     cIdx,
                     tttIdx)
{
  const Centroid &centroid = centroidList.at(cIdx);
  bool nllStations         = (tttList.at(tttIdx).type == "Native:NonLinLoc");

  if (nllStations && 
   (centroid.lat != nllGridCentroidLat || centroid.lon != nllGridCentroidLon))
  {
    BOOST_TEST_MESSAGE(
        "Skipping NonLinLoc test in a location without grid files");
    BOOST_CHECK(true); // Needed to suppress "Test case [...] did not check any
                       // assertions"
    return;
  }

  unique_ptr<HDD::TravelTimeTable> ttt = createTTT(tttList.at(tttIdx));

  const HDD::Catalog baseCat =
      buildCatalog(*ttt, nllStations, UTCClock::fromDate(1914, 11, 23, 11, 34, 23, 1234),
                   centroid.lat, centroid.lon, centroid.depth, 42, 1.5);

  // Test 2: all events at centroid location, with random perturbation of time
  HDD::Catalog cat(baseCat);
  HDD::NormalRandomer timeDist(0.2, 0.400, 0x1004); // sec
  for (const auto &kv : cat.getEvents())
  {
    Event ev = kv.second;
    ev.time += secToDur(timeDist.next());
    ev.latitude  = centroid.lat;
    ev.longitude = centroid.lon;
    ev.depth     = centroid.depth;
    cat.updateEvent(ev);
  }

  string workingDir =
      strf("./data/test_dd_multi_event2_%lu_%lu_reloc", cIdx, tttIdx);
  const HDD::Catalog relocCat =
      relocateCatalog(cat, std::move(ttt), workingDir);
  testCatalogEqual(baseCat, relocCat);
}

BOOST_DATA_TEST_CASE(test_dd_multi_event3,
                     bdata::xrange(centroidList.size()) *
                         bdata::xrange(tttList.size()),
                     cIdx,
                     tttIdx)
{
  const Centroid &centroid = centroidList.at(cIdx);
  bool nllStations         = (tttList.at(tttIdx).type == "Native:NonLinLoc");

  if (nllStations && 
   (centroid.lat != nllGridCentroidLat || centroid.lon != nllGridCentroidLon))
  {
    BOOST_TEST_MESSAGE(
        "Skipping NonLinLoc test in a location without grid files");
    BOOST_CHECK(true); // Needed to suppress "Test case [...] did not check any
                       // assertions"
    return;
  }

  unique_ptr<HDD::TravelTimeTable> ttt = createTTT(tttList.at(tttIdx));

  const HDD::Catalog baseCat =
      buildCatalog(*ttt, nllStations, UTCClock::fromDate(2041, 6, 14, 23, 46, 3, 65),
                   centroid.lat, centroid.lon, centroid.depth, 42, 1.0);

  // Test 3: random changes to all events (mean of all changes is != 0)
  HDD::Catalog cat(baseCat);
  randomPerturbation(cat);

  string workingDir =
      strf("./data/test_dd_multi_event3_%lu_%lu_reloc", cIdx, tttIdx);
  const HDD::Catalog relocCat =
      relocateCatalog(cat, std::move(ttt), workingDir);
  testCatalogEqual(baseCat, relocCat);
}

BOOST_DATA_TEST_CASE(test_dd_single_event1,
                     bdata::xrange(centroidList.size()) *
                         bdata::xrange(tttList.size()),
                     cIdx,
                     tttIdx)
{
  const Centroid &centroid = centroidList.at(cIdx);
  bool nllStations         = (tttList.at(tttIdx).type == "Native:NonLinLoc");

  if (nllStations && 
   (centroid.lat != nllGridCentroidLat || centroid.lon != nllGridCentroidLon))
  {
    BOOST_TEST_MESSAGE(
        "Skipping NonLinLoc test in a location without grid files");
    BOOST_CHECK(true); // Needed to suppress "Test case [...] did not check any
                       // assertions"
    return;
  }

  unique_ptr<HDD::TravelTimeTable> ttt = createTTT(tttList.at(tttIdx));

  const UTCTime clusterTime =
      UTCClock::fromDate(1973, 12, 4, 14, 53, 32, 69854);
  const HDD::Catalog backgroundCat = buildBackgroundCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 21, 2.0);
  const HDD::Catalog realTimeCat = buildCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 24, 2.0);

  // Test 1: no event changes
  string workingDir =
      strf("./data/test_dd_single_event1_%lu_%lu_reloc", cIdx, tttIdx);
  HDD::Catalog relocCat = relocateSingleEvent(backgroundCat, std::move(ttt),
                                              workingDir, realTimeCat);
  testCatalogEqual(realTimeCat, relocCat);
}

BOOST_DATA_TEST_CASE(test_dd_single_event2,
                     bdata::xrange(centroidList.size()) *
                         bdata::xrange(tttList.size()),
                     cIdx,
                     tttIdx)
{
  const Centroid &centroid = centroidList.at(cIdx);
  bool nllStations         = (tttList.at(tttIdx).type == "Native:NonLinLoc");

  if (nllStations &&  
   (centroid.lat != nllGridCentroidLat || centroid.lon != nllGridCentroidLon))
  {
    BOOST_TEST_MESSAGE(
        "Skipping NonLinLoc test in a location without grid files");
    BOOST_CHECK(true); // Needed to suppress "Test case [...] did not check any
                       // assertions"
    return;
  }

  unique_ptr<HDD::TravelTimeTable> ttt = createTTT(tttList.at(tttIdx));

  const UTCTime clusterTime =
      UTCClock::fromDate(1989, 7, 18, 23, 59, 59, 999999);
  const HDD::Catalog backgroundCat = buildBackgroundCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 21, 1.5);
  const HDD::Catalog realTimeCat = buildCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 24, 1.5);

  // Test 2: all events at centroid location, with random perturbation of time
  HDD::Catalog cat(realTimeCat);
  HDD::NormalRandomer timeDist(0.3, 0.500, 0x1004); // sec
  for (const auto &kv : cat.getEvents())
  {
    Event ev = kv.second;
    ev.time += secToDur(timeDist.next());
    ev.latitude  = centroid.lat;
    ev.longitude = centroid.lon;
    ev.depth     = centroid.depth;
    cat.updateEvent(ev);
  }

  string workingDir =
      strf("./data/test_dd_single_event2_%lu_%lu_reloc", cIdx, tttIdx);
  const HDD::Catalog relocCat = relocateSingleEvent(
      backgroundCat, std::move(ttt), workingDir, realTimeCat);
  testCatalogEqual(realTimeCat, relocCat);
}

BOOST_DATA_TEST_CASE(test_dd_single_event3,
                     bdata::xrange(centroidList.size()) *
                         bdata::xrange(tttList.size()),
                     cIdx,
                     tttIdx)
{
  const Centroid &centroid = centroidList.at(cIdx);
  bool nllStations         = (tttList.at(tttIdx).type == "Native:NonLinLoc");

  if (nllStations &&
   (centroid.lat != nllGridCentroidLat || centroid.lon != nllGridCentroidLon))
  {
    BOOST_TEST_MESSAGE(
        "Skipping NonLinLoc test in a location without grid files");
    BOOST_CHECK(true); // Needed to suppress "Test case [...] did not check any
                       // assertions"
    return;
  }

  unique_ptr<HDD::TravelTimeTable> ttt = createTTT(tttList.at(tttIdx));

  const UTCTime clusterTime  = UTCClock::fromDate(2034, 4, 8, 3, 24, 54, 2354);
  const HDD::Catalog backgroundCat = buildBackgroundCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 21, 1.0);
  const HDD::Catalog realTimeCat = buildCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 24, 1.0);

  // Test 3: random changes to all events (mean of all changes is != 0)
  HDD::Catalog cat(realTimeCat);
  randomPerturbation(cat);

  string workingDir =
      strf("./data/test_dd_single_event3_%lu_%lu_reloc", cIdx, tttIdx);
  const HDD::Catalog relocCat = relocateSingleEvent(
      backgroundCat, std::move(ttt), workingDir, realTimeCat);
  testCatalogEqual(realTimeCat, relocCat);
}

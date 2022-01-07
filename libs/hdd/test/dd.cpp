#define SEISCOMP_TEST_MODULE hdd
#include <seiscomp/unittest/unittests.h>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>

#include "catalog.h"
#include "dd.h"
#include "ttt.h"

#include <seiscomp/logging/log.h>
#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/math.h>
#include <vector>

#include "common.ipp"

using namespace std;
using namespace Seiscomp;
using HDD::strf;
using Event     = HDD::Catalog::Event;
using Phase     = HDD::Catalog::Phase;
using Station   = HDD::Catalog::Station;
namespace bdata = boost::unit_test::data;

namespace {

void addStationsToCatalog(HDD::Catalog& cat, int maxStations)
{
  int numStations = 0;
  for (const auto &sta : stationList)
  {
    if (numStations++ == maxStations) break;
    cat.addStation(sta);
  }
}

void addEventToCatalog(HDD::Catalog &cat,
                       HDD::TravelTimeTable &ttt,
                       Core::Time time,
                       double lat,
                       double lon,
                       double depth)
{
  Event ev{0};
  ev.time                = time;
  ev.latitude            = lat;
  ev.longitude           = lon;
  ev.depth               = depth;
  ev.magnitude           = 1.0;
  const unsigned eventId = cat.addEvent(ev);

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
    ph.isManual         = true;

    double travelTime;
    ttt.compute(ev, sta, "P", travelTime);
    ph.time = ev.time + Core::TimeSpan(travelTime);
    ph.type = "P";
    cat.addPhase(ph);

    ttt.compute(ev, sta, "S", travelTime);
    ph.time = ev.time + Core::TimeSpan(travelTime);
    ph.type = "S";
    cat.addPhase(ph);
  }
}

void addEvents1ToCatalog(HDD::Catalog &cat,
                         HDD::TravelTimeTable &ttt,
                         Core::Time time,
                         double lat,
                         double lon,
                         double depth,
                         int numEvents,
                         double extent)
{
  double distance = Math::Geo::km2deg(extent / 2);
  double startLon, endLon, dummy;
  Math::Geo::delandaz2coord(distance, -90, lat, lon, &dummy, &startLon);
  Math::Geo::delandaz2coord(distance, 90, lat, lon, &dummy, &endLon);
  const double lonStep = (endLon - startLon) / (numEvents - 1);
  for (int evn = -(numEvents / 2); evn < std::ceil(numEvents / 2.0); evn++)
  {
    addEventToCatalog(cat, ttt, time + Core::TimeSpan(evn * 61), lat,
                      lon + evn * lonStep, depth);
  }
}

void addEvents2ToCatalog(HDD::Catalog &cat,
                         HDD::TravelTimeTable &ttt,
                         Core::Time time,
                         double lat,
                         double lon,
                         double depth,
                         int numEvents,
                         double extent)
{
  double distance = Math::Geo::km2deg(extent / 2);
  double startLat, endLat, dummy;
  Math::Geo::delandaz2coord(distance, 180, lat, lon, &startLat, &dummy);
  Math::Geo::delandaz2coord(distance, 0, lat, lon, &endLat, &dummy);
  double latStep = (endLat - startLat) / (numEvents - 1);
  for (int evn = -(numEvents / 2); evn < std::ceil(numEvents / 2.0); evn++)
  {
    addEventToCatalog(cat, ttt, time + Core::TimeSpan(evn * 122),
                      lat + evn * latStep, lon, depth);
  }
}

void addEvents3ToCatalog(HDD::Catalog &cat,
                         HDD::TravelTimeTable &ttt,
                         Core::Time time,
                         double lat,
                         double lon,
                         double depth,
                         int numEvents,
                         double extent)
{
  double depthStep = extent / (numEvents - 1);
  for (int evn = -(numEvents / 2); evn < std::ceil(numEvents / 2.0); evn++)
  {
    addEventToCatalog(cat, ttt, time + Core::TimeSpan(evn * 183), lat, lon,
                      depth + evn * depthStep);
  }
}

std::unique_ptr<HDD::Catalog> buildCatalog(HDD::TravelTimeTable &ttt,
                             int maxStations,
                             Core::Time time,
                             double lat,
                             double lon,
                             double depth,
                             int numEvents,
                             double extent)
{
  std::unique_ptr<HDD::Catalog>cat(new HDD::Catalog());
  addStationsToCatalog(*cat, maxStations);
  double distance = Math::Geo::km2deg(extent * 2.5);
  double clusterLat, clusterLon;
  // cluster 1
  Math::Geo::delandaz2coord(distance, 135, lat, lon, &clusterLat, &clusterLon);
  addEvents1ToCatalog(*cat, ttt, time + Core::TimeSpan(1), clusterLat,
                      clusterLon, depth, numEvents / 6, extent);
  // cluster 2
  Math::Geo::delandaz2coord(distance, 315, lat, lon, &clusterLat, &clusterLon);
  addEvents2ToCatalog(*cat, ttt, time + Core::TimeSpan(2), clusterLat,
                      clusterLon, depth, numEvents / 6, extent);
  // cluster 3
  Math::Geo::delandaz2coord(distance, 45, lat, lon, &clusterLat, &clusterLon);
  addEvents3ToCatalog(*cat, ttt, time + Core::TimeSpan(3), clusterLat,
                      clusterLon, depth, numEvents / 6, extent);
  // cluster 4
  Math::Geo::delandaz2coord(distance, 225, lat, lon, &clusterLat, &clusterLon);
  addEvents1ToCatalog(*cat, ttt, time + Core::TimeSpan(4), clusterLat,
                      clusterLon, depth, numEvents / 6, extent);
  addEvents2ToCatalog(*cat, ttt, time + Core::TimeSpan(4), clusterLat,
                      clusterLon, depth, numEvents / 6, extent);
  addEvents3ToCatalog(*cat, ttt, time + Core::TimeSpan(4), clusterLat,
                      clusterLon, depth, numEvents / 6, extent);
  return cat;
}

std::unique_ptr<HDD::Catalog> buildBackgroundCatalog(HDD::TravelTimeTable &ttt,
                                       int maxStations,
                                       Core::Time time,
                                       double lat,
                                       double lon,
                                       double depth,
                                       int numEvents,
                                       double extent)
{
  std::unique_ptr<HDD::Catalog>cat(new HDD::Catalog());
  addStationsToCatalog(*cat, maxStations);
  addEvents1ToCatalog(*cat, ttt, time + Core::TimeSpan(1), lat, lon, depth,
                      numEvents / 3, extent);
  addEvents2ToCatalog(*cat, ttt, time + Core::TimeSpan(2), lat, lon, depth,
                      numEvents / 3, extent);
  addEvents3ToCatalog(*cat, ttt, time + Core::TimeSpan(3), lat, lon, depth,
                      numEvents / 3, extent);
  return cat;
}

std::unique_ptr<HDD::Catalog> relocateCatalog(const HDD::Catalog& cat,
                                HDD::TravelTimeTable &ttt,
                                const string &workingDir)
{
  HDD::Config ddCfg;
  ddCfg.ttt.type  = ttt.type;
  ddCfg.ttt.model = ttt.model;

  HDD::DD dd(cat, ddCfg, workingDir);
  dd.setSaveProcessing(true);
  dd.setUseCatalogWaveformDiskCache(false);
  dd.setWaveformCacheAll(false);
  dd.setWaveformDebug(false);
  dd.setUseArtificialPhases(false);

  HDD::ClusteringOptions clusterCfg;
  clusterCfg.numEllipsoids    = 0;
  clusterCfg.maxEllipsoidSize = 100;
  // disable cross-correlation
  clusterCfg.xcorrMaxEvStaDist   = 0;
  clusterCfg.xcorrMaxInterEvDist = 0;

  HDD::SolverOptions solverCfg;
  solverCfg.algoIterations               = 20;
  solverCfg.absLocConstraintStart        = 1.0;
  solverCfg.absLocConstraintEnd          = 1.0;
  solverCfg.dampingFactorStart           = 0.;
  solverCfg.dampingFactorEnd             = 0.;
  solverCfg.downWeightingByResidualStart = 0;
  solverCfg.downWeightingByResidualEnd   = 0;

  std::unique_ptr<HDD::Catalog> relocCat = dd.relocateMultiEvents(clusterCfg, solverCfg);

  // comment this for debugging
  boost::filesystem::remove_all(workingDir);

  return relocCat;
}

unique_ptr<const HDD::Catalog> relocateSingleEvent(const HDD::Catalog& bgCat,
                                     HDD::TravelTimeTable &ttt,
                                     const string &workingDir,
                                     const HDD::Catalog &realTimeCat)
{
  HDD::Config ddCfg;
  ddCfg.ttt.type  = ttt.type;
  ddCfg.ttt.model = ttt.model;

  HDD::DD dd(bgCat, ddCfg, workingDir);
  dd.setSaveProcessing(true);
  dd.setUseCatalogWaveformDiskCache(false);
  dd.setWaveformCacheAll(false);
  dd.setWaveformDebug(false);
  dd.setUseArtificialPhases(false);

  HDD::ClusteringOptions clusterCfg;
  clusterCfg.numEllipsoids    = 5;
  clusterCfg.maxEllipsoidSize = 15;
  clusterCfg.maxNumNeigh      = 40;
  // disable cross-correlation
  clusterCfg.xcorrMaxEvStaDist   = 0;
  clusterCfg.xcorrMaxInterEvDist = 0;

  HDD::SolverOptions solverCfg;
  solverCfg.algoIterations               = 20;
  solverCfg.absLocConstraintStart        = 1.0;
  solverCfg.absLocConstraintEnd          = 1.0;
  solverCfg.dampingFactorStart           = 0.;
  solverCfg.dampingFactorEnd             = 0.;
  solverCfg.downWeightingByResidualStart = 0;
  solverCfg.downWeightingByResidualEnd   = 0;

  unique_ptr<HDD::Catalog> relocCat(new HDD::Catalog());
  for (const auto &kv : realTimeCat.getEvents())
  {
    const Event &ev                = kv.second;
    unique_ptr<HDD::Catalog> orgToRelocate  = realTimeCat.extractEvent(ev.id, false);
    unique_ptr<HDD::Catalog> relocatedEvent = dd.relocateSingleEvent(
        *orgToRelocate, clusterCfg, clusterCfg, solverCfg);
    relocCat->add(*relocatedEvent, false);
  }

  // comment this for debugging
  boost::filesystem::remove_all(workingDir);

  return relocCat;
}

void testCatalogEqual(const HDD::Catalog& cat1, const HDD::Catalog& cat2)
{
  for (const auto &kv : cat1.getEvents())
  {
    const Event &ev1 = kv.second;
    BOOST_CHECK_EQUAL(cat2.getEvents().count(ev1.id), 1);
    if (cat2.getEvents().count(ev1.id) != 1)
    {
      continue;
    }
    const Event &ev2 = cat2.getEvents().at(ev1.id);
    BOOST_CHECK_SMALL((ev1.time - ev2.time).length(),
                      0.01); // tolerance in seconds
    BOOST_CHECK_CLOSE(ev1.latitude, ev2.latitude, 1);
    BOOST_CHECK_CLOSE(ev1.longitude, ev2.longitude, 1);
    BOOST_CHECK_CLOSE(ev1.depth, ev2.depth, 1);
  }
}

} // namespace

BOOST_DATA_TEST_CASE(test_dd_multi_event, bdata::xrange(tttList.size()), tttIdx)
{
  // Logging::enableConsoleLogging(Logging::getAll());

  unique_ptr<HDD::TravelTimeTable> ttt =
      HDD::TravelTimeTable::create(tttList[tttIdx].type, tttList[tttIdx].model);

  const Core::Time clusterTime = Core::Time::FromString("2001-01-02", "%F");
  const double clusterLat      = 47.0;
  const double clusterLon      = 8.5;
  const double clusterDepth    = 5;

  unique_ptr<const HDD::Catalog> baseCat = buildCatalog(
      *ttt, 8, clusterTime, clusterLat, clusterLon, clusterDepth, 66, 1.0);

  // Test 1: no event changes
  string workingDir = strf("./data/test_dd_multi_event_%d_reloc1", tttIdx);
  unique_ptr<const HDD::Catalog> relocCat = relocateCatalog(*baseCat, *ttt, workingDir);
  testCatalogEqual(*baseCat, *relocCat);

  // Test 2: all events at the center location, depth 1 km
  unique_ptr<HDD::Catalog> cat(new HDD::Catalog(*baseCat));
  HDD::NormalRandomer timeDist(0.1, 0.400, 0x1004); // sec
  for (const auto &kv : cat->getEvents())
  {
    Event ev = kv.second;
    ev.time += Core::TimeSpan(timeDist.next());
    ev.latitude  = clusterLat;
    ev.longitude = clusterLon;
    ev.depth     = 1;
    cat->updateEvent(ev);
  }

  workingDir = strf("./data/test_dd_multi_event_%d_reloc2", tttIdx);
  relocCat   = relocateCatalog(*cat, *ttt, workingDir);
  testCatalogEqual(*baseCat, *relocCat);

  // Test 3: all events at the center location, depth 10 km
  cat.reset(new HDD::Catalog(*baseCat));
  for (const auto &kv : cat->getEvents())
  {
    Event ev = kv.second;
    ev.time += Core::TimeSpan(timeDist.next());
    ev.latitude  = clusterLat;
    ev.longitude = clusterLon;
    ev.depth     = 10;
    cat->updateEvent(ev);
  }

  workingDir = strf("./data/test_dd_multi_event_%d_reloc3", tttIdx);
  relocCat   = relocateCatalog(*cat, *ttt, workingDir);
  testCatalogEqual(*baseCat, *relocCat);

  // Test 4: random changes to all events (mean of all changes is != 0)
  cat.reset(new HDD::Catalog(*baseCat));
  HDD::NormalRandomer latDist(0.0055, 0.02, 0x1001);
  HDD::NormalRandomer lonDist(-0.011, 0.04, 0x1002);
  HDD::NormalRandomer depthDist(-0.6, 2.0, 0x1003); // km
  for (const auto &kv : cat->getEvents())
  {
    Event ev = kv.second;
    ev.time += Core::TimeSpan(timeDist.next());
    ev.latitude += latDist.next();
    ev.longitude += lonDist.next();
    ev.depth += depthDist.next();
    cat->updateEvent(ev);
  }

  workingDir = strf("./data/test_dd_multi_event_%d_reloc4", tttIdx);
  relocCat   = relocateCatalog(*cat, *ttt, workingDir);
  testCatalogEqual(*baseCat, *relocCat);
}

BOOST_DATA_TEST_CASE(test_dd_single_event,
                     bdata::xrange(tttList.size()),
                     tttIdx)
{
  // Logging::enableConsoleLogging(Logging::getAll());

  unique_ptr<HDD::TravelTimeTable> ttt =
      HDD::TravelTimeTable::create(tttList[tttIdx].type, tttList[tttIdx].model);

  const Core::Time clusterTime = Core::Time::FromString("2001-01-02", "%F");
  const double clusterLat      = 47.0;
  const double clusterLon      = 8.5;
  const double clusterDepth    = 5;

  unique_ptr<const HDD::Catalog> backgroundCat = buildBackgroundCatalog(
      *ttt, 8, clusterTime, clusterLat, clusterLon, clusterDepth, 21, 1.0);


  unique_ptr<const HDD::Catalog> realTimeCat = buildCatalog(
      *ttt, 8, clusterTime, clusterLat, clusterLon, clusterDepth, 66, 1.0);

  // Test 1: no event changes
  string workingDir =
      strf("./data/test_dd_single_event_%d_reloc1", tttIdx);
  unique_ptr<const HDD::Catalog> relocCat =
      relocateSingleEvent(*backgroundCat, *ttt, workingDir, *realTimeCat);
  testCatalogEqual(*realTimeCat, *relocCat);

  // Test 2: all events at the center location, depth 1 km
  std::unique_ptr<HDD::Catalog> cat(new HDD::Catalog(*realTimeCat));
  HDD::NormalRandomer timeDist(0.2, 0.500, 0x1004); // sec
  for (const auto &kv : cat->getEvents())
  {
    Event ev = kv.second;
    ev.time += Core::TimeSpan(timeDist.next());
    ev.latitude  = clusterLat;
    ev.longitude = clusterLon;
    ev.depth     = 1;
    cat->updateEvent(ev);
  }

  workingDir = strf("./data/test_dd_single_event_%d_reloc2", tttIdx);
  relocCat   = relocateSingleEvent(*backgroundCat, *ttt, workingDir, *realTimeCat);
  testCatalogEqual(*realTimeCat, *relocCat);

  // Test 3: all events at the center location, depth 10 km
  cat.reset(new HDD::Catalog(*realTimeCat));
  for (const auto &kv : cat->getEvents())
  {
    Event ev = kv.second;
    ev.time += Core::TimeSpan(timeDist.next());
    ev.latitude  = clusterLat;
    ev.longitude = clusterLon;
    ev.depth     = 10;
    cat->updateEvent(ev);
  }

  workingDir = strf("./data/test_dd_single_event_%d_reloc3", tttIdx);
  relocCat   = relocateSingleEvent(*backgroundCat, *ttt, workingDir, *realTimeCat);
  testCatalogEqual(*realTimeCat, *relocCat);

  // Test 4: random changes to all events (mean of all changes is != 0)
  cat.reset(new HDD::Catalog(*realTimeCat));
  HDD::NormalRandomer latDist(0.006, 0.02, 0x1001);
  HDD::NormalRandomer lonDist(-0.012, 0.04, 0x1002);
  HDD::NormalRandomer depthDist(-0.6, 2.0, 0x1003); // km
  for (const auto &kv : cat->getEvents())
  {
    Event ev = kv.second;
    ev.time += Core::TimeSpan(timeDist.next());
    ev.latitude += latDist.next();
    ev.longitude += lonDist.next();
    ev.depth += depthDist.next();
    cat->updateEvent(ev);
  }

  workingDir = strf("./data/test_dd_single_event_%d_reloc4", tttIdx);
  relocCat   = relocateSingleEvent(*backgroundCat, *ttt, workingDir, *realTimeCat);
  testCatalogEqual(*realTimeCat, *relocCat);
}

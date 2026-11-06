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
                             const string &workingDir,
                             double downWeightingStart = 10,
                             double downWeightingEnd   = 5,
                             double absLocConstraint   = 0.3)
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
  solverOpt.absLocConstraintStart        = absLocConstraint;
  solverOpt.absLocConstraintEnd          = absLocConstraint;
  solverOpt.dampingFactorStart           = 0.01;
  solverOpt.dampingFactorEnd             = 0.01;
  solverOpt.downWeightingByResidualStart = downWeightingStart;
  solverOpt.downWeightingByResidualEnd   = downWeightingEnd;

  std::list<std::unordered_map<unsigned, Neighbours>> clusters;
  XCorrCache xcorrData;
  return dd.relocateMultiEvents(clusterOpt, clusters, xcorrData, xcorrOpt,
                                solverOpt, saveProcessing);
}

// Residual down-weighting is disabled by default: with noise-free picks the
// residuals of a single event are not random noise and the down-weighting
// prevents some events from converging to the true location
HDD::Catalog relocateSingleEvent(const HDD::Catalog &bgCat,
                                 unique_ptr<HDD::TravelTimeTable> ttt,
                                 const string &workingDir,
                                 const HDD::Catalog &realTimeCat,
                                 double downWeightingStart = 0,
                                 double downWeightingEnd   = 0)
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
  solverOpt.downWeightingByResidualStart = downWeightingStart;
  solverOpt.downWeightingByResidualEnd   = downWeightingEnd;

  HDD::Catalog relocCat;
  for (const auto &kv : realTimeCat.getEvents())
  {
    const Event &ev            = kv.second;
    HDD::Catalog orgToRelocate = realTimeCat.extractEvent(ev.id, false);
    HDD::Catalog relocatedEvent =
        dd.relocateSingleEvent(orgToRelocate, true, clusterOpt, xcorrOpt,
                               solverOpt, saveProcessing);
    if (relocatedEvent.empty()) continue;

    // keep the original event id, so that the relocated events can be
    // compared with the input ones even if some relocations fail
    Event relocEv = relocatedEvent.getEvents().begin()->second;
    unordered_multimap<unsigned, Phase> phases;
    for (const auto &kv2 : relocatedEvent.getPhases())
    {
      Phase ph   = kv2.second;
      ph.eventId = ev.id;
      phases.emplace(ev.id, ph);
    }
    relocEv.id = ev.id;
    relocCat.add(HDD::Catalog(relocatedEvent.getStations(),
                              {{ev.id, relocEv}}, std::move(phases)),
                 true);
  }

  return relocCat;
}

HDD::Catalog relocateCatalogAgainstReference(const HDD::Catalog &refCat,
                                             unique_ptr<HDD::TravelTimeTable> ttt,
                                             const string &workingDir,
                                             const HDD::Catalog &cat,
                                             double downWeightingStart = 10,
                                             double downWeightingEnd   = 5,
                                             double absLocConstraint   = 0.3,
                                             bool saveProcessing = false)
{
  HDD::Config ddCfg;

  HDD::DD dd(refCat, ddCfg, std::move(ttt));

  HDD::ClusteringOptions clusterOpt;
  clusterOpt.numEllipsoids    = 0;
  clusterOpt.maxNeighbourDist = 15;
  clusterOpt.minEvStaToInterEvRatio   = 0;

  HDD::XcorrOptions xcorrOpt;
  xcorrOpt.enable = false; // disable cross-correlation

  HDD::SolverOptions solverOpt;
  solverOpt.algoIterations               = 20;
  solverOpt.absLocConstraintStart        = absLocConstraint;
  solverOpt.absLocConstraintEnd          = absLocConstraint;
  solverOpt.dampingFactorStart           = 0.01;
  solverOpt.dampingFactorEnd             = 0.01;
  solverOpt.downWeightingByResidualStart = downWeightingStart;
  solverOpt.downWeightingByResidualEnd   = downWeightingEnd;

  std::list<std::unordered_map<unsigned, Neighbours>> clusters;
  XCorrCache xcorrData;
  HDD::Catalog relocCat =
      dd.relocateMultiEvents(cat, clusterOpt, clusters, xcorrData, xcorrOpt,
                             solverOpt, saveProcessing, workingDir);
  if (saveProcessing) removePath(workingDir);

  // reference catalog events are fixed: they must not be relocated
  BOOST_CHECK(!clusters.empty());
  for (const auto &cluster : clusters)
  {
    for (const auto &kv : cluster)
    {
      BOOST_CHECK_EQUAL(refCat.getEvents().count(kv.first), 0);
    }
  }
  return relocCat;
}

// Move all events by the same amount: `north`, `east` and `down` [km] and
// `time` [sec]
void shiftCatalog(HDD::Catalog &cat,
                  double north,
                  double east,
                  double down,
                  double time)
{
  for (const auto &kv : cat.getEvents())
  {
    Event ev = kv.second;
    ev.time += secToDur(time);
    ev.depth += down;
    computeCoordinates(std::abs(north), degToRad(north >= 0 ? 0 : 180),
                       ev.latitude, ev.longitude, ev.latitude, ev.longitude,
                       ev.depth);
    computeCoordinates(std::abs(east), degToRad(east >= 0 ? 90 : 270),
                       ev.latitude, ev.longitude, ev.latitude, ev.longitude,
                       ev.depth);
    cat.updateEvent(ev);
  }
}

// Add gaussian noise (`stdDev` sec) to all pick times. In addition, a
// fraction `outlierRatio` of the picks is shifted by +/-`outlierError` sec
HDD::Catalog addPickNoise(const HDD::Catalog &cat,
                          double stdDev,
                          double outlierRatio,
                          double outlierError,
                          unsigned seed)
{
  HDD::NormalRandomer noise(0, stdDev, seed);
  HDD::UniformRandomer outlier(0, 9999, seed + 1);
  unordered_multimap<unsigned, Phase> phases;
  for (const auto &kv : cat.getEvents())
  {
    auto eqlrng = cat.getPhases().equal_range(kv.first);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      Phase ph = it->second;
      double error = noise.next();
      size_t draw  = outlier.next();
      if (draw < outlierRatio * 10000)
      {
        error += (draw % 2 == 0) ? outlierError : -outlierError;
      }
      ph.time += secToDur(error);
      phases.emplace(ph.eventId, ph);
    }
  }
  return HDD::Catalog(cat.getStations(), cat.getEvents(), std::move(phases));
}

struct RelocationError
{
  size_t numEvents = 0; // events in `cat2` matched in `cat1`
  struct Stats
  {
    double median, p90, max;
  } time, hdist, depth; // [sec], [km], [km]
};

// Compare the events in `cat2` with the events with same id in `cat1`
RelocationError computeRelocationError(const HDD::Catalog &cat1,
                                       const HDD::Catalog &cat2)
{
  vector<double> times, hdists, depths;
  for (const auto &kv : cat2.getEvents())
  {
    const Event &ev2 = kv.second;
    auto it          = cat1.getEvents().find(ev2.id);
    if (it == cat1.getEvents().end()) continue;
    const Event &ev1 = it->second;
    times.push_back(std::abs(durToSec(ev1.time - ev2.time)));
    hdists.push_back(computeDistance(ev1.latitude, ev1.longitude,
                                     ev2.latitude, ev2.longitude));
    depths.push_back(std::abs(ev1.depth - ev2.depth));
  }

  auto stats = [](vector<double> values) -> RelocationError::Stats {
    std::sort(values.begin(), values.end());
    size_t p90Idx = std::min(values.size() - 1, size_t(0.9 * values.size()));
    return {computeMedian(values), values.at(p90Idx), values.back()};
  };

  RelocationError err;
  err.numEvents = times.size();
  if (err.numEvents == 0) return err;
  err.time  = stats(times);
  err.hdist = stats(hdists);
  err.depth = stats(depths);
  return err;
}

struct NoiseTolerance
{
  double medianTime, medianHDist, medianDepth; // [sec], [km], [km]
  double p90Time, p90HDist, p90Depth;          // [sec], [km], [km]
};

//
// `relocCat` is the relocation of a noisy catalog with residual
// down-weighting enabled, `relocCatNoDw` the same relocation without
// down-weighting. Check that the relocation errors are compatible with the
// pick noise and that the down-weighting reduced the effect of the outliers
//
void testNoisyRelocation(const HDD::Catalog &baseCat,
                         const HDD::Catalog &relocCat,
                         const HDD::Catalog &relocCatNoDw,
                         const NoiseTolerance &tol)
{
  const RelocationError err     = computeRelocationError(baseCat, relocCat);
  const RelocationError errNoDw = computeRelocationError(baseCat, relocCatNoDw);

  BOOST_TEST_MESSAGE(strf(
      "Relocation error with (without) down-weighting: time median %.1f (%.1f) "
      "p90 %.1f (%.1f) max %.1f (%.1f) msec, horizontal median %.1f (%.1f) "
      "p90 %.1f (%.1f) max %.1f (%.1f) m, depth median %.1f (%.1f) p90 %.1f "
      "(%.1f) max %.1f (%.1f) m",
      err.time.median * 1e3, errNoDw.time.median * 1e3, err.time.p90 * 1e3,
      errNoDw.time.p90 * 1e3, err.time.max * 1e3, errNoDw.time.max * 1e3,
      err.hdist.median * 1e3, errNoDw.hdist.median * 1e3,
      err.hdist.p90 * 1e3, errNoDw.hdist.p90 * 1e3, err.hdist.max * 1e3,
      errNoDw.hdist.max * 1e3, err.depth.median * 1e3,
      errNoDw.depth.median * 1e3, err.depth.p90 * 1e3,
      errNoDw.depth.p90 * 1e3, err.depth.max * 1e3, errNoDw.depth.max * 1e3));

  // all events must be relocated
  BOOST_REQUIRE_EQUAL(err.numEvents, baseCat.getEvents().size());
  BOOST_REQUIRE_EQUAL(errNoDw.numEvents, baseCat.getEvents().size());

  // Only robust statistics are checked: occasionally a single event can
  // diverge (e.g. the homogeneous travel time table allows air-quakes)
  BOOST_CHECK_LE(err.time.median, tol.medianTime);
  BOOST_CHECK_LE(err.hdist.median, tol.medianHDist);
  BOOST_CHECK_LE(err.depth.median, tol.medianDepth);
  BOOST_CHECK_LE(err.time.p90, tol.p90Time);
  BOOST_CHECK_LE(err.hdist.p90, tol.p90HDist);
  BOOST_CHECK_LE(err.depth.p90, tol.p90Depth);

  // the down-weighting must reduce the effect of the outliers
  BOOST_CHECK_GE(errNoDw.hdist.median, 2 * err.hdist.median);
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

// Noisy data: gaussian pick noise plus a fraction of outlier picks
const double pickNoise    = 0.010; // [sec] std dev
const double outlierRatio = 0.05;  // fraction of picks
const double outlierError = 0.5;   // [sec]

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

// All the following tests run for every centroid in `centroidList` (which
// includes locations close to the poles and to the antimeridian) and every
// travel time table in `tttList`. The NonLinLoc tables are tested only at the
// locations covered by the grid files. The catalogs are synthetic and their
// picks are computed with the same travel time table used for the
// relocation: unless noise is added to the picks, a correct relocation
// recovers the true event locations and origin times.

// Multi-event relocation of a catalog whose locations are already correct:
// verify that the relocation does not move the events
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

// Multi-event relocation starting with all events at the same location and
// with perturbed origin times: verify that the true locations and origin
// times are recovered
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

// Multi-event relocation of a catalog whose events have random location and
// time errors (with a non-zero mean): verify that the true locations and
// origin times are recovered
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

// Single-event relocation, one event at the time, against the background
// catalog of events whose locations are already correct: verify that the
// relocation does not move the events
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

// Single-event relocation, one event at the time, against the background
// catalog of events starting at the same location and with perturbed origin
// times: verify that the true locations and origin times are recovered
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
      backgroundCat, std::move(ttt), workingDir, cat);
  testCatalogEqual(realTimeCat, relocCat);
}

// Single-event relocation, one event at the time, against the background
// catalog of events with random location and time errors (with a non-zero
// mean): verify that the true locations and origin times are recovered
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
      backgroundCat, std::move(ttt), workingDir, cat);
  testCatalogEqual(realTimeCat, relocCat);
}

// Multi-event relocation against a fixed reference catalog of events whose
// locations are already correct: verify that the relocation does not move the
// events. As in all the reference catalog tests, verify also that the
// reference events are never relocated and that the relocated events keep
// their original ids (they overlap with the reference ids, so they are
// renumbered internally)
BOOST_DATA_TEST_CASE(test_dd_multi_event_reference1,
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
      UTCClock::fromDate(1981, 2, 13, 7, 21, 11, 54321);
  const HDD::Catalog referenceCat = buildBackgroundCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 21, 2.0);
  const HDD::Catalog baseCat = buildCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 24, 2.0);

  // Test 1: no event changes
  string workingDir =
      strf("./data/test_dd_multi_event_reference1_%lu_%lu_reloc", cIdx, tttIdx);
  const HDD::Catalog relocCat = relocateCatalogAgainstReference(
      referenceCat, std::move(ttt), workingDir, baseCat);
  testCatalogEqual(baseCat, relocCat);
}

// As test_dd_multi_event_reference1, with a reference catalog whose events are
// flagged as relocated (e.g. the output of a previous relocation or a catalog
// loaded with its relocation information): the reference events must not be
// returned, with and without saving the processing files
BOOST_DATA_TEST_CASE(test_dd_multi_event_reference_relocated,
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
      UTCClock::fromDate(1981, 2, 13, 7, 21, 11, 54321);
  HDD::Catalog referenceCat = buildBackgroundCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 21, 2.0);
  for (const auto &kv : referenceCat.getEvents())
  {
    Event ev                 = kv.second;
    ev.relocInfo.isRelocated = true;
    referenceCat.updateEvent(ev);
  }
  const HDD::Catalog baseCat = buildCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 24, 2.0);

  for (bool saveProcessing : {false, true})
  {
    string workingDir =
        strf("./data/test_dd_multi_event_reference_relocated_%lu_%lu_reloc",
             cIdx, tttIdx);
    const HDD::Catalog relocCat = relocateCatalogAgainstReference(
        referenceCat, createTTT(tttList.at(tttIdx)), workingDir, baseCat, 10,
        5, 0.3, saveProcessing);
    BOOST_CHECK_EQUAL(relocCat.getEvents().size(),
                      baseCat.getEvents().size());
    testCatalogEqual(baseCat, relocCat);
  }
}

// Multi-event relocation against a fixed reference catalog of events starting
// at the same location and with perturbed origin times: verify that the true
// locations and origin times are recovered
BOOST_DATA_TEST_CASE(test_dd_multi_event_reference2,
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
      UTCClock::fromDate(1996, 9, 3, 18, 0, 41, 900001);
  const HDD::Catalog referenceCat = buildBackgroundCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 21, 1.5);
  const HDD::Catalog baseCat = buildCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 24, 1.5);

  // Test 2: all events at centroid location, with random perturbation of time
  HDD::Catalog cat(baseCat);
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
      strf("./data/test_dd_multi_event_reference2_%lu_%lu_reloc", cIdx, tttIdx);
  const HDD::Catalog relocCat = relocateCatalogAgainstReference(
      referenceCat, std::move(ttt), workingDir, cat);
  testCatalogEqual(baseCat, relocCat);
}

// Multi-event relocation against a fixed reference catalog of events with
// random location and time errors (with a non-zero mean): verify that the
// true locations and origin times are recovered
BOOST_DATA_TEST_CASE(test_dd_multi_event_reference3,
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

  const UTCTime clusterTime  = UTCClock::fromDate(2027, 1, 30, 9, 12, 5, 777);
  const HDD::Catalog referenceCat = buildBackgroundCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 21, 1.0);
  const HDD::Catalog baseCat = buildCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 24, 1.0);

  // Test 3: random changes to all events (mean of all changes is != 0)
  HDD::Catalog cat(baseCat);
  randomPerturbation(cat);

  string workingDir =
      strf("./data/test_dd_multi_event_reference3_%lu_%lu_reloc", cIdx, tttIdx);
  const HDD::Catalog relocCat = relocateCatalogAgainstReference(
      referenceCat, std::move(ttt), workingDir, cat);
  testCatalogEqual(baseCat, relocCat);
}


// Multi-event relocation of events with random location and time errors and
// noisy picks (gaussian noise plus outliers): verify that the median and 90th
// percentile of the relocation errors are compatible with the pick noise and
// that the residual down-weighting reduces the effect of the outliers
BOOST_DATA_TEST_CASE(test_dd_multi_event_noise,
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
      buildCatalog(*ttt, nllStations, UTCClock::fromDate(2011, 3, 1, 4, 5, 6, 789),
                   centroid.lat, centroid.lon, centroid.depth, 42, 1.0);

  // Test noisy data: random changes to all events and noisy picks
  HDD::Catalog cat =
      addPickNoise(baseCat, pickNoise, outlierRatio, outlierError, 0x2001);
  randomPerturbation(cat);

  string workingDir =
      strf("./data/test_dd_multi_event_noise_%lu_%lu_reloc", cIdx, tttIdx);
  const HDD::Catalog relocCat =
      relocateCatalog(cat, std::move(ttt), workingDir, 10, 5);
  const HDD::Catalog relocCatNoDw =
      relocateCatalog(cat, createTTT(tttList.at(tttIdx)), workingDir, 0, 0);
  testNoisyRelocation(baseCat, relocCat, relocCatNoDw,
                      {0.020, 0.100, 0.500, 0.040, 0.150, 1.000});
}

// As test_dd_multi_event_noise, but relocating the events against a fixed
// reference catalog with noisy picks
BOOST_DATA_TEST_CASE(test_dd_multi_event_reference_noise,
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

  const UTCTime clusterTime = UTCClock::fromDate(2011, 3, 1, 4, 5, 6, 789);
  const HDD::Catalog referenceCat = addPickNoise(
      buildBackgroundCatalog(*ttt, nllStations, clusterTime, centroid.lat,
                             centroid.lon, centroid.depth, 21, 1.0),
      pickNoise, outlierRatio, outlierError, 0x3001);
  const HDD::Catalog baseCat = buildCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 24, 1.0);

  // Test noisy data: random changes to all events and noisy picks (the
  // reference catalog picks are noisy too)
  HDD::Catalog cat =
      addPickNoise(baseCat, pickNoise, outlierRatio, outlierError, 0x2001);
  randomPerturbation(cat);

  string workingDir =
      strf("./data/test_dd_multi_event_reference_noise_%lu_%lu_reloc", cIdx, tttIdx);
  const HDD::Catalog relocCat =
      relocateCatalogAgainstReference(referenceCat, std::move(ttt), workingDir, cat, 10, 5);
  const HDD::Catalog relocCatNoDw = relocateCatalogAgainstReference(
      referenceCat, createTTT(tttList.at(tttIdx)), workingDir, cat, 0, 0);
  testNoisyRelocation(baseCat, relocCat,
                      relocCatNoDw, {0.015, 0.050, 0.300, 0.030, 0.075, 0.800});
}

// As test_dd_multi_event_noise, but relocating the events one at the time
// (single-event) against a background catalog with noisy picks
BOOST_DATA_TEST_CASE(test_dd_single_event_noise,
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

  const UTCTime clusterTime = UTCClock::fromDate(2011, 3, 1, 4, 5, 6, 789);
  const HDD::Catalog referenceCat = addPickNoise(
      buildBackgroundCatalog(*ttt, nllStations, clusterTime, centroid.lat,
                             centroid.lon, centroid.depth, 21, 1.0),
      pickNoise, outlierRatio, outlierError, 0x3001);
  const HDD::Catalog baseCat = buildCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 24, 1.0);

  // Test noisy data: random changes to all events and noisy picks (the
  // reference catalog picks are noisy too)
  HDD::Catalog cat =
      addPickNoise(baseCat, pickNoise, outlierRatio, outlierError, 0x2001);
  randomPerturbation(cat);

  string workingDir =
      strf("./data/test_dd_single_event_noise_%lu_%lu_reloc", cIdx, tttIdx);
  const HDD::Catalog relocCat =
      relocateSingleEvent(referenceCat, std::move(ttt), workingDir, cat, 10, 5);
  const HDD::Catalog relocCatNoDw = relocateSingleEvent(
      referenceCat, createTTT(tttList.at(tttIdx)), workingDir, cat, 0, 0);
  testNoisyRelocation(baseCat, relocCat,
                      relocCatNoDw, {0.015, 0.050, 0.300, 0.030, 0.075, 0.800});
}

// Verify that the reference catalog anchors the absolute location of the
// relocated events: a common origin time shift, which cannot be resolved by
// the double-differences alone, is recovered when relocating against the
// reference catalog, but not by a standard multi-event relocation
BOOST_DATA_TEST_CASE(test_dd_multi_event_reference_anchor,
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

  const UTCTime clusterTime = UTCClock::fromDate(2019, 8, 7, 6, 5, 4, 321);
  const HDD::Catalog referenceCat = buildBackgroundCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 21, 1.0);
  const HDD::Catalog baseCat = buildCatalog(
      *ttt, nllStations, clusterTime, centroid.lat, centroid.lon, centroid.depth, 24, 1.0);

  // Test the reference catalog acts as an anchor: all events are shifted by
  // the same amount and the absolute location constraint is disabled, so only
  // the double-differences are used. A common origin time shift cancels out
  // in the double-differences and cannot be recovered by a standard
  // multi-event relocation, while the fixed reference events allow the
  // relocation against the reference catalog to recover it
  const double timeShift = 0.15; // sec
  HDD::Catalog cat(baseCat);
  shiftCatalog(cat, 1.0, -0.8, 1.0, timeShift);

  string workingDir =
      strf("./data/test_dd_multi_event_reference_anchor_%lu_%lu_reloc", cIdx, tttIdx);
  const HDD::Catalog relocCat = relocateCatalogAgainstReference(
      referenceCat, std::move(ttt), workingDir, cat, 10, 5, 0);
  testCatalogEqual(baseCat, relocCat);

  const HDD::Catalog relocCatNoRef = relocateCatalog(
      cat, createTTT(tttList.at(tttIdx)), workingDir, 10, 5, 0);
  const RelocationError errNoRef = computeRelocationError(baseCat, relocCatNoRef);
  BOOST_REQUIRE_EQUAL(errNoRef.numEvents, baseCat.getEvents().size());
  BOOST_CHECK_GE(errNoRef.time.median, 0.9 * timeShift);
}

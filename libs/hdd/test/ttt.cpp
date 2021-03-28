#define SEISCOMP_TEST_MODULE hdd
#include <seiscomp/unittest/unittests.h>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>

#include "catalog.h"
#include "ttt.h"

#include <seiscomp3/math/math.h>
#include <vector>

using namespace std;
using namespace Seiscomp;
using Seiscomp::Core::stringify;
namespace bdata = boost::unit_test::data;

namespace {

vector<HDD::TravelTimeTablePtr> tttList = {
    HDD::TravelTimeTable::create("LOCSAT", "iasp91"),
    HDD::TravelTimeTable::create("libtau", "iasp91"),
    HDD::TravelTimeTable::create(
        "NonLinLoc",
        "./data/nll/iasp91_2D_simple/model/iasp91.PHASE.mod;"
        "./data/nll/iasp91_2D_simple/time/iasp91.PHASE.STATION.time;"
        "./data/nll/iasp91_2D_simple/time/iasp91.PHASE.STATION.angle"),
    HDD::TravelTimeTable::create(
        "NonLinLoc",
        "./data/nll/iasp91_2D_sdc/model/iasp91.PHASE.mod;"
        "./data/nll/iasp91_2D_sdc/time/iasp91.PHASE.STATION.time;"
        "./data/nll/iasp91_2D_sdc/time/iasp91.PHASE.STATION.angle"),
    HDD::TravelTimeTable::create(
        "NonLinLoc",
        "./data/nll/iasp91_2D_global/model/iasp91.PHASE.mod;"
        "./data/nll/iasp91_2D_global/time/iasp91.PHASE.STATION.time;"
        "./data/nll/iasp91_2D_global/time/iasp91.PHASE.STATION.angle"),
    HDD::TravelTimeTable::create(
        "NonLinLoc",
        "./data/nll/iasp91_3D_simple/model/iasp91.PHASE.mod;"
        "./data/nll/iasp91_3D_simple/time/iasp91.PHASE.STATION.time;"
        "./data/nll/iasp91_3D_simple/time/iasp91.PHASE.STATION.angle"),
    HDD::TravelTimeTable::create(
        "NonLinLoc",
        "./data/nll/iasp91_3D_sdc/model/iasp91.PHASE.mod;"
        "./data/nll/iasp91_3D_sdc/time/iasp91.PHASE.STATION.time;"
        "./data/nll/iasp91_3D_sdc/time/iasp91.PHASE.STATION.angle")};

vector<HDD::Catalog::Station> stationList = {
    {"NET.ST01A", 47.1, 8.6, 250, "NET", "ST01A", ""},
    {"NET.ST02A", 47.1, 8.4, 295, "NET", "ST02A", ""},
    {"NET.ST03A", 46.9, 8.4, 301, "NET", "ST03A", ""},
    {"NET.ST04A", 46.9, 8.6, 395, "NET", "ST04A", ""},
    {"NET.ST01B", 47.0, 8.7, 212, "NET", "ST01B", ""},
    {"NET.ST02B", 47.0, 8.3, 346, "NET", "ST02B", ""},
    {"NET.ST03B", 47.2, 8.5, 351, "NET", "ST03B", ""},
    {"NET.ST04B", 46.8, 8.5, 268, "NET", "ST04B", ""},
};

struct Delta
{
  double lat;
  double lon;
  double depth;
};

vector<Delta> deltaList = {
  //
  {0.05, 0, 1}, {0.05, 0, 2}, {0.05, 0, 4}, {0.05, 0, 8}, {0.05, 0, 16},
  {0, 0.05, 1}, {0, 0.05, 2}, {0, 0.05, 4}, {0, 0.05, 8}, {0, 0.05, 16},
  {-0.05, 0, 1}, {-0.05, 0, 2}, {-0.05, 0, 4}, {-0.05, 0, 8}, {-0.05, 0,
  16}, {0, -0.05, 1}, {0, -0.05, 2}, {0, -0.05, 4}, {0, -0.05, 8}, {0,
  -0.05, 16}, {0.05, 0.05, 1}, {0.05, 0.05, 2}, {0.05, 0.05, 4}, {0.05,
  0.05, 8}, {0.05, 0.05, 16},
  {-0.05, -0.05, 1}, {-0.05, -0.05, 2}, {-0.05, -0.05, 4}, {-0.05, -0.05,
  8}, {-0.05, -0.05, 16},
  //
  {0.1, 0, 1}, {0.1, 0, 2}, {0.1, 0, 4}, {0.1, 0, 8}, {0.1, 0, 16},
  {0, 0.1, 1}, {0, 0.1, 2}, {0, 0.1, 4}, {0, 0.1, 8}, {0, 0.1, 16},
  {-0.1, 0, 1}, {-0.1, 0, 2}, {-0.1, 0, 4}, {-0.1, 0, 8}, {-0.1, 0, 16},
  {0, -0.1, 1}, {0, -0.1, 2}, {0, -0.1, 4}, {0, -0.1, 8}, {0, -0.1, 16},
  //
  {0.4, 0.4, 1}, {0.4, 0.4, 2}, {0.4, 0.4, 4}, {0.4, 0.4, 8}, {0.4, 0.4,
  16},
  {-0.4, -0.4, 1}, {-0.4, -0.4, 2}, {-0.4, -0.4, 4}, {-0.4, -0.4, 8}, {-0.4,
  -0.4, 16},
  //
  {0.9, 0, 1}, {0.9, 0, 2}, {0.9, 0, 4}, {0.9, 0, 8}, {0.9, 0, 16},
  {0, 0.9, 1}, {0, 0.9, 2}, {0, 0.9, 4}, {0, 0.9, 8}, {0, 0.9, 16},
  {-0.9, 0, 1}, {-0.9, 0, 2}, {-0.9, 0, 4}, {-0.9, 0, 8}, {-0.9, 0, 16},
  {0, -0.9, 1}, {0, -0.9, 2}, {0, -0.9, 4}, {0, -0.9, 8}, {0, -0.9, 16},
};

} // namespace

BOOST_DATA_TEST_CASE(test_ttt, bdata::xrange(deltaList.size()), deltaIdx)
{
  const Delta &delta = deltaList[deltaIdx];

  BOOST_TEST_MESSAGE(stringify("Testing DELTA lat %.1f lon %.1f depth %.1f",
                               delta.lat, delta.lon, delta.depth));

  for (auto station : stationList)
  {
    const double stationDepth = -(station.elevation / 1000.);
    const double lat          = station.latitude + delta.lat;
    const double lon          = station.longitude + delta.lon;
    const double depth        = stationDepth + delta.depth;

    BOOST_TEST_MESSAGE(stringify(
        "Testing station %s lat %.1f lon %.1f depth %.1f", station.id.c_str(),
        station.latitude, station.longitude, stationDepth));

    vector<double> travelTimeP(tttList.size(), 0);
    vector<double> takeOffAngleAzimP(tttList.size(), 0);
    vector<double> takeOffAngleDipP(tttList.size(), 0);
    vector<double> velocityAtSrcP(tttList.size(), 0);

    vector<double> travelTimeS(tttList.size(), 0);
    vector<double> takeOffAngleAzimS(tttList.size(), 0);
    vector<double> takeOffAngleDipS(tttList.size(), 0);
    vector<double> velocityAtSrcS(tttList.size(), 0);

    for (size_t i = 0; i < tttList.size(); i++)
    {
      HDD::TravelTimeTablePtr ttt = tttList[i];
      BOOST_REQUIRE(ttt);
      BOOST_CHECK_NO_THROW(ttt->compute(
          lat, lon, depth, station, "P", travelTimeP[i], takeOffAngleAzimP[i],
          takeOffAngleDipP[i], velocityAtSrcP[i]));
      BOOST_CHECK_NO_THROW(ttt->compute(
          lat, lon, depth, station, "S", travelTimeS[i], takeOffAngleAzimS[i],
          takeOffAngleDipS[i], velocityAtSrcS[i]));
    }

    for (size_t i = 0; i < tttList.size(); i++)
    {
      const HDD::TravelTimeTableCPtr &ttt = tttList[i];
      BOOST_TEST_MESSAGE(stringify(
          "TTT type %-9s [Travel time, Velocity, Take-Off Angle Azimuth and "
          "Dip] "
          "P [%5.2f %5.2f %4.f %4.f] S [%5.2f %5.2f %4.f %4.f]",
          ttt->type.c_str(), travelTimeP[i], velocityAtSrcP[i],
          rad2deg(takeOffAngleAzimP[i]), rad2deg(takeOffAngleDipP[i]),
          travelTimeS[i], velocityAtSrcS[i], rad2deg(takeOffAngleAzimS[i]),
          rad2deg(takeOffAngleDipS[i])));
    }

    for (size_t i = 0; i < tttList.size() - 1; i++)
    {
      for (size_t j = i + 1; j < tttList.size(); j++)
      {
        BOOST_CHECK_CLOSE(travelTimeP[i], travelTimeP[j], 15);
        BOOST_CHECK_CLOSE(travelTimeS[i], travelTimeS[j], 15);
        BOOST_CHECK_CLOSE(velocityAtSrcP[i], velocityAtSrcP[j], 1);
        BOOST_CHECK_CLOSE(velocityAtSrcS[i], velocityAtSrcS[j], 1);

        auto checkCloseAngles = [](double x, double y,
                                   double degreeTol) {
          double diffAngle = atan2(sin(x-y), cos(x-y));
          diffAngle = rad2deg(diffAngle);
          BOOST_TEST_MESSAGE(" x " << rad2deg(x) << " y " << rad2deg(y) << " delta " << diffAngle);
          BOOST_CHECK_SMALL(diffAngle, degreeTol);
        };

        checkCloseAngles(takeOffAngleAzimP[i], takeOffAngleAzimP[j], 30);
        checkCloseAngles(takeOffAngleDipP[i], takeOffAngleDipP[j], 30);
        checkCloseAngles(takeOffAngleAzimS[i], takeOffAngleAzimS[j], 30);
        checkCloseAngles(takeOffAngleDipS[i], takeOffAngleDipS[j], 30);
      }
    }
  }
}

#define SEISCOMP_TEST_MODULE hdd
#include <seiscomp/unittest/unittests.h>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp> 

#include "catalog.h"
#include "ttt.h"

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
        "./data/nll/iasp91_3D_sdc/model/iasp91.PHASE.mod;"
        "./data/nll/iasp91_3D_sdc/time/iasp91.PHASE.STATION.time;"
        "./data/nll/iasp91_3D_sdc/time/iasp91.PHASE.STATION.angle"),
    HDD::TravelTimeTable::create(
        "NonLinLoc",
        "./data/nll/iasp91_3D_simple/model/iasp91.PHASE.mod;"
        "./data/nll/iasp91_3D_simple/time/iasp91.PHASE.STATION.time;"
        "./data/nll/iasp91_3D_simple/time/iasp91.PHASE.STATION.angle"),
    HDD::TravelTimeTable::create(
        "NonLinLoc",
        "./data/nll/iasp91_2D_global/model/iasp91.PHASE.mod;"
        "./data/nll/iasp91_2D_global/time/iasp91.PHASE.STATION.time;"
        "./data/nll/iasp91_2D_global/time/iasp91.PHASE.STATION.angle"),
    HDD::TravelTimeTable::create(
        "NonLinLoc",
        "./data/nll/iasp91_3D_global/model/iasp91.PHASE.mod;"
        "./data/nll/iasp91_3D_global/time/iasp91.PHASE.STATION.time;"
        "./data/nll/iasp91_3D_global/time/iasp91.PHASE.STATION.angle")
};

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
    {0, 0, 1}, {0, 0, 2}, {0, 0, 4}, {0, 0, 8}, {0, 0, 16},
    //
    {0.05, 0, 1}, {0.05, 0, 2}, {0.05, 0, 4}, {0.05, 0, 8}, {0.05, 0, 16},
    {0, 0.05, 1}, {0, 0.05, 2}, {0, 0.05, 4}, {0, 0.05, 8}, {0, 0.05, 16},
    {-0.05, 0, 1}, {-0.05, 0, 2}, {-0.05, 0, 4}, {-0.05, 0, 8}, {-0.05, 0, 16},
    {0, -0.05, 1}, {0, -0.05, 2}, {0, -0.05, 4}, {0, -0.05, 8}, {0, -0.05, 16},
    {0.05, 0.05, 1}, {0.05, 0.05, 2}, {0.05, 0.05, 4}, {0.05, 0.05, 8}, {0.05, 0.05, 16},
    {-0.05, -0.05, 1}, {-0.05, -0.05, 2}, {-0.05, -0.05, 4}, {-0.05, -0.05, 8}, {-0.05, -0.05, 16},
    //
    {0.1, 0, 1}, {0.1, 0, 2}, {0.1, 0, 4}, {0.1, 0, 8}, {0.1, 0, 16},
    {0, 0.1, 1}, {0, 0.1, 2}, {0, 0.1, 4}, {0, 0.1, 8}, {0, 0.1, 16},
    {-0.1, 0, 1}, {-0.1, 0, 2}, {-0.1, 0, 4}, {-0.1, 0, 8}, {-0.1, 0, 16},
    {0, -0.1, 1}, {0, -0.1, 2}, {0, -0.1, 4}, {0, -0.1, 8}, {0, -0.1, 16},
    //
    {0.4, 0.4, 1}, {0.4, 0.4, 2}, {0.4, 0.4, 4}, {0.4, 0.4, 8}, {0.4, 0.4, 16},
    {-0.4, -0.4, 1}, {-0.4, -0.4, 2}, {-0.4, -0.4, 4}, {-0.4, -0.4, 8}, {-0.4, -0.4, 16},
    //
    {0.9, 0, 1}, {0.9, 0, 2}, {0.9, 0, 4}, {0.9, 0, 8}, {0.9, 0, 16},
    {0, 0.9, 1}, {0, 0.9, 2}, {0, 0.9, 4}, {0, 0.9, 8}, {0, 0.9, 16},
    {-0.9, 0, 1}, {-0.9, 0, 2}, {-0.9, 0, 4}, {-0.9, 0, 8}, {-0.9, 0, 16},
    {0, -0.9, 1}, {0, -0.9, 2}, {0, -0.9, 4}, {0, -0.9, 8}, {0, -0.9, 16},
};

} // namespace

BOOST_DATA_TEST_CASE(test_ttt,  
                     bdata::xrange(deltaList.size()), deltaIdx)
{
  const Delta& delta = deltaList[deltaIdx];

  BOOST_TEST_MESSAGE(
   stringify("--------------------------------------------DELTA lat %.1f lon %.1f depth %.1f", 
            delta.lat, delta.lon, delta.depth));

  vector<double> travelTimeP(tttList.size(), 0);
  vector<double> takeOfAngleAzimP(tttList.size(), 0);
  vector<double> takeOfAngleDipP(tttList.size(), 0);
  vector<double> velocityAtSrcP(tttList.size(), 0);

  vector<double> travelTimeS(tttList.size(), 0);
  vector<double> takeOfAngleAzimS(tttList.size(), 0);
  vector<double> takeOfAngleDipS(tttList.size(), 0);
  vector<double> velocityAtSrcS(tttList.size(), 0);

  for (auto station : stationList)
  {
    double stationDepth = -(station.elevation / 1000.);
    double lat   = station.latitude + delta.lat;
    double lon   = station.longitude + delta.lon;
    double depth = stationDepth + delta.depth;

    for (size_t i = 0; i < tttList.size(); i++)
    {
      HDD::TravelTimeTablePtr ttt = tttList[i];

      BOOST_TEST_MESSAGE(
          stringify("Station %s lat %.1f lon %.1f depth %.1f -----------  TTT %-9s",
                    station.id.c_str(), station.latitude,
                    station.longitude, stationDepth, ttt->type.c_str()));

      BOOST_REQUIRE(ttt);
      BOOST_CHECK_NO_THROW(ttt->compute(lat, lon, depth, station, "P",
                                        travelTimeP[i], takeOfAngleAzimP[i],
                                        takeOfAngleDipP[i], velocityAtSrcP[i]));
      BOOST_CHECK_NO_THROW(ttt->compute(lat, lon, depth, station, "S",
                                        travelTimeS[i], takeOfAngleAzimS[i],
                                        takeOfAngleDipS[i], velocityAtSrcS[i]));
    }
  }

  for (size_t i = 0; i< tttList.size(); i++)
  {
    const HDD::TravelTimeTableCPtr& ttt = tttList[i];
    BOOST_TEST_MESSAGE(
        stringify("  TTT %-9s Travel time/Take-Off Angle Azimuth/Take-Off "
                  "Angle Dip/Velocity "
                  "P %5.2f %5.1f %5.1f %5.2f S %5.2f %5.1f %5.1f %5.2f",
                  ttt->type.c_str(), travelTimeP[i], takeOfAngleAzimP[i],
                  takeOfAngleDipP[i], velocityAtSrcP[i], travelTimeS[i],
                  takeOfAngleAzimS[i], takeOfAngleDipS[i], velocityAtSrcS[i]));
  }
}


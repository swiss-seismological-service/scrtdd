#define BOOST_TEST_MODULE libhdd
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>

#include "hdd/catalog.h"
#include "hdd/ttt.h"
#include "hdd/utils.h"
#include "hdd/test/common.ipp"
#include "ttt.h"

using namespace std;
using namespace HDD;
namespace bdata = boost::unit_test::data;

namespace {

const std::vector<TTTParams> tttListSC = {
    {"LOCSAT", {"iasp91"}},
    {"libtau", {"iasp91"}},
    {"homogeneous", {"iasp91"}},
    {"tttnll", {"iasp91_2D_simple"}},
    {"tttnll", {"iasp91_2D_sdc"}},
    {"tttnll", {"iasp91_2D_global"}},
    {"tttnll", {"iasp91_2D_azimuthal_equidist"}},
    {"tttnll", {"iasp91_2D_merc"}},
    {"tttnll", {"iasp91_2D_lambert"}},
    {"tttnll", {"iasp91_3D_simple"}},
    {"tttnll", {"iasp91_3D_sdc"}},
    {"tttnll", {"iasp91_3D_azimuthal_equidist"}},
    {"tttnll", {"iasp91_3D_merc"}},
    {"tttnll", {"iasp91_3D_lambert"}},
};


struct Delta
{
  double lat;
  double lon;
  double depth;
};

vector<Delta> deltaList = {
    //
    {   0,    0, -1},
    {   0,    0,  1},
    {0.01,    0,  0},
    {   0, 0.01,  0},
    //
    {0.05, 0, -1},
    {0.05, 0, 1},
    {0.05, 0, 2},
    {0.05, 0, 4},
    {0.05, 0, 8},
    {0.05, 0, 16},
    {0, 0.1, -1},
    {0, 0.1, 1},
    {0, 0.1, 2},
    {0, 0.1, 4},
    {0, 0.1, 8},
    {0, 0.1, 16},
    {-0.05, 0, -1},
    {-0.05, 0, 1},
    {-0.05, 0, 2},
    {-0.05, 0, 4},
    {-0.05, 0, 8},
    {-0.05, 0, 16},
    {0, -0.1, -1},
    {0, -0.1, 1},
    {0, -0.1, 2},
    {0, -0.1, 4},
    {0, -0.1, 8},
    {0, -0.1, 16},
    {0.05, 0.1, -1},
    {0.05, 0.1, 1},
    {0.05, 0.1, 2},
    {0.05, 0.1, 4},
    {0.05, 0.1, 8},
    {0.05, 0.1, 16},
    {-0.05, -0.1, -1},
    {-0.05, -0.1, 1},
    {-0.05, -0.1, 2},
    {-0.05, -0.1, 4},
    {-0.05, -0.1, 8},
    {-0.05, -0.1, 16},
};

vector<Delta> deltaListZeroElevation = {
    {0.05, 0, 4},
    {0.05, 0, 8},
    {0.05, 0, 16},
    {0.05, 0, 25},
    {0.05, 0, 35},
    {0, 0.05, 4},
    {0, 0.05, 8},
    {0, 0.05, 16},
    {0, 0.05, 25},
    {0, 0.05, 35},
    {-0.05, 0, 4},
    {-0.05, 0, 8},
    {-0.05, 0, 16},
    {-0.05, 0, 25},
    {-0.05, 0, 35},
    {0, -0.05, 4},
    {0, -0.05, 8},
    {0, -0.05, 16},
    {0, -0.05, 25},
    {0, -0.05, 35},
    //
    {0.1, 0, 4},
    {0.1, 0, 8},
    {0.1, 0, 16},
    {0.1, 0, 25},
    {0.1, 0, 35},
    {0, 0.1, 4},
    {0, 0.1, 8},
    {0, 0.1, 16},
    {0, 0.1, 25},
    {0, 0.1, 35},
    {-0.1, 0, 4},
    {-0.1, 0, 8},
    {-0.1, 0, 16},
    {-0.1, 0, 25},
    {-0.1, 0, 35},
    {0, -0.1, 4},
    {0, -0.1, 8},
    {0, -0.1, 16},
    {0, -0.1, 25},
    {0, -0.1, 35},
    //
    {0.4, 0, 4},
    {0.4, 0, 8},
    {0.4, 0, 16},
    {0.4, 0, 25},
    {0.4, 0, 35},
    {0, 0.4, 4},
    {0, 0.4, 8},
    {0, 0.4, 16},
    {0, 0.4, 25},
    {0, 0.4, 35},
    {-0.4, 0, 4},
    {-0.4, 0, 8},
    {-0.4, 0, 16},
    {-0.4, 0, 25},
    {-0.4, 0, 35},
    {0, -0.4, 4},
    {0, -0.4, 8},
    {0, -0.4, 16},
    {0, -0.4, 25},
    {0, -0.4, 35},
};

void test_ttt_internal(const std::vector<TTTParams>& tttPrms,
                       const std::vector<HDD::Catalog::Station>& stations,
                       const Delta &delta) {

 std::vector<shared_ptr<HDD::TravelTimeTable>> ttts;
 for (TTTParams prm : tttPrms)
  {
    shared_ptr<HDD::TravelTimeTable> ttt =
      prm.type.rfind("Native:", 0) == 0 
      ? createTTT(prm)
      : std::unique_ptr<HDD::TravelTimeTable>(
          new HDD::SCAdapter::TravelTimeTable(prm.type, prm.args[0]));
    BOOST_REQUIRE(ttt);
    ttts.push_back(ttt);
  }

  for (auto station : stations)
  {
    const double stationDepth = -(station.elevation / 1000.);
    const double lat          = station.latitude + delta.lat;
    const double lon          = station.longitude + delta.lon;
    const double depth        = stationDepth + delta.depth;

    BOOST_TEST_MESSAGE(strf("Testing station %s lat %.1f lon %.1f depth %.1f",
                            station.id.c_str(), station.latitude,
                            station.longitude, stationDepth));

    vector<double> travelTimeP(ttts.size(), 0);
    vector<double> takeOffAngleAzimP(ttts.size(), 0);
    vector<double> takeOffAngleDipP(ttts.size(), 0);
    vector<double> velocityAtSrcP(ttts.size(), 0);

    vector<double> travelTimeS(ttts.size(), 0);
    vector<double> takeOffAngleAzimS(ttts.size(), 0);
    vector<double> takeOffAngleDipS(ttts.size(), 0);
    vector<double> velocityAtSrcS(ttts.size(), 0);

    for (size_t i = 0; i < ttts.size(); i++)
    {
      BOOST_CHECK_NO_THROW(ttts[i]->compute(
          lat, lon, depth, station, "P", travelTimeP[i], takeOffAngleAzimP[i],
          takeOffAngleDipP[i], velocityAtSrcP[i]));
      BOOST_CHECK_NO_THROW(ttts[i]->compute(
          lat, lon, depth, station, "S", travelTimeS[i], takeOffAngleAzimS[i],
          takeOffAngleDipS[i], velocityAtSrcS[i]));

      BOOST_CHECK_EQUAL(travelTimeP[i],
                        ttts[i]->compute(lat, lon, depth, station, "P"));
      BOOST_CHECK_EQUAL(travelTimeS[i],
                        ttts[i]->compute(lat, lon, depth, station, "S"));
    }

    for (size_t i = 0; i < ttts.size(); i++)
    {
      BOOST_TEST_MESSAGE(strf(
          "TTT type %-9s [Travel time, Velocity, Take-Off Angle Azimuth and "
          "Dip] P [%5.2f %5.2f %4.f %4.f] S [%5.2f %5.2f %4.f %4.f]",
          tttPrms[i].type.c_str(), travelTimeP[i], velocityAtSrcP[i],
          radToDeg(takeOffAngleAzimP[i]), radToDeg(takeOffAngleDipP[i]),
          travelTimeS[i], velocityAtSrcS[i], radToDeg(takeOffAngleAzimS[i]),
          radToDeg(takeOffAngleDipS[i])));
    }

    for (size_t i = 0; i < ttts.size() - 1; i++)
    {
      for (size_t j = i + 1; j < ttts.size(); j++)
      {
        BOOST_CHECK_CLOSE(travelTimeP[i], travelTimeP[j], 12);
        BOOST_CHECK_CLOSE(travelTimeS[i], travelTimeS[j], 12);
        BOOST_CHECK_CLOSE(velocityAtSrcP[i], velocityAtSrcP[j], 5);
        BOOST_CHECK_CLOSE(velocityAtSrcS[i], velocityAtSrcS[j], 5);

        auto checkCloseAngles = [](double x, double y, double degreeTol) {
          double diffAngle = atan2(sin(x - y), cos(x - y));
          diffAngle        = radToDeg(diffAngle);
          BOOST_TEST_MESSAGE(" x " << radToDeg(x) << " y " << radToDeg(y)
                                   << " delta " << diffAngle);
          BOOST_CHECK_SMALL(diffAngle, degreeTol);
        };

        checkCloseAngles(takeOffAngleAzimP[i], takeOffAngleAzimP[j], 3);
        checkCloseAngles(takeOffAngleDipP[i], takeOffAngleDipP[j], 3);
        checkCloseAngles(takeOffAngleAzimS[i], takeOffAngleAzimS[j], 3);
        checkCloseAngles(takeOffAngleDipS[i], takeOffAngleDipS[j], 3);
      }
    }
  }
}

} // namespace

BOOST_DATA_TEST_CASE(test_ttt, bdata::xrange(deltaList.size()), deltaIdx)
{
  const Delta &delta = deltaList[deltaIdx];

  BOOST_TEST_MESSAGE(strf("Testing DELTA lat %.1f lon %.1f depth %.1f",
                          delta.lat, delta.lon, delta.depth));

  std::vector<TTTParams> tttPrms(tttList);
  for (const TTTParams& prms : tttListSC){
    if (prms.type == "homogeneous" || prms.type == "tttnll")
    {
      tttPrms.push_back(prms);
    }
  }
  test_ttt_internal(tttPrms, stationList, delta);
}

BOOST_DATA_TEST_CASE(test_ttt_zero_elevation, bdata::xrange(deltaListZeroElevation.size()), deltaIdx)
{
  const Delta &delta = deltaList[deltaIdx];

  BOOST_TEST_MESSAGE(strf("Testing DELTA lat %.1f lon %.1f depth %.1f",
                          delta.lat, delta.lon, delta.depth));

  std::vector<TTTParams> tttPrms;
  for (const TTTParams& prms : tttListSC){
    if (prms.type == "libtau" || prms.type == "LOCSAT" || prms.type == "tttnll")
    {
      tttPrms.push_back(prms);
    }
  }
  test_ttt_internal(tttPrms, stationListZeroElevation, delta);
}


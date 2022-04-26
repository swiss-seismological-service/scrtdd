#define SEISCOMP_TEST_MODULE hdd
#include <seiscomp/unittest/unittests.h>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>

#include "hdd/catalog.h"
#include "hdd/ttt.h"
#include "hdd/utils.h"
#include "common.ipp"

using namespace std;
using namespace HDD;
namespace bdata = boost::unit_test::data;

namespace {

struct Delta
{
  double lat;
  double lon;
  double depth;
};

vector<Delta> deltaList = {
    //
    {0.05, 0, 1},
    {0.05, 0, 2},
    {0.05, 0, 4},
    {0.05, 0, 8},
    {0.05, 0, 16},
    {0, 0.1, 1},
    {0, 0.1, 2},
    {0, 0.1, 4},
    {0, 0.1, 8},
    {0, 0.1, 16},
    {-0.05, 0, 1},
    {-0.05, 0, 2},
    {-0.05, 0, 4},
    {-0.05, 0, 8},
    {-0.05, 0, 16},
    {0, -0.1, 1},
    {0, -0.1, 2},
    {0, -0.1, 4},
    {0, -0.1, 8},
    {0, -0.1, 16},
    {0.05, 0.1, 1},
    {0.05, 0.1, 2},
    {0.05, 0.1, 4},
    {0.05, 0.1, 8},
    {0.05, 0.1, 16},
    {-0.05, -0.1, 1},
    {-0.05, -0.1, 2},
    {-0.05, -0.1, 4},
    {-0.05, -0.1, 8},
    {-0.05, -0.1, 16},
    //
    {0.1, 0, 1},
    {0.1, 0, 2},
    {0.1, 0, 4},
    {0.1, 0, 8},
    {0.1, 0, 16},
    {0, 0.1, 1},
    {0, 0.1, 2},
    {0, 0.1, 4},
    {0, 0.1, 8},
    {0, 0.1, 16},
    {-0.1, 0, 1},
    {-0.1, 0, 2},
    {-0.1, 0, 4},
    {-0.1, 0, 8},
    {-0.1, 0, 16},
    {0, -0.1, 1},
    {0, -0.1, 2},
    {0, -0.1, 4},
    {0, -0.1, 8},
    {0, -0.1, 16},
    //
    {0.4, 0.4, 1},
    {0.4, 0.4, 2},
    {0.4, 0.4, 4},
    {0.4, 0.4, 8},
    {0.4, 0.4, 16},
    {-0.4, -0.4, 1},
    {-0.4, -0.4, 2},
    {-0.4, -0.4, 4},
    {-0.4, -0.4, 8},
    {-0.4, -0.4, 16},
    //
    {0.7, 0, 1},
    {0.7, 0, 2},
    {0.7, 0, 4},
    {0.7, 0, 8},
    {0.7, 0, 16},
    {0, 0.7, 1},
    {0, 0.7, 2},
    {0, 0.7, 4},
    {0, 0.7, 8},
    {0, 0.7, 16},
    {-0.7, 0, 1},
    {-0.7, 0, 2},
    {-0.7, 0, 4},
    {-0.7, 0, 8},
    {-0.7, 0, 16},
    {0, -0.7, 1},
    {0, -0.7, 2},
    {0, -0.7, 4},
    {0, -0.7, 8},
    {0, -0.7, 16},
};

} // namespace

BOOST_DATA_TEST_CASE(test_ttt, bdata::xrange(deltaList.size()), deltaIdx)
{
  const Delta &delta = deltaList[deltaIdx];

  BOOST_TEST_MESSAGE(strf("Testing DELTA lat %.1f lon %.1f depth %.1f",
                          delta.lat, delta.lon, delta.depth));

  for (auto station : stationList)
  {
    const double stationDepth = -(station.elevation / 1000.);
    const double lat          = station.latitude + delta.lat;
    const double lon          = station.longitude + delta.lon;
    const double depth        = stationDepth + delta.depth;

    BOOST_TEST_MESSAGE(strf("Testing station %s lat %.1f lon %.1f depth %.1f",
                            station.id.c_str(), station.latitude,
                            station.longitude, stationDepth));

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
      unique_ptr<HDD::TravelTimeTable> ttt = createTTT(tttList[i]);
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
      BOOST_TEST_MESSAGE(strf(
          "TTT type %-9s [Travel time, Velocity, Take-Off Angle Azimuth and "
          "Dip] P [%5.2f %5.2f %4.f %4.f] S [%5.2f %5.2f %4.f %4.f]",
          tttList[i].type.c_str(), travelTimeP[i], velocityAtSrcP[i],
          radToDeg(takeOffAngleAzimP[i]), radToDeg(takeOffAngleDipP[i]),
          travelTimeS[i], velocityAtSrcS[i], radToDeg(takeOffAngleAzimS[i]),
          radToDeg(takeOffAngleDipS[i])));
    }

    for (size_t i = 0; i < tttList.size() - 1; i++)
    {
      for (size_t j = i + 1; j < tttList.size(); j++)
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

        checkCloseAngles(takeOffAngleAzimP[i], takeOffAngleAzimP[j], 15);
        checkCloseAngles(takeOffAngleDipP[i], takeOffAngleDipP[j], 15);
        checkCloseAngles(takeOffAngleAzimS[i], takeOffAngleAzimS[j], 15);
        checkCloseAngles(takeOffAngleDipS[i], takeOffAngleDipS[j], 15);
      }
    }
  }
}

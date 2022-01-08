#define SEISCOMP_TEST_MODULE hdd
#include <seiscomp/unittest/unittests.h>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>

#include "catalog.h"
#include "ellipsoid.h"
#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/math.h>

using namespace std;
using namespace Seiscomp;
using Seiscomp::Core::stringify;
using Event     = HDD::Catalog::Event;
using Phase     = HDD::Catalog::Phase;
using Station   = HDD::Catalog::Station;
namespace bdata = boost::unit_test::data;

namespace {

struct Point
{
  double latitude;
  double longitude;
  double depth;
};

vector<Point> orgList = {
    {0, 0, 10},        {15, 90, 7},      {15, -90, 7},      {-60, 120, 15},
    {-60, -120, 15},   {45, 45, 3},      {-45, -45, 2},     {-45, 45, 3},
    {45, -45, 2},      {45, 180, 10},    {-45, 180, 10},    {60, -180, 5},
    {-60, -180, 5},    {30, 179.99, 10}, {-30, 179.99, 10}, {-25, -179.99, 10},
    {25, -179.99, 10}, {0, 178, 4},      {0, -178, 1},      {0, -179, 8},
    {0, 179, 11}};

} // namespace

BOOST_DATA_TEST_CASE(test_ellipsoid, bdata::xrange(orgList.size()), orgIdx)
{
  const Point &org = orgList[orgIdx];

  for (double distance : {0.5, 1.5, 13.5, 25.})
  {
    BOOST_TEST_MESSAGE(
        stringify("Ellipsoid lat %.3f lon %.3f depth %.3f - Distance %f",
                  org.latitude, org.longitude, org.depth, distance));

    HDD::HddEllipsoid ellip(distance * 2, org.latitude, org.longitude,
                            org.depth);

    double pointLat, pointLon, pointDepth;

    // Inside ellipsoid
    distance = Math::Geo::km2deg(distance - 0.01);

    // upper quadrants
    pointDepth = org.depth + 0.001;

    Math::Geo::delandaz2coord(distance, 45, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 1));
    Math::Geo::delandaz2coord(distance, 315, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 2));
    Math::Geo::delandaz2coord(distance, 225, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 3));
    Math::Geo::delandaz2coord(distance, 135, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 4));

    // bottom quadrants
    pointDepth = org.depth - 0.001;

    Math::Geo::delandaz2coord(distance, 45, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 5));
    Math::Geo::delandaz2coord(distance, 315, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 6));
    Math::Geo::delandaz2coord(distance, 225, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 7));
    Math::Geo::delandaz2coord(distance, 135, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 8));

    pointDepth = org.depth;

    // inside the inner ellipsoid
    distance = Math::Geo::km2deg(distance / 2 - 0.01);

    Math::Geo::delandaz2coord(distance, 45, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
    Math::Geo::delandaz2coord(distance, 315, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
    Math::Geo::delandaz2coord(distance, 225, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
    Math::Geo::delandaz2coord(distance, 135, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));

    // outside the outer ellipsoid
    distance = Math::Geo::km2deg(distance + 0.01);

    Math::Geo::delandaz2coord(distance, 45, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
    Math::Geo::delandaz2coord(distance, 315, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
    Math::Geo::delandaz2coord(distance, 225, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
    Math::Geo::delandaz2coord(distance, 135, org.latitude, org.longitude,
                              &pointLat, &pointLon);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
  }
}

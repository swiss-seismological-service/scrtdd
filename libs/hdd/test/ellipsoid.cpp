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

#include "hdd/catalog.h"
#include "hdd/ellipsoid.h"
#include "hdd/utils.h"

using namespace std;
using namespace HDD;
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
        strf("Ellipsoid lat %.3f lon %.3f depth %.3f - Distance %f",
             org.latitude, org.longitude, org.depth, distance));

    HDD::HddEllipsoid ellip(distance * 2, org.latitude, org.longitude,
                            org.depth);

    double pointLat, pointLon, pointDepth;

    // Inside ellipsoid
    distance -= 0.01;

    // upper quadrants
    pointDepth = org.depth + 0.001;

    computeCoordinates(distance, degToRad(45), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 1));
    computeCoordinates(distance, degToRad(315), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 2));
    computeCoordinates(distance, degToRad(225), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 3));
    computeCoordinates(distance, degToRad(135), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 4));

    // bottom quadrants
    pointDepth = org.depth - 0.001;

    computeCoordinates(distance, degToRad(45), org.latitude, org.longitude, pointLat,
                       pointLon);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 5));
    computeCoordinates(distance, degToRad(315), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 6));
    computeCoordinates(distance, degToRad(225), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 7));
    computeCoordinates(distance, degToRad(135), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(ellip.isInside(pointLat, pointLon, pointDepth, 8));

    pointDepth = org.depth;

    // inside the inner ellipsoid
    distance = distance / 2 - 0.01;

    computeCoordinates(distance, degToRad(45), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
    computeCoordinates(distance, degToRad(315), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
    computeCoordinates(distance, degToRad(225), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
    computeCoordinates(distance, degToRad(135), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));

    // outside the outer ellipsoid
    distance = distance + 0.01;

    computeCoordinates(distance, degToRad(45), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
    computeCoordinates(distance, degToRad(315), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
    computeCoordinates(distance, degToRad(225), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
    computeCoordinates(distance, degToRad(135), org.latitude, org.longitude, pointLat,
                       pointLon, org.depth);
    BOOST_CHECK(!ellip.isInside(pointLat, pointLon, pointDepth));
  }
}

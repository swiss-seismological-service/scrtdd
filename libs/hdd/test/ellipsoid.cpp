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

#include <map>

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

//
// Verify the ellipsoidal layers used by the Waldhauser (2009) neighbour
// selection. HddEllipsoid(2 * d, ...) is the volume inside an outer ellipsoid
// (horizontal semi-axes d, vertical semi-axis 2d) and outside an inner
// ellipsoid (half those semi-axes), divided in 8 quadrants: 1-4 are deeper
// than the center, 5-8 shallower. Points are tested close to the layer
// boundaries, along the horizontal diagonal of each quadrant and along the
// vertical axis.
//
BOOST_DATA_TEST_CASE(test_ellipsoid, bdata::xrange(orgList.size()), orgIdx)
{
  const Point &org = orgList[orgIdx];

  // azimuth of the horizontal diagonal of each quadrant
  const map<int, double> quadrantAzimuth = {{1, 45},  {2, 315}, {3, 225},
                                            {4, 135}, {5, 45},  {6, 315},
                                            {7, 225}, {8, 135}};

  for (const double d : {0.5, 1.5, 13.5, 25.})
  {
    BOOST_TEST_MESSAGE(
        strf("Ellipsoid lat %.3f lon %.3f depth %.3f - Distance %f",
             org.latitude, org.longitude, org.depth, d));

    const HDD::HddEllipsoid ellip(d * 2, org.latitude, org.longitude,
                                  org.depth);

    for (const auto &kv : quadrantAzimuth)
    {
      const int quadrant   = kv.first;
      const double azimuth = degToRad(kv.second);
      const double depth = quadrant <= 4 ? org.depth + 0.001 : org.depth - 0.001;
      double lat, lon;

      // inside the layer: inside the outer ellipsoid and outside the inner
      // one, in this quadrant only
      computeCoordinates(d - 0.01, azimuth, org.latitude, org.longitude, lat,
                         lon, org.depth);
      BOOST_CHECK(ellip.isInside(lat, lon, depth));
      for (int q = 1; q <= 8; q++)
      {
        BOOST_CHECK_EQUAL(ellip.isInside(lat, lon, depth, q), q == quadrant);
      }

      // inside the inner ellipsoid
      computeCoordinates(d / 2 - 0.01, azimuth, org.latitude, org.longitude,
                         lat, lon, org.depth);
      BOOST_CHECK(!ellip.isInside(lat, lon, depth));

      // outside the outer ellipsoid
      computeCoordinates(d + 0.01, azimuth, org.latitude, org.longitude, lat,
                         lon, org.depth);
      BOOST_CHECK(!ellip.isInside(lat, lon, depth));
    }

    // along the vertical axis: inside the layer, inside the inner ellipsoid
    // and outside the outer ellipsoid, below and above the center
    for (const double sign : {1., -1.})
    {
      BOOST_CHECK(
          ellip.isInside(org.latitude, org.longitude, org.depth + sign * 1.5 * d));
      BOOST_CHECK(
          !ellip.isInside(org.latitude, org.longitude, org.depth + sign * 0.5 * d));
      BOOST_CHECK(!ellip.isInside(org.latitude, org.longitude,
                                  org.depth + sign * (2 * d + 0.01)));
    }
  }
}

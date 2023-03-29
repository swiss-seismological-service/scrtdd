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

#ifndef __HDD_ELLIPSOID_H__
#define __HDD_ELLIPSOID_H__

#include "utils.h"
#include <set>

namespace HDD {

/*
 * Ellipsoid standard equation:
 *
 *      (x-xo)^2 / axix_a^2 + (y-yo)^2 / axis_b^2 + (z-zo)^2 / axis_c^2 = 1
 *
 */
struct Ellipsoid
{
  bool isInside(double lat, double lon, double depth) const
  {
    double az;
    double distance = computeDistance(lat, lon, this->lat, this->lon, &az,
                                      nullptr, this->depth);

    double dist_x = distance * std::sin(az);
    double dist_y = distance * std::cos(az);
    double dist_z = depth - this->depth;

    double one = square(dist_x / axis_a) + square(dist_y / axis_b) +
                 square(dist_z / axis_c);
    return one <= 1;
  }

  double axis_a = 0, axis_b = 0, axis_c = 0; // axis in km
  double lat = 0, lon = 0, depth = 0;        // origin
};

/*
 * Helper class to implement Waldhauser's paper method of neighboring events
 * selection based on 5 concentric ellipsoidal layers.
 *
 * Quadrants (1-4 above depth, 5-8 below depth):
 *
 *        lat
 *         ^
 *         |
 *    2/6  |   1/5
 *         |
 * -----------------> lon
 *         |
 *    3/7  |   4/8
 *         |
 */
class HddEllipsoid
{

public:
  HddEllipsoid(const HddEllipsoid &other)            = default;
  HddEllipsoid &operator=(const HddEllipsoid &other) = default;

  HddEllipsoid(HddEllipsoid &&other)            = default;
  HddEllipsoid &operator=(HddEllipsoid &&other) = default;

  ~HddEllipsoid() = default;

  HddEllipsoid(double vertical_axis_len, double lat, double lon, double depth)
      : HddEllipsoid(vertical_axis_len / 2., vertical_axis_len, lat, lon, depth)
  {}

  HddEllipsoid(double inner_vertical_axis_len,
               double outer_vertical_axis_len,
               double lat,
               double lon,
               double depth)
  {
    _outerEllipsoid.axis_a = outer_vertical_axis_len / 2.;
    _outerEllipsoid.axis_b = _outerEllipsoid.axis_a;
    _outerEllipsoid.axis_c = outer_vertical_axis_len;

    _innerEllipsoid.axis_a = inner_vertical_axis_len / 2.;
    _innerEllipsoid.axis_b = _innerEllipsoid.axis_a;
    _innerEllipsoid.axis_c = inner_vertical_axis_len;

    _outerEllipsoid.lat = _innerEllipsoid.lat = lat;
    _outerEllipsoid.lon = _innerEllipsoid.lon = lon;
    _outerEllipsoid.depth = _innerEllipsoid.depth = depth;
  }

  const Ellipsoid &getInnerEllipsoid() const { return _innerEllipsoid; }
  const Ellipsoid &getOuterEllipsoid() const { return _outerEllipsoid; }

  // Returns if the coordinate provided is located within the correct quadrant,
  // and is both inside the outer layer and outside of the inner one.
  bool
  isInside(double lat, double lon, double depth, int quadrant /* 1-8 */) const
  {
    return isInQuadrant(_innerEllipsoid, lat, lon, depth, quadrant) &&
           isInside(lat, lon, depth);
  }

  // Returns if the coordinate provided is located both inside the outer layer
  // and outside of the inner one.
  bool isInside(double lat, double lon, double depth) const
  {
    return _outerEllipsoid.isInside(lat, lon, depth) &&
           !_innerEllipsoid.isInside(lat, lon, depth);
  }

  static bool isInQuadrant(const Ellipsoid &ellip,
                           double lat,
                           double lon,
                           double depth,
                           int quadrant /* 1-8 */)
  {
    if (quadrant < 1 || quadrant > 8)
      throw Exception("quadrant should be between 1 and 8");

    if (depth < ellip.depth && std::set<int>{1, 2, 3, 4}.count(quadrant) != 0)
      return false;
    if (depth > ellip.depth && std::set<int>{5, 6, 7, 8}.count(quadrant) != 0)
      return false;

    if (lat < ellip.lat && std::set<int>{1, 2, 5, 6}.count(quadrant) != 0)
      return false;
    if (lat > ellip.lat && std::set<int>{3, 4, 7, 8}.count(quadrant) != 0)
      return false;

    double lonDelta = normalizeLon(lon - ellip.lon);

    if (lonDelta < 0 && std::set<int>{1, 4, 5, 8}.count(quadrant) != 0)
      return false;
    if (lonDelta > 0 && std::set<int>{2, 3, 6, 7}.count(quadrant) != 0)
      return false;

    return true;
  }

private:
  Ellipsoid _outerEllipsoid;
  Ellipsoid _innerEllipsoid;
};

} // namespace HDD

#endif

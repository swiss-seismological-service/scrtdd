/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU Affero General Public License as published*
 * by the Free Software Foundation, either version 3 of the License, or    *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This program is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#ifndef __RTDD_APPLICATIONS_ELLIPSOID_H__
#define __RTDD_APPLICATIONS_ELLIPSOID_H__

#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/math.h>
#include <set>

namespace Seiscomp {
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
        double distance, az, baz;
        Math::Geo::delazi(lat, lon, this->lat, this->lon, &distance, &az, &baz);

        distance = Math::Geo::deg2km(distance); // distance to km
        az += orientation; // correct azimuth by the orientation of a and b axes
        az = deg2rad(az);

        double dist_x = distance * std::cos(az);
        double dist_y = distance * std::sin(az);
        double dist_z = depth - this->depth;

        double one = std::pow( dist_x / axis_a, 2) +
                     std::pow( dist_y / axis_b, 2) +
                     std::pow( dist_z / axis_c, 2);
        return one <= 1;
    }

    double axis_a=0, axis_b=0, axis_c=0; // axis in km
    double lat=0, lon=0, depth=0;        // origin
    double orientation = 0;              // degree: when 0 -> axis_a is East-West and axis_b is North-South
};

DEFINE_SMARTPOINTER(Ellipsoid);

/*
 * Helper class to implement Waldhauser's paper method of neighboring events selection
 * based on 5 concentric ellipsoidal layers
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
class HddEllipsoid : public Core::BaseObject
{

public:

    HddEllipsoid(double vertical_axis_len, double lat, double lon, double depth)
        : HddEllipsoid(vertical_axis_len/2., vertical_axis_len, lat, lon, depth)
    { }

    HddEllipsoid(double inner_vertical_axis_len, double outer_vertical_axis_len,
                 double lat, double lon, double depth)
    {
        _outerEllipsoid.axis_a = outer_vertical_axis_len / 2.;
        _outerEllipsoid.axis_b = _outerEllipsoid.axis_a;
        _outerEllipsoid.axis_c = outer_vertical_axis_len;

        _innerEllipsoid.axis_a = inner_vertical_axis_len / 2.;
        _innerEllipsoid.axis_b = _innerEllipsoid.axis_a;
        _innerEllipsoid.axis_c = inner_vertical_axis_len;

        _outerEllipsoid.lat   = _innerEllipsoid.lat   = lat;
        _outerEllipsoid.lon   = _innerEllipsoid.lon   = lon;
        _outerEllipsoid.depth = _innerEllipsoid.depth = depth;
    }

    const Ellipsoid& getInnerEllipsoid() const { return _innerEllipsoid; }
    const Ellipsoid& getOuterEllipsoid() const { return _outerEllipsoid; }

    bool isInside(double lat, double lon, double depth, int quadrant/* 1-8 */) const
    {
        // be in the right quadrant and inside the outer layer and outside of the inner one
        return isInQuadrant(_innerEllipsoid, lat, lon, depth, quadrant) &&  isInside(lat, lon, depth);
    }

    bool isInside(double lat, double lon, double depth) const
    {
        // be inside the outer layer and outside of the inner one
        return _outerEllipsoid.isInside(lat, lon, depth) && ! _innerEllipsoid.isInside(lat, lon, depth);
    }

    static bool
    isInQuadrant(const Ellipsoid& ellip, double lat, double lon, double depth, int quadrant /* 1-8 */)
    {
        if (quadrant < 1 || quadrant > 8)
            throw std::invalid_argument( "quadrant should be between 1 and 8");

        if (depth < ellip.depth && std::set<int>{1,2,3,4}.count(quadrant) != 0) return false;
        if (depth > ellip.depth && std::set<int>{5,6,7,8}.count(quadrant) != 0) return false;

        if (lon < ellip.lon && std::set<int>{1,4,5,8}.count(quadrant) != 0) return false;
        if (lon > ellip.lon && std::set<int>{2,3,6,7}.count(quadrant) != 0) return false;

        if (lat < ellip.lat && std::set<int>{1,2,5,6}.count(quadrant) != 0) return false;
        if (lat > ellip.lat && std::set<int>{3,4,7,8}.count(quadrant) != 0) return false;

        return true;
    }

private:

    Ellipsoid _outerEllipsoid;
    Ellipsoid _innerEllipsoid;
};

DEFINE_SMARTPOINTER(HddEllipsoid);


}
}

#endif

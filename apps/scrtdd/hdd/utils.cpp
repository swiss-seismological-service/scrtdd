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

#include "utils.h"
#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/math.h>

using namespace std;

namespace Seiscomp {
namespace HDD {

/*
 * Compute distance in km between two points and optionally
 * azimuth and backazimuth
 */
double computeDistance(double lat1, double lon1, double depth1,
                       double lat2, double lon2, double depth2,
                       double *azimuth, double *backAzimuth)
{
    double Hdist, az, baz;
    Math::Geo::delazi(lat1, lon1, lat2, lon2, &Hdist, &az, &baz);
    Hdist = Math::Geo::deg2km(Hdist);

    if (azimuth) *azimuth = az;
    if (backAzimuth) *backAzimuth = baz;

    if ( depth1 == depth2 )
        return Hdist;

    // this is an approximation that works when the distance is small
    // and the Earth curvature can be assumed flat 
    double Vdist = abs(depth1 - depth2);
    return std::sqrt( std::pow(Hdist,2) + std::pow(Vdist,2) );
}

double computeDistance(const Catalog::Event& ev1, const Catalog::Event& ev2,
                       double *azimuth, double *backAzimuth)
{
    return computeDistance(ev1.latitude, ev1.longitude, ev1.depth,
                           ev2.latitude, ev2.longitude, ev2.depth,
                           azimuth, backAzimuth);
}

double computeDistance(const Catalog::Event& event, const Catalog::Station& station,
                       double *azimuth, double *backAzimuth)
{
    return computeDistance(event.latitude, event.longitude, event.depth,
                           station.latitude, station.longitude, -(station.elevation/1000.),
                           azimuth, backAzimuth);
}


double computeMedian(const std::vector<double>& values)
{
    if (values.size() == 0) return 0;

    vector<double> tmp(values);
    const auto middleItr = tmp.begin() + tmp.size() / 2;
    std::nth_element(tmp.begin(), middleItr, tmp.end());
    double median = *middleItr;
    if (tmp.size() % 2 == 0)
    {
        const auto leftMiddleItr = std::max_element(tmp.begin(), middleItr);
        median = (*leftMiddleItr + *middleItr) / 2;
    }
    return median;
}


double computeMedianAbsoluteDeviation(const std::vector<double>& values, const double median)
{
    vector<double> absoluteDeviations( values.size() );
    for ( unsigned i = 0; i < values.size(); i++ )
    {
        absoluteDeviations[i] = std::abs(values[i] - median);
    }
    return computeMedian(absoluteDeviations);
}


double computeMean(const vector<double>& values)
{
    if (values.size() == 0) return 0;
    return std::accumulate( values.begin(), values.end(), 0.0) / values.size();
}


double computeMeanAbsoluteDeviation(const std::vector<double>& values, const double mean)
{
    vector<double> absoluteDeviations( values.size() );
    for ( unsigned i = 0; i < values.size(); i++ )
    {
        absoluteDeviations[i] = std::abs(values[i] - mean);
    }
    return computeMean(absoluteDeviations);
}


}
}

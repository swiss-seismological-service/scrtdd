/***************************************************************************
 * Copyright (C) gempa GmbH                                                *
 * All rights reserved.                                                    *
 * Contact: gempa GmbH (seiscomp-dev@gempa.de)                             *
 *                                                                         *
 * GNU Affero General Public License Usage                                 *
 * This file may be used under the terms of the GNU Affero                 *
 * Public License version 3.0 as published by the Free Software Foundation *
 * and appearing in the file LICENSE included in the packaging of this     *
 * file. Please review the following information to ensure the GNU Affero  *
 * Public License version 3.0 requirements will be met:                    *
 * https://www.gnu.org/licenses/agpl-3.0.html.                             *
 *                                                                         *
 * Other Usage                                                             *
 * Alternatively, this file may be used in accordance with the terms and   *
 * conditions contained in a signed written agreement between you and      *
 * gempa GmbH.                                                             *
 ***************************************************************************/


#ifndef __HDD_GEO__
#define __HDD_GEO__

#include <vector>

namespace HDD
{

/**
 * For two points (lat1, lon1) and (lat2, lon2),
 * the angular distance 'dist' in degrees,
 * the azimuth 'azi1' (azimuth of point 2 seen from point 1) and
 * the azimuth 'azi2' (azimuth of point 1 seen from point 2)
 * are computed.
 */
void delazi(double lat1, double lon1, double lat2, double lon2,
            double *out_dist, double *out_azi1, double *out_azi2);


/**
 * For two points (lat1, lon1) and (lat2, lon2),
 * the angular distance 'dist' in degrees,
 * the azimuth 'azi1' (azimuth of point 2 seen from point 1) and
 * the azimuth 'azi2' (azimuth of point 1 seen from point 2)
 * are computed.
 */
void delazi_wgs84(double lat1, double lon1, double lat2, double lon2,
                  double *out_dist, double *out_azi1, double *out_azi2);


/**
 * Computes the coordinates (lat, lon) of the point which
 * is at an azimuth of 'azi' and a distance of 'dist' as seen
 * from the point (lat0, lon0).
 *      -> lat, lon
 */
void delandaz2coord(double dist, double azi, double lat0, double lon0,
                    double *out_lat, double *out_lon);

#define KM_OF_DEGREE 111.1329149013519096

template<typename T>
T deg2km(T deg) { return deg * (T)KM_OF_DEGREE; }

template<typename T>
T km2deg(T km)  { return km / (T)KM_OF_DEGREE; }

}

#endif

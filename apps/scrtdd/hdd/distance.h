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

#ifndef __RTDD_APPLICATIONS_DISTANCE_H__
#define __RTDD_APPLICATIONS_DISTANCE_H__

#include "catalog.h"
#include <seiscomp3/core/strings.h>
#include <unordered_map>

namespace Seiscomp {
namespace HDD { 

double computeDistance(double lat1, double lon1, double depth1,
                       double lat2, double lon2, double depth2,
                      double *azimuth = nullptr, double *backAzimuth = nullptr);

double computeDistance(const Catalog::Event& ev1, const Catalog::Event& ev2,
                       double *azimuth = nullptr, double *backAzimuth = nullptr);

double computeDistance(const Catalog::Event& event, const Catalog::Station& station,
                       double *azimuth = nullptr, double *backAzimuth = nullptr);

class DistanceCache
{

public:

    double compute(double lat1, double lon1, double depth1,
                   double lat2, double lon2, double depth2,
                   double *azimuth = nullptr, double *backAzimuth = nullptr)
    {
        double distance, _azimuth, _backAzimuth;

        std::string key;

        key = Seiscomp::Core::stringify("%.6f-%.6f-%.6f-%.6f-%.6f-%.6f",
                                        lat1, lon1, depth1, lat2, lon2, depth2);

        if ( _cache.find(key) == _cache.end() )
        {
            distance = computeDistance(lat1, lon1, depth1, lat2, lon2, depth2,
                                       &_azimuth, &_backAzimuth);
            _cache[key] = Entry( {distance, _azimuth, _backAzimuth} );
        }
        else
        {
            const Entry& entry = _cache.at(key);
            distance = entry.distance;
            _azimuth = entry.azimuth;
            _backAzimuth = entry.backAzimuth;
        }

        if (azimuth) *azimuth = _azimuth;
        if (backAzimuth) *backAzimuth = _backAzimuth;

        return distance;
    }

    void clear() { _cache.clear(); }

private:

    struct Entry {
        double distance, azimuth, backAzimuth;
    };
    std::unordered_map<std::string,Entry> _cache;

};

}
}

#endif 


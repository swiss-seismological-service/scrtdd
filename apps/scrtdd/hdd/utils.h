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
#include <vector>
#include <random>

namespace Seiscomp {
namespace HDD { 

double computeDistance(double lat1, double lon1, double depth1,
                       double lat2, double lon2, double depth2,
                      double *azimuth = nullptr, double *backAzimuth = nullptr);

double computeDistance(const Catalog::Event& ev1, const Catalog::Event& ev2,
                       double *azimuth = nullptr, double *backAzimuth = nullptr);

double computeDistance(const Catalog::Event& event, const Catalog::Station& station,
                       double *azimuth = nullptr, double *backAzimuth = nullptr);


double computeMedian(const std::vector<double>& values);

double computeMedianAbsoluteDeviation(const std::vector<double>& values, const double median);

double computeMean(const std::vector<double>& values);

double computeMeanAbsoluteDeviation(const std::vector<double>& values, const double mean);


class Randomer {

public:

    Randomer(size_t min, size_t max, unsigned int seed = std::random_device{}())
        : gen_{seed}, dist_{min, max}
    { }

    // if you want predictable numbers
    void setSeed(unsigned int seed)
    {
        gen_.seed(seed);
    }

    size_t next()
    {
        return dist_(gen_);
    }

private:

    // random seed by default
    std::mt19937 gen_;
    std::uniform_int_distribution<size_t> dist_;
}; 

}
}

#endif 


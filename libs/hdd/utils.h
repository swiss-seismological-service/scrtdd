/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as          *
 * published by the Free Software Foundation, either version 3 of the      *
 * License, or (at your option) any later version.                         *
 *                                                                         *
 * This software is distributed in the hope that it will be useful,        *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#ifndef __HDD_UTILS_H__
#define __HDD_UTILS_H__

#include "catalog.h"
#include "timewindow.h"
#include "xcorrcache.h"
#include <cmath>
#include <initializer_list>
#include <random>
#include <regex>
#include <stdexcept>
#include <vector>

namespace HDD {

class Exception : public std::runtime_error
{
public:
  Exception(const std::string &message) : std::runtime_error(message) {}
};

std::string strf(const char *fmt, ...) __attribute__((format(printf, 1, 2)));

std::vector<std::string> splitString(const std::string &str,
                                     const std::regex &regex);

inline double radToDeg(double r) { return 180.0 * r / M_PI; }

inline double degToRad(double d) { return M_PI * d / 180.0; }

template <typename T> T square(T x) { return x * x; }

static const double EARTH_MEAN_RADIUS_METER = 6371008.77141506;

inline double kmOfDegree(double kmDepth = 0)
{
  return ((EARTH_MEAN_RADIUS_METER / 1000.) - kmDepth) * M_PI / 180;
}

inline double deg2km(double deg, double kmDepth = 0)
{
  return deg * kmOfDegree(kmDepth);
}

inline double km2deg(double km, double kmDepth = 0)
{
  return km / kmOfDegree(kmDepth);
}

inline double rad2km(double rad, double kmDepth = 0)
{
  return rad * ((EARTH_MEAN_RADIUS_METER / 1000.) - kmDepth);
}

inline double km2rad(double km, double kmDepth = 0)
{
  return km / ((EARTH_MEAN_RADIUS_METER / 1000.) - kmDepth);
}

inline double normalizeLon(double lon)
{
  while (lon < -180.0) lon += 360.0;
  while (lon > 180.0) lon -= 360.0;
  return lon;
}

inline double normalizeAzimuth(double az)
{
  while (az < 0) az += 360.0;
  while (az > 360) az -= 360.0;
  return az;
}

void computeCoordinates(double distance,
                        double azimuth,
                        double clat,
                        double clon,
                        double &lat,
                        double &lon,
                        double atKmDepth = 0);

double computeDistance(double lat1,
                       double lon1,
                       double lat2,
                       double lon2,
                       double *azimuth     = nullptr,
                       double *backAzimuth = nullptr,
                       double atKmDepth    = 0);

double computeDistance(double lat1,
                       double lon1,
                       double kmDepth1,
                       double lat2,
                       double lon2,
                       double kmDepth2,
                       double *azimuth     = nullptr,
                       double *backAzimuth = nullptr);

double computeDistance(const Catalog::Event &ev1,
                       const Catalog::Event &ev2,
                       double *azimuth     = nullptr,
                       double *backAzimuth = nullptr);

double computeDistance(const Catalog::Event &event,
                       const Catalog::Station &station,
                       double *azimuth     = nullptr,
                       double *backAzimuth = nullptr);

double computeMedian(const std::vector<double> &values);

double computeMedianAbsoluteDeviation(const std::vector<double> &values,
                                      const double median);

double computeMean(const std::vector<double> &values);

double computeMeanAbsoluteDeviation(const std::vector<double> &values,
                                    const double mean);

double computeCircularMean(const std::vector<double> &angles);

std::string joinPath(const std::string &path1, const std::string &path2);

bool pathExists(const std::string &path);

bool directoryEmpty(const std::string &path);

bool createDirectories(const std::string &path);

bool removePath(const std::string &path);

void writeXCorrToFile(const XCorrCache &xcorr,
                      const Catalog &cat,
                      const std::string &file);

XCorrCache readXCorrFromFile(const Catalog &cat, const std::string &file);

} // namespace HDD

#endif

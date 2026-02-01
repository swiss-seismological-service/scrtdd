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

#ifndef __HDD_UTILS_H__
#define __HDD_UTILS_H__

#include "catalog.h"
#include "timewindow.h"
#include "xcorrcache.h"
#include <cmath>
#include <cstdarg>
#include <initializer_list>
#include <random>
#include <regex>
#include <stdexcept>
#include <vector>

namespace HDD {

class Exception : public std::runtime_error
{
public:
  explicit Exception(const std::string &message) : std::runtime_error(message)
  {}
};

std::string strf(const char *fmt, ...) __attribute__((format(printf, 1, 2)));

std::string strf(const char *fmt, va_list args);

std::vector<std::string> splitString(const std::string &str,
                                     const std::regex &regex);

inline constexpr double radToDeg(double r) { return 180.0 * r / M_PI; }

inline constexpr double degToRad(double d) { return M_PI * d / 180.0; }

template <typename T> inline constexpr T square(T x) { return x * x; }

constexpr double EARTH_MEAN_RADIUS_METER = 6371008.77141506;

inline constexpr double kmOfDegree(double kmDepth = 0)
{
  return ((EARTH_MEAN_RADIUS_METER / 1000.) - kmDepth) * M_PI / 180.0;
}

inline constexpr double deg2km(double deg, double kmDepth = 0)
{
  return deg * kmOfDegree(kmDepth);
}

inline constexpr double km2deg(double km, double kmDepth = 0)
{
  return km / kmOfDegree(kmDepth);
}

inline constexpr double rad2km(double rad, double kmDepth = 0)
{
  return rad * ((EARTH_MEAN_RADIUS_METER / 1000.) - kmDepth);
}

inline constexpr double km2rad(double km, double kmDepth = 0)
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

double computeAzimuth(double lat1, double lon1, double lat2, double lon2);

void computeCoordinates(double distance,
                        double azimuth,
                        double clat,
                        double clon,
                        double &lat,
                        double &lon,
                        double atKmDepth     = 0,
                        bool angularDistance = false);

double computeDistance(double lat1,
                       double lon1,
                       double lat2,
                       double lon2,
                       double *azimuth      = nullptr,
                       double *backAzimuth  = nullptr,
                       double atKmDepth     = 0,
                       bool angularDistance = false);

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

void compute5numberSummary(const std::vector<double> &values,
                           double &min,
                           double &max,
                           double &q1,
                           double &q2,
                           double &q3);

double computeMedian(const std::vector<double> &values);

double computeMedianAbsoluteDeviation(const std::vector<double> &values,
                                      const double median);

double computeMean(const std::vector<double> &values);

double computeMeanAbsoluteDeviation(const std::vector<double> &values,
                                    const double mean);

double computeCircularMean(const std::vector<double> &angles);

double interpolateSquareLagrange(double xdiff,
                                 double zdiff,
                                 double vval00,
                                 double vval01,
                                 double vval10,
                                 double vval11);

double interpolateCubeLagrange(double xdiff,
                               double ydiff,
                               double zdiff,
                               double vval000,
                               double vval001,
                               double vval010,
                               double vval011,
                               double vval100,
                               double vval101,
                               double vval110,
                               double vval111);

double interpolateSquareAngles(
    double xdiff, double zdiff, double v00, double v01, double v10, double v11);

double interpolateCubeAngles(double xdiff,
                             double ydiff,
                             double zdiff,
                             double v000,
                             double v001,
                             double v010,
                             double v011,
                             double v100,
                             double v101,
                             double v110,
                             double v111);

std::string joinPath(const std::string &path1, const std::string &path2);

bool pathExists(const std::string &path);

bool directoryEmpty(const std::string &path);

bool createDirectories(const std::string &path);

bool removePath(const std::string &path);

} // namespace HDD

#endif

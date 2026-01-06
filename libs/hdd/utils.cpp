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

#include "utils.h"
#include "csvreader.h"
#include "log.h"

#ifdef USE_BOOST_FS
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

using namespace std;
using namespace HDD::Logger;

namespace HDD {

std::string strf(const char *fmt, va_list args)
{
  constexpr size_t STACK_BUFFER_SIZE = 512;
  char stackBuffer[STACK_BUFFER_SIZE];

  va_list args_copy;
  va_copy(args_copy, args);
  int needed = vsnprintf(stackBuffer, STACK_BUFFER_SIZE, fmt, args_copy);
  va_end(args_copy);

  if (needed < 0)
  {
    return "Log formatting error";
  }

  if (static_cast<size_t>(needed) < STACK_BUFFER_SIZE)
  {
    // Fits in buffer, return directly
    return std::string(stackBuffer, needed);
  }

  // Otherwise, fallback to heap allocation with vector
  std::vector<char> heapBuffer(needed + 1);

  va_copy(args_copy, args);
  needed = vsnprintf(heapBuffer.data(), heapBuffer.size(), fmt, args_copy);
  va_end(args_copy);

  if (needed < 0)
  {
    return "Log formatting error";
  }

  return std::string(heapBuffer.data(), needed);
}

// Eventually this can be replaced with C++20 std::format
std::string strf(const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  std::string result = strf(fmt, args);
  va_end(args);
  return result;
}

std::vector<std::string> splitString(const std::string &str,
                                     const std::regex &regex)
{
  return {std::sregex_token_iterator{str.begin(), str.end(), regex, -1},
          std::sregex_token_iterator()};
}

/*
 * Computes the azimuth [radians] of the point lat2,lon2 as seen from
 * lat1,lon1
 *
 * All these formulas are for calculations on the basis of a spherical
 * earth (ignoring ellipsoidal effects)
 *
 */
double computeAzimuth(double lat1, double lon1, double lat2, double lon2)
{
  const double deltaLon = degToRad(lon2 - lon1);
  const double Fi1      = degToRad(lat1);
  const double Fi2      = degToRad(lat2);
  const double cosLat2  = cos(Fi2);

  // Handle identity cases
  if (std::abs(lat1 - lat2) < 1e-9 && std::abs(deltaLon) < 1e-9)
  {
    return 0.;
  }

  const double y = sin(deltaLon) * cosLat2;
  const double x = cos(Fi1) * sin(Fi2) - sin(Fi1) * cosLat2 * cos(deltaLon);

  const double azimuth = std::atan2(y, x);

  if (!std::isfinite(azimuth))
  {
    throw Exception("Internal logic error: computeAzimuth failed");
  }
  return azimuth;
}

/*
 * Computes the coordinates (lat, lon) of the point which is at the
 * passed azimuth [radians, clockwise from North (0 = North)] and distance
 * [radians] (isAngularDist true) or [km] (isAngularDist false and set
 * `atKmDepth` approprieatly) as seen from the point at latitude, longitude
 * [clat, clon]
 *
 * All these formulas are for calculations on the basis of a spherical
 * earth (ignoring ellipsoidal effects)
 *
 */
void computeCoordinates(double distance,
                        double azimuth,
                        double clat,
                        double clon,
                        double &lat,
                        double &lon,
                        double atKmDepth,
                        bool isAngularDist)
{
  if (distance == 0)
  {
    lat = clat;
    lon = clon;
    return;
  }

  if (!isAngularDist)
  {
    distance = km2rad(distance, atKmDepth);
  }
  clat = degToRad(clat);
  clon = degToRad(clon);

  const double cosCLat = cos(clat);
  const double sinCLat = sin(clat);
  const double cosDist = cos(distance);
  const double sinDist = sin(distance);

  lat = asin(sinCLat * cosDist + cosCLat * sinDist * cos(azimuth));
  lon = clon +
        atan2(sin(azimuth) * sinDist * cosCLat, cosDist - sinCLat * sin(lat));

  if (!std::isfinite(lat) || !std::isfinite(lon))
  {
    throw Exception("Internal logic error: computeCoordinates failed");
  }

  lat = radToDeg(lat);
  lon = normalizeLon(radToDeg(lon));
}

/*
 * Compute distance [radians] (`isAngularDist` true) or [km]
 * (`isAngularDist` false and set `atKmDepth` approprieatly) between two
 * points and optionally `azimuth` and `backazimuth` [radianss]
 *
 * All these formulas are for calculations on the basis of a spherical
 * earth (ignoring ellipsoidal effects)
 *
 */
double computeDistance(double lat1,
                       double lon1,
                       double lat2,
                       double lon2,
                       double *azimuth,
                       double *backAzimuth,
                       double atKmDepth,
                       bool isAngularDist)
{
  const double deltaLat = degToRad(lat2 - lat1);
  const double deltaLon = degToRad(lon2 - lon1);
  const double Fi1      = degToRad(lat1);
  const double Fi2      = degToRad(lat2);
  const double cosLat1  = cos(Fi1);
  const double cosLat2  = cos(Fi2);
  const double sinLat1  = sin(Fi1);
  const double sinLat2  = sin(Fi2);

  // Handle identity cases
  if (std::abs(lat1 - lat2) < 1e-9 && std::abs(deltaLon) < 1e-9)
  {
    if (azimuth) *azimuth = 0.;
    if (backAzimuth) *backAzimuth = 0.;
    return 0.;
  }

  double a = square(sin(deltaLat / 2.)) +
             cosLat1 * cosLat2 * square(sin(deltaLon / 2.));

  // Clamp 'a' to [0, 1] to avoid NaN in sqrt due to floating point precision
  a = std::max(0.0, std::min(1.0, a));

  double angularDist = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
  if (!std::isfinite(angularDist))
  {
    throw Exception("Internal logic error: computeDistance failed");
  }

  auto computeAzimuthLocal = [](double deltaLon, double cosLat1, double cosLat2,
                                double sinLat1, double sinLat2) {
    double y = sin(deltaLon) * cosLat2;
    double x = cosLat1 * sinLat2 - sinLat1 * cosLat2 * cos(deltaLon);
    return atan2(y, x);
  };

  if (azimuth)
  {
    *azimuth =
        computeAzimuthLocal(deltaLon, cosLat1, cosLat2, sinLat1, sinLat2);
    if (!std::isfinite(*azimuth))
    {
      throw Exception("Internal logic error: computeDistance failed");
    }
  }

  if (backAzimuth)
  {
    *backAzimuth = computeAzimuthLocal(degToRad(lon1 - lon2), cosLat2, cosLat1,
                                       sinLat2, sinLat1);
    if (!std::isfinite(*backAzimuth))
    {
      throw Exception("Internal logic error: computeDistance failed");
    }
  }
  return isAngularDist ? angularDist : rad2km(angularDist, atKmDepth);
}

/*
 * 3D Distance (Pythagorean approximation for depth)
 */
double computeDistance(double lat1,
                       double lon1,
                       double kmDepth1,
                       double lat2,
                       double lon2,
                       double kmDepth2,
                       double *azimuth,
                       double *backAzimuth)
{
  double avgDepth = (kmDepth1 + kmDepth2) / 2.0;
  double hDist = computeDistance(lat1, lon1, lat2, lon2, azimuth, backAzimuth,
                                 avgDepth, false);
  double vDist = abs(kmDepth1 - kmDepth2);
  if (vDist < 1e-7) return hDist;
  return std::sqrt(square(hDist) + square(vDist));
}

double computeDistance(const Catalog::Event &ev1,
                       const Catalog::Event &ev2,
                       double *azimuth,
                       double *backAzimuth)
{
  return computeDistance(ev1.latitude, ev1.longitude, ev1.depth, ev2.latitude,
                         ev2.longitude, ev2.depth, azimuth, backAzimuth);
}

double computeDistance(const Catalog::Event &event,
                       const Catalog::Station &station,
                       double *azimuth,
                       double *backAzimuth)
{
  return computeDistance(event.latitude, event.longitude, event.depth,
                         station.latitude, station.longitude,
                         -(station.elevation / 1000.), azimuth, backAzimuth);
}

std::string joinPath(const std::string &path1, const std::string &path2)
{
  return (fs::path(path1) / path2).string();
}

bool pathExists(const std::string &path)
{
  try
  {
    return fs::exists(path);
  }
  catch (exception &e)
  {
    logErrorF("%s", e.what());
    return false;
  }
}

bool directoryEmpty(const std::string &path)
{
  try
  {
    return !fs::exists(path) || (fs::is_directory(path) && fs::is_empty(path));
  }
  catch (exception &e)
  {
    logErrorF("%s", e.what());
    return false;
  }
}

bool createDirectories(const std::string &path)
{
  try
  {
    fs::create_directories(path);
    return true;
  }
  catch (exception &e)
  {
    logErrorF("%s", e.what());
    return false;
  }
}

bool removePath(const std::string &path)
{
  try
  {
    fs::remove_all(path);
    return true;
  }
  catch (exception &e)
  {
    logErrorF("%s", e.what());
    return false;
  }
}

void compute5numberSummary(const std::vector<double> &values,
                           double &min,
                           double &max,
                           double &q1,
                           double &q2,
                           double &q3)
{
  if (values.size() == 0)
  {
    min = max = q1 = q2 = q3 = 0;
    return;
  }
  vector<double> tmp(values);
  auto const Q1 = tmp.size() / 4;
  auto const Q2 = tmp.size() / 2;
  auto const Q3 = Q1 + Q2;
  std::nth_element(tmp.begin(), tmp.begin() + Q1, tmp.end());
  std::nth_element(tmp.begin() + Q1 + 1, tmp.begin() + Q2, tmp.end());
  std::nth_element(tmp.begin() + Q2 + 1, tmp.begin() + Q3, tmp.end());
  q1  = tmp[Q1];
  q2  = tmp[Q2];
  q3  = tmp[Q3];
  min = *std::min_element(tmp.begin(), tmp.begin() + Q1);
  max = *std::max_element(tmp.begin() + Q3, tmp.end());
}

double computeMedian(const std::vector<double> &values)
{
  if (values.size() == 0) return 0;

  vector<double> tmp(values);
  const auto middleItr = tmp.begin() + tmp.size() / 2;
  std::nth_element(tmp.begin(), middleItr, tmp.end());
  double median = *middleItr;
  if (tmp.size() % 2 == 0)
  {
    const auto leftMiddleItr = std::max_element(tmp.begin(), middleItr);
    median                   = (*leftMiddleItr + *middleItr) / 2;
  }
  return median;
}

double computeMedianAbsoluteDeviation(const std::vector<double> &values,
                                      const double meanOrMedian)
{
  vector<double> absoluteDeviations(values.size());
  for (unsigned i = 0; i < values.size(); i++)
  {
    absoluteDeviations[i] = std::abs(values[i] - meanOrMedian);
  }
  return computeMedian(absoluteDeviations);
}

double computeMean(const vector<double> &values)
{
  if (values.size() == 0) return 0;
  return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

double computeMeanAbsoluteDeviation(const std::vector<double> &values,
                                    const double meanOrMedian)
{
  vector<double> absoluteDeviations(values.size());
  for (unsigned i = 0; i < values.size(); i++)
  {
    absoluteDeviations[i] = std::abs(values[i] - meanOrMedian);
  }
  return computeMean(absoluteDeviations);
}

double computeCircularMean(const std::vector<double> &angles)
{
  double y = 0, x = 0;
  for (size_t i = 0; i < angles.size(); i++)
  {
    x += cos(angles[i]);
    y += sin(angles[i]);
  }
  return atan2(y / angles.size(), x / angles.size());
}

} // namespace HDD

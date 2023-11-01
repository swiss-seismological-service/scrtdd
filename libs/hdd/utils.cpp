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
#include <fstream>
#include <stdarg.h>

#ifdef USE_BOOST_FS
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

using namespace std;

namespace HDD {

// Eventually this can be replaced with C++20 std::format
// https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprint
std::string strf(const char *fmt, ...)
{
  // A static buffer that hopefully covers 99% of all use cases
  char staticBuffer[128];

  // The dynamic buffer that will be used if the static buffer is
  // not large enough
  unique_ptr<char[]> dynamicBuffer = nullptr;

  // The buffer actually written to
  char *buffer = staticBuffer;
  size_t size  = sizeof(staticBuffer);
  va_list params;
  int maxIterations = 10;

  va_start(params, fmt);
  int r = vsnprintf(buffer, size, fmt, params);
  if (r < 0)
  {
    va_end(params);
    Logger::logError("strf error");
    return std::string();
  }

  size_t requiredSize = size_t(r) + 1; // +1 for \0

  while (requiredSize > size && // create dynamic buffer with more space
         (--maxIterations >= 0))
  {
    dynamicBuffer = unique_ptr<char[]>(new char[requiredSize]);
    size          = requiredSize;
    buffer        = dynamicBuffer.get();
    *buffer       = '\0';

    va_end(params);
    va_start(params, fmt);

    r = vsnprintf(buffer, size, fmt, params);
    if (r < 0)
    {
      Logger::logError("strf error");
      break;
    }

    requiredSize = size_t(r) + 1; // +1 for \0
  }

  std::string ret(buffer);
  va_end(params);

  if (maxIterations < 0)
  {
    Logger::logError(
        "strf failed after 10 iterations: buffer still not large enough");
  }
  return ret;
}

std::vector<std::string> splitString(const std::string &str,
                                     const std::regex &regex)
{
  return {std::sregex_token_iterator{str.begin(), str.end(), regex, -1},
          std::sregex_token_iterator()};
}

/*
 * Computes the azimuth [radiant] of the point lat2,lon2 as seen from
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

  if (lat1 == lat2 && (lat1 == 90 || lat1 == -90 || deltaLon == 0))
  {
    return 0.;
  }

  const double y = sin(deltaLon) * cosLat2;
  const double x = cos(Fi1) * sin(Fi2) - sin(Fi1) * cosLat2 * cos(deltaLon);

  const double azimuth = atan2(y, x);
  if (!std::isfinite(azimuth))
  {
    throw Exception("Internal logic error: computeAzimuth failed");
  }

  return azimuth;
}

/*
 * Computes the coordinates (lat, lon) of the point which is at the
 * passed azimuth [radiant] and distance [radiant] (angularDistance true)
 * or [km] (angularDistance false and set `atKmDepth` approprieatly) as seen
 * from the point at latitude, longitude [clat, clon]
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
                        bool angularDistance)
{
  if (distance == 0)
  {
    lat = clat;
    lon = clon;
    return;
  }

  if (!angularDistance)
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
 * Compute distance [radiant] (`angularDistance` true) or [km]
 * (`angularDistance` false and set `atKmDepth` approprieatly) between two
 * points and optionally `azimuth` and `backazimuth` [radiants]
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
                       bool angularDistance)
{
  const double deltaLat = degToRad(lat2 - lat1);
  const double deltaLon = degToRad(lon2 - lon1);
  const double Fi1      = degToRad(lat1);
  const double Fi2      = degToRad(lat2);
  const double cosLat1  = cos(Fi1);
  const double cosLat2  = cos(Fi2);
  const double sinLat1  = sin(Fi1);
  const double sinLat2  = sin(Fi2);

  if (lat1 == lat2 && (lat1 == 90 || lat1 == -90 || deltaLon == 0))
  {
    if (azimuth) *azimuth = 0.;
    if (backAzimuth) *backAzimuth = 0.;
    return 0.;
  }

  double a = square(sin(deltaLat / 2.)) +
             cosLat1 * cosLat2 * square(sin(deltaLon / 2.));

  double distance = 2. * atan2(sqrt(a), sqrt(1. - a));
  if (!std::isfinite(distance))
  {
    throw Exception("Internal logic error: computeDistance failed");
  }

  auto computeAzimuth = [](double deltaLon, double cosLat1, double cosLat2,
                           double sinLat1, double sinLat2) {
    double y = sin(deltaLon) * cosLat2;
    double x = cosLat1 * sinLat2 - sinLat1 * cosLat2 * cos(deltaLon);
    return atan2(y, x);
  };

  if (azimuth)
  {
    *azimuth = computeAzimuth(deltaLon, cosLat1, cosLat2, sinLat1, sinLat2);
    if (!std::isfinite(*azimuth))
    {
      throw Exception("Internal logic error: computeDistance failed");
    }
  }

  if (backAzimuth)
  {
    *backAzimuth = computeAzimuth(degToRad(lon1 - lon2), cosLat2, cosLat1,
                                  sinLat2, sinLat1);
    if (!std::isfinite(*backAzimuth))
    {
      throw Exception("Internal logic error: computeDistance failed");
    }
  }

  return angularDistance ? distance : rad2km(distance, atKmDepth);
}

double computeDistance(double lat1,
                       double lon1,
                       double depth1,
                       double lat2,
                       double lon2,
                       double depth2,
                       double *azimuth,
                       double *backAzimuth)
{
  double atKmDepth = (depth1 + depth2) / 2.;
  double Hdist =
      computeDistance(lat1, lon1, lat2, lon2, azimuth, backAzimuth, atKmDepth);
  if (depth1 == depth2) return Hdist;
  double Vdist = abs(depth1 - depth2);
  return std::sqrt(square(Hdist) + square(Vdist));
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
    logError("%s", e.what());
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
    logError("%s", e.what());
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
    logError("%s", e.what());
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
    logError("%s", e.what());
    return false;
  }
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

void writeXCorrToFile(const XCorrCache &xcorr,
                      const Catalog &cat,
                      const std::string &file)
{
  ofstream os(file);
  os << "eventId1,eventId2,networkCode,stationCode,locationCode,component,"
        "phaseType,valid,coefficient,lag"
     << endl;

  auto callback = [&os, &cat](unsigned ev1, unsigned ev2,
                              const std::string &stationId,
                              const Catalog::Phase::Type &type,
                              const XCorrCache::Entry &e) {
    const Catalog::Station &sta = cat.getStations().at(stationId);

    os << strf("%u,%u,%s,%s,%s,%s,%c,%s,%f,%f", ev1, ev2,
               sta.networkCode.c_str(), sta.stationCode.c_str(),
               sta.locationCode.c_str(), e.component.c_str(),
               static_cast<char>(type), e.valid ? "true" : "false", e.coeff,
               e.lag)
       << endl;
  };

  xcorr.forEach(callback);
}

XCorrCache readXCorrFromFile(const Catalog &cat, const std::string &file)
{
  auto strToBool = [](const std::string &s) -> bool {
    return s == "1" || s == "true" || s == "True" || s == "TRUE";
  };
  auto strToPhaseType = [](const std::string &s) -> Catalog::Phase::Type {
    return (s == "P" || s == "p") ? Catalog::Phase::Type::P
                                  : Catalog::Phase::Type::S;
  };
  XCorrCache xcorr;
  vector<unordered_map<string, string>> phases = CSV::readWithHeader(file);
  for (const auto &row : phases)
  {
    unsigned ev1              = std::stoul(row.at("eventId1"));
    unsigned ev2              = std::stoul(row.at("eventId2"));
    string networkCode        = row.at("networkCode");
    string stationCode        = row.at("stationCode");
    string locationCode       = row.at("locationCode");
    string component          = row.at("component");
    Catalog::Phase::Type type = strToPhaseType(row.at("phaseType"));
    bool valid                = strToBool(row.at("valid"));
    double coeff              = std::stod(row.at("coefficient"));
    double lag                = std::stod(row.at("lag"));

    std::string stationId =
        cat.searchStation(networkCode, stationCode, locationCode)->second.id;

    if (!xcorr.has(ev1, ev2, stationId, type))
      xcorr.add(ev1, ev2, stationId, type, valid, coeff, lag, component);
  }
  return xcorr;
}

} // namespace HDD

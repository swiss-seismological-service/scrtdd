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

#include "utils.h"
#include "csvreader.h"
#include "geo.h"
#include "log.h"
#include <boost/filesystem.hpp>
#include <fstream>
#include <stdarg.h>

using namespace std;

namespace fs = boost::filesystem;

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
 * Computes the coordinates (lat, lon) of the point which is at the
 * passed azimuth degree and km distance as seen from the centroid
 * (clat, clon)
 */
void computeCoordinates(double distance,
                        double azimuth,
                        double clat,
                        double clon,
                        double &lat,
                        double &lon)
{
  delandaz2coord(km2deg(distance), azimuth, clat, clon, &lat, &lon);
  lon = normalizeLon(lon);
}

/*
 * Compute distance in km between two points and optionally
 * `azimuth` and `backazimuth`.
 */
double computeDistance(double lat1,
                       double lon1,
                       double depth1,
                       double lat2,
                       double lon2,
                       double depth2,
                       double *azimuth,
                       double *backAzimuth)
{
  double Hdist = computeDistance(lat1, lon1, lat2, lon2, azimuth, backAzimuth);

  if (depth1 == depth2) return Hdist;

  // Use the Euclidean distance. This approximation is sufficient when the
  // distance is small and the Earth curvature can be assumed flat.
  double Vdist = abs(depth1 - depth2);
  return std::sqrt(square(Hdist) + square(Vdist));
}

double computeDistance(double lat1,
                       double lon1,
                       double lat2,
                       double lon2,
                       double *azimuth,
                       double *backAzimuth)
{
  double dist, az, baz;
  delazi(lat1, lon1, lat2, lon2, &dist, &az, &baz);
  dist = deg2km(dist);

  if (azimuth) *azimuth = az;
  if (backAzimuth) *backAzimuth = baz;

  return dist;
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
    return !boost::filesystem::exists(path) ||
           (fs::is_directory(path) && fs::is_empty(path));
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

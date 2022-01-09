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
#include "log.h"
#include <boost/filesystem.hpp>
#include <seiscomp3/client/inventory.h>
#include <seiscomp3/math/geo.h>

using namespace std;

namespace fs = boost::filesystem;

namespace HDD {

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
  Seiscomp::Math::Geo::delandaz2coord(Seiscomp::Math::Geo::km2deg(distance),
                                      azimuth, clat, clon, &lat, &lon);
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
  Seiscomp::Math::Geo::delazi(lat1, lon1, lat2, lon2, &dist, &az, &baz);
  dist = Seiscomp::Math::Geo::deg2km(dist);

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
    return ! boost::filesystem::exists(path) ||
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
                                      const double median)
{
  vector<double> absoluteDeviations(values.size());
  for (unsigned i = 0; i < values.size(); i++)
  {
    absoluteDeviations[i] = std::abs(values[i] - median);
  }
  return computeMedian(absoluteDeviations);
}

double computeMean(const vector<double> &values)
{
  if (values.size() == 0) return 0;
  return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

double computeMeanAbsoluteDeviation(const std::vector<double> &values,
                                    const double mean)
{
  vector<double> absoluteDeviations(values.size());
  for (unsigned i = 0; i < values.size(); i++)
  {
    absoluteDeviations[i] = std::abs(values[i] - mean);
  }
  return computeMean(absoluteDeviations);
}

Seiscomp::DataModel::SensorLocation *
findSensorLocation(const std::string &networkCode,
                   const std::string &stationCode,
                   const std::string &locationCode,
                   const Time &atTime)
{
  Seiscomp::DataModel::Inventory *inv =
      Seiscomp::Client::Inventory::Instance()->inventory();
  if (!inv)
  {
    logDebug("Inventory not available");
    return nullptr;
  }

  Seiscomp::DataModel::InventoryError error;
  Seiscomp::DataModel::SensorLocation *loc =
      Seiscomp::DataModel::getSensorLocation(inv, networkCode, stationCode,
                                             locationCode, atTime, &error);

  if (!loc)
  {
    logDebug("Unable to fetch SensorLocation information (%s.%s.%s at %s): %s",
             networkCode.c_str(), stationCode.c_str(), locationCode.c_str(),
             atTime.iso().c_str(), error.toString());
  }
  return loc;
}

} // namespace HDD

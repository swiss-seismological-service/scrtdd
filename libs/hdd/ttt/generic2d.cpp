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

#include "generic2d.h"
#include "../csvreader.h"
#include "../log.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

#ifdef USE_BOOST_FS
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

namespace HDD {
namespace TTT {

using namespace std;
using namespace HDD::Logger;

const string Generic2D::COL_ELEVATION = "StationElevation[m]";
const string Generic2D::COL_DISTDEG   = "Distance[deg]";
const string Generic2D::COL_DISTKM    = "Distance[km]";
const string Generic2D::COL_DEPTH     = "Depth[km]";
const string Generic2D::COL_TTIME     = "TravelTime[sec]";
const string Generic2D::COL_TAKEOFF   = "Takeoff[deg]"; // 0(down):180(up)
const string Generic2D::COL_DTDD      = "dtdd[sec/deg]";
const string Generic2D::COL_DTDH      = "dtdh[sec/km]";
const string Generic2D::COL_VELOCITY  = "SourceVelocity[km/sec]";

Generic2D::Generic2D(const std::string &gridPath, const std::string &gridModel)
    : _gridPath(gridPath), _gridModel(gridModel)
{
  logDebugF("Searching for travel times in %s model name %s", _gridPath.c_str(),
            _gridModel.c_str());

  const double nan = std::numeric_limits<double>::quiet_NaN();

  struct RowEntry
  {
    double dist;
    double depth;
    double ttime;
    double takeoff;
    double dtdd;
    double dtdh;
  };
  // keys phase, elevation
  unordered_map<string, map<double, vector<RowEntry>>> elevationMap;

  //
  // Read all model files and store the data in elevationMap
  //
  for (const auto &entry : fs::directory_iterator(_gridPath))
  {
    const string filename = entry.path().filename().string();
    if (filename.rfind(_gridModel + ".", 0) == 0)
    {
      std::size_t pos = filename.find_last_of('.');
      if (pos == std::string::npos || pos == 0)
      {
        // No extension or hidden file like ".filename"
        logWarning("No phase in file name: " + filename);
        continue;
      }
      const string phase = filename.substr(pos + 1);

      vector<unordered_map<string, string>> csvData =
          CSV::readWithHeader(entry.path().string());
      if (csvData.empty())
      {
        logWarning("CSV file is empty: " + filename);
        continue;
      }

      //
      // Group rows by station elevation, and verify mandatory fields
      //
      unsigned row_count = 0;
      for (const auto &row : csvData)
      {
        row_count++;
        try
        {
          RowEntry re;
          auto it = row.find(COL_DISTDEG);
          if (it != row.end()) // Check if "Distance[deg]" exists
          {
            re.dist = std::stod(it->second);
          }
          else
          {
            re.dist = km2deg(std::stod(row.at(COL_DISTKM)));
          }
          re.depth = std::stod(row.at(COL_DEPTH));
          re.ttime = std::stod(row.at(COL_TTIME));
          if (re.ttime < 0)
          {
            // Special case: negative ttime means no data for this dist/depth
            re.ttime   = nan;
            re.takeoff = nan;
            re.dtdd    = nan;
            re.dtdh    = nan;
          }
          else
          {
            // Read other columns and fill in the data
            re.takeoff = std::stod(row.at(COL_TAKEOFF));
            it         = row.find(COL_VELOCITY);
            if (it == row.end())
            {
              re.dtdd = std::stod(row.at(COL_DTDD));
              re.dtdh = std::stod(row.at(COL_DTDH));
            }
            else
            {
              // if velocity is provided, we convert it to partial derivatives
              const double velocity =
                  std::stod(row.at(COL_VELOCITY)); // [km/sec]
              const double dip = degToRad(re.takeoff - 90);
              re.dtdd          = std::cos(dip) / km2deg(velocity); // [sec/deg]
              re.dtdh          = std::sin(dip) / velocity;         // [sec/km]
            }
          }
          double elevation = std::stod(row.at(COL_ELEVATION));
          elevationMap[phase][elevation].push_back(std::move(re));
        }
        catch (const std::exception &e)
        {
          logWarningF("Error while parsing file '%s' at row %d: %s",
                      filename.c_str(), row_count, e.what());
          continue;
        }
      }
    }
  }

  for (auto const &[phase, emap] : elevationMap)
  {
    for (auto const &[elevation, rows] : emap)
    {
      GridData gridData;
      set<double> distancesSet, depthsSet;

      for (const auto &re : rows)
      {
        distancesSet.insert(re.dist);
        depthsSet.insert(re.depth);
      }

      gridData.distances.assign(distancesSet.begin(), distancesSet.end());
      gridData.depths.assign(depthsSet.begin(), depthsSet.end());

      size_t nDist  = gridData.distances.size();
      size_t nDepth = gridData.depths.size();

      if (nDist == 0 || nDepth == 0) continue;

      gridData.travelTimes.resize(nDist, vector<double>(nDepth, nan));
      gridData.takeoffDips.resize(nDist, vector<double>(nDepth, nan));
      gridData.dtdds.resize(nDist, vector<double>(nDepth, nan));
      gridData.dtdhs.resize(nDist, vector<double>(nDepth, nan));

      map<double, int> distMap, depthMap;
      for (size_t i = 0; i < nDist; ++i) distMap[gridData.distances[i]] = i;
      for (size_t i = 0; i < nDepth; ++i) depthMap[gridData.depths[i]] = i;

      for (const auto &re : rows)
      {
        int distIdx  = distMap.at(re.dist);
        int depthIdx = depthMap.at(re.depth);

        gridData.travelTimes[distIdx][depthIdx] = re.ttime;
        gridData.takeoffDips[distIdx][depthIdx] = re.takeoff;
        gridData.dtdds[distIdx][depthIdx]       = re.dtdd;
        gridData.dtdhs[distIdx][depthIdx]       = re.dtdh;
      }
      _phaseData[phase][elevation] = std::move(gridData);
    }
  }
  logDebugF("Loaded %zu phase files", _phaseData.size());
}

std::unordered_set<std::string> Generic2D::availablePhases() const
{
  std::unordered_set<string> phases;
  phases.reserve(_phaseData.size());
  for (const auto &kv : _phaseData)
  {
    phases.insert(kv.first);
  }
  return phases;
}

void Generic2D::findIndexes(const std::string &phaseType,
                            double stationElevation,
                            double distance_deg,
                            double eventDepth,
                            const GridData *&grid,
                            size_t &i_dist1,
                            double &dist_frac,
                            size_t &i_depth1,
                            double &depth_frac) const
{
  auto it_phase = _phaseData.find(phaseType);
  if (it_phase == _phaseData.end())
  {
    throw Exception("No data for phase: " + phaseType);
  }

  const auto &elevation_map = it_phase->second;
  if (elevation_map.empty())
  {
    throw Exception("No data for phase: " + phaseType);
  }

  // Find closest elevation
  auto it_grid = elevation_map.lower_bound(stationElevation);
  if (it_grid == elevation_map.end())
  {
    // stationElevation is greater than all keys
    it_grid--;
  }
  else if (it_grid != elevation_map.begin())
  {
    // find the closest elevation between the greater and smaller neighbours
    auto prev = it_grid;
    --prev;
    if ((it_grid->first - stationElevation) > (stationElevation - prev->first))
    {
      it_grid = prev;
    }
  }
  grid = &it_grid->second;

  // Find indices and diff for distance
  auto it_dist = std::lower_bound(grid->distances.begin(),
                                  grid->distances.end(), distance_deg);
  if (it_dist == grid->distances.end() ||
      (it_dist == grid->distances.begin() && distance_deg < *it_dist))
  {
    throw OutOfRange("Distance out of bounds.");
  }

  if (it_dist == grid->distances.begin())
  {
    // We are at the first grid point.
    // The above check ensures distance_deg == *it_dist
    i_dist1 = 0;
  }
  else
  {
    i_dist1 = std::distance(grid->distances.begin(), it_dist) - 1;
  }
  if (i_dist1 + 1 >= grid->distances.size())
  {
    throw OutOfRange("Distance out of bounds.");
  }
  size_t i_dist2   = i_dist1 + 1;
  double dist_diff = grid->distances[i_dist2] - grid->distances[i_dist1];
  if (dist_diff <= 0)
  {
    dist_frac = 0;
  }
  else
  {
    dist_frac = (distance_deg - grid->distances[i_dist1]) / dist_diff;
  }

  // Find indices and diff for depth
  auto it_depth =
      std::lower_bound(grid->depths.begin(), grid->depths.end(), eventDepth);
  if (it_depth == grid->depths.end() ||
      (it_depth == grid->depths.begin() && eventDepth < *it_depth))
  {
    throw OutOfRange("Depth out of bounds.");
  }

  if (it_depth == grid->depths.begin())
  {
    i_depth1 = 0;
  }
  else
  {
    i_depth1 = std::distance(grid->depths.begin(), it_depth) - 1;
  }
  if (i_depth1 + 1 >= grid->depths.size())
  {
    throw OutOfRange("Depth out of bounds.");
  }
  size_t i_depth2   = i_depth1 + 1;
  double depth_diff = grid->depths[i_depth2] - grid->depths[i_depth1];
  if (depth_diff <= 0)
  {
    depth_frac = 0;
  }
  else
  {
    depth_frac = (eventDepth - grid->depths[i_depth1]) / depth_diff;
  }
}

double Generic2D::compute(double eventLat,
                          double eventLon,
                          double eventDepth,
                          double stationLat,
                          double stationLon,
                          double stationElevation,
                          const std::string &phaseType)
{
  const GridData *grid;
  size_t i_dist1, i_depth1;
  double dist_frac, depth_frac;
  double distance_deg = radToDeg(computeDistance(
      eventLat, eventLon, stationLat, stationLon, nullptr, nullptr, 0, true));

  findIndexes(phaseType, stationElevation, distance_deg, eventDepth, grid,
              i_dist1, dist_frac, i_depth1, depth_frac);

  const size_t i_dist2  = i_dist1 + 1;
  const size_t i_depth2 = i_depth1 + 1;

  const double tt11 = grid->travelTimes[i_dist1][i_depth1];
  const double tt12 = grid->travelTimes[i_dist1][i_depth2];
  const double tt21 = grid->travelTimes[i_dist2][i_depth1];
  const double tt22 = grid->travelTimes[i_dist2][i_depth2];

  if (std::isnan(tt11) || std::isnan(tt12) || std::isnan(tt21) ||
      std::isnan(tt22))
  {
    throw OutOfRange("Grid data is sparse. Cannot interpolate.");
  }

  return interpolateBilinear(dist_frac, depth_frac, tt11, tt12, tt21, tt22);
}

void Generic2D::compute(double eventLat,
                        double eventLon,
                        double eventDepth,
                        double stationLat,
                        double stationLon,
                        double stationElevation,
                        const std::string &phaseType,
                        double &travelTime,
                        double &takeOffAzi,
                        double &takeOffDip,
                        double &dtdd,
                        double &dtdh)
{
  const GridData *grid;
  size_t i_dist1, i_depth1;
  double dist_frac, depth_frac;
  double distance_deg = radToDeg(computeDistance(
      eventLat, eventLon, stationLat, stationLon, nullptr, nullptr, 0, true));

  findIndexes(phaseType, stationElevation, distance_deg, eventDepth, grid,
              i_dist1, dist_frac, i_depth1, depth_frac);

  const size_t i_dist2  = i_dist1 + 1;
  const size_t i_depth2 = i_depth1 + 1;

  const double tt11 = grid->travelTimes[i_dist1][i_depth1];
  const double tt12 = grid->travelTimes[i_dist1][i_depth2];
  const double tt21 = grid->travelTimes[i_dist2][i_depth1];
  const double tt22 = grid->travelTimes[i_dist2][i_depth2];

  const double td11 = grid->takeoffDips[i_dist1][i_depth1];
  const double td12 = grid->takeoffDips[i_dist1][i_depth2];
  const double td21 = grid->takeoffDips[i_dist2][i_depth1];
  const double td22 = grid->takeoffDips[i_dist2][i_depth2];

  const double dd11 = grid->dtdds[i_dist1][i_depth1];
  const double dd12 = grid->dtdds[i_dist1][i_depth2];
  const double dd21 = grid->dtdds[i_dist2][i_depth1];
  const double dd22 = grid->dtdds[i_dist2][i_depth2];

  const double dh11 = grid->dtdhs[i_dist1][i_depth1];
  const double dh12 = grid->dtdhs[i_dist1][i_depth2];
  const double dh21 = grid->dtdhs[i_dist2][i_depth1];
  const double dh22 = grid->dtdhs[i_dist2][i_depth2];

  if (std::isnan(tt11) || std::isnan(tt12) || std::isnan(tt21) ||
      std::isnan(tt22) || std::isnan(td11) || std::isnan(td12) ||
      std::isnan(td21) || std::isnan(td22) || std::isnan(dd11) ||
      std::isnan(dd12) || std::isnan(dd21) || std::isnan(dd22) ||
      std::isnan(dh11) || std::isnan(dh12) || std::isnan(dh21) ||
      std::isnan(dh22))
  {
    throw OutOfRange("Grid data is sparse. Cannot interpolate.");
  }

  travelTime =
      interpolateBilinear(dist_frac, depth_frac, tt11, tt12, tt21, tt22);

  takeOffDip =
      interpolateBilinear(dist_frac, depth_frac, td11, td12, td21, td22);

  dtdd = interpolateBilinear(dist_frac, depth_frac, dd11, dd12, dd21, dd22);

  dtdh = interpolateBilinear(dist_frac, depth_frac, dh11, dh12, dh21, dh22);

  takeOffAzi =
      radToDeg(computeAzimuth(eventLat, eventLon, stationLat, stationLon));
}

} // namespace TTT
} // namespace HDD

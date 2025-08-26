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

#include "nllgrid.h"
#include "../log.h"
#include "../utils.h"

#include <regex>

#ifdef USE_BOOST_FS
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

using namespace std;
using namespace HDD::NLL;

namespace HDD {

namespace TTT {

NLLGrid::GridCloser::GridCloser(const shared_ptr<VelGrid> &grid)
{
  velGrid = grid;
}

NLLGrid::GridCloser::GridCloser(const shared_ptr<TimeGrid> &grid)
{
  timeGrid = grid;
}

NLLGrid::GridCloser::GridCloser(const shared_ptr<AngleGrid> &grid)
{
  angleGrid = grid;
}

NLLGrid::GridCloser::~GridCloser()
{
  if (velGrid) velGrid->close();
  if (timeGrid) timeGrid->close();
  if (angleGrid) angleGrid->close();
}

NLLGrid::NLLGrid(const std::string &griddPath,
                 const std::string &gridModel,
                 bool swapBytes,
                 unsigned maxOpenFiles)
    : _gridPath(griddPath), _gridModel(gridModel), _swapBytes(swapBytes),
      _openGrids(maxOpenFiles)
{

  static const std::regex rePModel(_gridModel + R"(\.P\.mod\.hdr)",
                                   std::regex::optimize);
  static const std::regex reSModel(_gridModel + R"(\.S\.mod\.hdr)",
                                   std::regex::optimize);
  static const std::regex rePTime(_gridModel + R"(\.P\..+\.time\.hdr)",
                                  std::regex::optimize);
  static const std::regex reSTime(_gridModel + R"(\.S\..+\.time\.hdr)",
                                  std::regex::optimize);
  static const std::regex rePAngle(_gridModel + R"(\.P\..+\.angle\.hdr)",
                                   std::regex::optimize);
  static const std::regex reSAngle(_gridModel + R"(\.S\..+\.angle\.hdr)",
                                   std::regex::optimize);

  logDebugF("Searching grids in %s model name %s", _gridPath.c_str(),
            _gridModel.c_str());

  vector<TimeKDTree::Point> PTimeGrids;
  vector<TimeKDTree::Point> STimeGrids;
  vector<AngleKDTree::Point> PAngleGrids;
  vector<AngleKDTree::Point> SAngleGrids;

  for (const auto &entry : fs::directory_iterator(_gridPath))
  {

    const string filename = entry.path().filename().string();
    const string baseFilePath =
        fs::path(entry.path()).replace_extension().string();
    try
    {
      if (std::regex_match(filename, rePTime))
      {
        TimeKDTree::Point point;
        auto grid  = make_shared<TimeGrid>(baseFilePath, _swapBytes);
        point.data = grid;
        const Grid::Info &ginfo = grid->getInfo();
        ginfo.transform->toLatLon(ginfo.dx, ginfo.dy, point.latitude,
                                  point.longitude);
        point.elevation = -ginfo.dz;
        PTimeGrids.emplace_back(point);
      }
      else if (std::regex_match(filename, reSTime))
      {
        TimeKDTree::Point point;
        auto grid  = make_shared<TimeGrid>(baseFilePath, _swapBytes);
        point.data = grid;
        const Grid::Info &ginfo = grid->getInfo();
        ginfo.transform->toLatLon(ginfo.dx, ginfo.dy, point.latitude,
                                  point.longitude);
        point.elevation = -ginfo.dz;
        STimeGrids.emplace_back(point);
      }
      else if (std::regex_match(filename, rePAngle))
      {
        AngleKDTree::Point point;
        auto grid  = make_shared<AngleGrid>(baseFilePath, _swapBytes);
        point.data = grid;
        const Grid::Info &ginfo = grid->getInfo();
        ginfo.transform->toLatLon(ginfo.dx, ginfo.dy, point.latitude,
                                  point.longitude);
        point.elevation = -ginfo.dz;
        PAngleGrids.emplace_back(point);
      }
      else if (std::regex_match(filename, reSAngle))
      {
        AngleKDTree::Point point;
        auto grid  = make_shared<AngleGrid>(baseFilePath, _swapBytes);
        point.data = grid;
        const Grid::Info &ginfo = grid->getInfo();
        ginfo.transform->toLatLon(ginfo.dx, ginfo.dy, point.latitude,
                                  point.longitude);
        point.elevation = -ginfo.dz;
        SAngleGrids.emplace_back(point);
      }
      else if (std::regex_match(filename, rePModel))
      {
        _PVelGrid = make_shared<VelGrid>(baseFilePath, _swapBytes);
      }
      else if (std::regex_match(filename, reSModel))
      {
        _SVelGrid = make_shared<VelGrid>(baseFilePath, _swapBytes);
      }
    }
    catch (exception &e)
    {
      logWarningF("Cannot load grid file: %s", e.what());
    }
  }

  logDebugF("Read %d/%d models %zu/%zu time %zu/%zu angle P/S grid files",
            _PVelGrid ? 1 : 0, _SVelGrid ? 1 : 0, PTimeGrids.size(),
            STimeGrids.size(), PAngleGrids.size(), SAngleGrids.size());

  _PTimeKDTree  = KDTree<shared_ptr<TimeGrid>>(PTimeGrids);
  _STimeKDTree  = KDTree<shared_ptr<TimeGrid>>(STimeGrids);
  _PAngleKDTree = KDTree<shared_ptr<AngleGrid>>(PAngleGrids);
  _SAngleKDTree = KDTree<shared_ptr<AngleGrid>>(SAngleGrids);
}

void NLLGrid::closeOpenFiles() { _openGrids.clear(); }

double NLLGrid::compute(double eventLat,
                        double eventLon,
                        double eventDepth,
                        double stationLat,
                        double stationLon,
                        double stationElevation,
                        const std::string &phaseType)
{
  string key = "time" + phaseType + to_string(stationLat) + ":" +
               to_string(stationLon) + ":" + to_string(stationElevation);

  shared_ptr<NLL::TimeGrid> timeGrid;
  try
  {
    // try the cached grid first
    timeGrid = _openGrids.get(key).timeGrid;
  }
  catch (std::range_error &e)
  {
    try
    {
      // Find the correct grid
      if (phaseType == "P")
      {
        timeGrid =
            _PTimeKDTree.search(stationLat, stationLon, stationElevation).data;
      }
      else if (phaseType == "S")
      {
        timeGrid =
            _STimeKDTree.search(stationLat, stationLon, stationElevation).data;
      }
    }
    catch (std::range_error &e)
    {
      string msg =
          strf("Cannot find a suitable %s time grid (station lat %g lon "
               "%g elevation %g)",
               phaseType.c_str(), stationLat, stationLon, stationElevation);
      throw Exception(msg);
    }

    // cache the grid
    _openGrids.put(key, GridCloser(timeGrid));
  }
  return timeGrid->getTime(eventLat, eventLon, eventDepth);
}

void NLLGrid::compute(double eventLat,
                      double eventLon,
                      double eventDepth,
                      double stationLat,
                      double stationLon,
                      double stationElevation,
                      const std::string &phaseType,
                      double &travelTime,
                      double &azimuth,
                      double &takeOffAngle,
                      double &velocityAtSrc)
{
  //
  // get travelTime
  //
  travelTime = compute(eventLat, eventLon, eventDepth, stationLat, stationLon,
                       stationElevation, phaseType);

  //
  // Get velocityAtSrc
  //
  string key = "velocity" + phaseType;

  shared_ptr<NLL::VelGrid> velGrid;
  try
  {
    // try the cached grid first (mostly to refresh the lru cache)
    velGrid = _openGrids.get(key).velGrid;
  }
  catch (std::range_error &e)
  {
    // Fetch the correct grid
    if (phaseType == "P")
    {
      velGrid = _PVelGrid;
    }
    else if (phaseType == "S")
    {
      velGrid = _SVelGrid;
    }

    if (!velGrid)
    {
      string msg =
          strf("Cannot find a suitable %s velocity grid (station lat %g "
               "lon %g elevation %g)",
               phaseType.c_str(), stationLat, stationLon, stationElevation);
      throw Exception(msg);
    }

    // cache the grid (to keep track of open files)
    _openGrids.put(key, GridCloser(velGrid));
  }

  // get the value from the grid now that it is loaded
  velocityAtSrc = velGrid->getVel(eventLat, eventLon, eventDepth);

  //
  // Get takeOffAngles
  //
  key = "angle" + phaseType + to_string(stationLat) + ":" +
        to_string(stationLon) + ":" + to_string(stationElevation);

  shared_ptr<NLL::AngleGrid> angleGrid;
  try
  {
    // try the cached grid first
    angleGrid = _openGrids.get(key).angleGrid;
  }
  catch (std::range_error &e)
  {
    try
    {
      // Find the correct grid
      if (phaseType == "P")
      {
        angleGrid =
            _PAngleKDTree.search(stationLat, stationLon, stationElevation).data;
      }
      else if (phaseType == "S")
      {
        angleGrid =
            _SAngleKDTree.search(stationLat, stationLon, stationElevation).data;
      }
    }
    catch (std::range_error &e)
    {
      string msg =
          strf("Cannot find a suitable %s angle grid (station lat %g "
               "lon %g elevation %g)",
               phaseType.c_str(), stationLat, stationLon, stationElevation);
      throw Exception(msg);
    }

    // cache the grid
    _openGrids.put(key, GridCloser(angleGrid));
  }

  angleGrid->getAngles(eventLat, eventLon, eventDepth, azimuth, takeOffAngle);

  if (!std::isfinite(azimuth)) // if not 3D model compute azimuth
  {
    azimuth = computeAzimuth(eventLat, eventLon, stationLat, stationLon);
  }
}

} // namespace TTT

} // namespace HDD

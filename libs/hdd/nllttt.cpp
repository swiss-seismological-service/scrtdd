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

#include "nllttt.h"
#include "log.h"
#include "utils.h"

using namespace std;

namespace HDD {
namespace NLL {

TravelTimeTable::TravelTimeTable(const std::string &velGridPath,
                                 const std::string &timeGridPath,
                                 const std::string &angleGridPath,
                                 bool swapBytes,
                                 unsigned cacheSize)
    : _velGridPath(velGridPath), _timeGridPath(timeGridPath),
      _angleGridPath(angleGridPath), _swapBytes(swapBytes),
      _velGrids(cacheSize), _timeGrids(cacheSize), _angleGrids(cacheSize)
{}

void TravelTimeTable::freeResources()
{
  _velGrids.clear();
  _timeGrids.clear();
  _angleGrids.clear();
}

std::string TravelTimeTable::filePath(const std::string &basePath,
                                      const Catalog::Station &station,
                                      const std::string &phaseType)
{
  static const std::regex reNet("@NETWORK@", std::regex::optimize);
  static const std::regex reSta("@STATION@", std::regex::optimize);
  static const std::regex reLoc("@LOCATION@", std::regex::optimize);
  static const std::regex rePha("@PHASE@", std::regex::optimize);

  string out = std::regex_replace(basePath, reNet, station.networkCode);
  out        = std::regex_replace(out, reSta, station.stationCode);
  out        = std::regex_replace(out, reLoc, station.locationCode);
  return std::regex_replace(out, rePha, phaseType);
}

double TravelTimeTable::compute(double eventLat,
                                double eventLon,
                                double eventDepth,
                                const Catalog::Station &station,
                                const std::string &phaseType)
{
  string timeGridFile = filePath(_timeGridPath, station, phaseType);
  string timeGId      = "timeGrid:" + timeGridFile;
  double travelTime;
  try
  {
    travelTime =
        _timeGrids.get(timeGId)->getTime(eventLat, eventLon, eventDepth);
  }
  catch (std::range_error &e)
  {
    // Check if we have already excluded the grid because we couldn't load it
    if (_unloadableGrids.find(timeGId) != _unloadableGrids.end())
    {
      string msg = strf("Time grid (%s) not avaliable", timeGId.c_str());
      throw Exception(msg.c_str());
    }

    // Load the grid
    try
    {
      _timeGrids.put(timeGId, shared_ptr<TimeGrid>(
                                  new TimeGrid(timeGridFile, _swapBytes)));
    }
    catch (exception &e)
    {
      _unloadableGrids.insert(timeGId);
      throw Exception(e.what());
    }

    travelTime =
        _timeGrids.get(timeGId)->getTime(eventLat, eventLon, eventDepth);
  }
  return travelTime;
}

void TravelTimeTable::compute(double eventLat,
                              double eventLon,
                              double eventDepth,
                              const Catalog::Station &station,
                              const std::string &phaseType,
                              double &travelTime,
                              double &azimuth,
                              double &takeOffAngle,
                              double &velocityAtSrc)
{
  // get travelTime
  travelTime = compute(eventLat, eventLon, eventDepth, station, phaseType);

  // Get velocityAtSrc
  string velGridFile = filePath(_velGridPath, station, phaseType);
  string velGId      = "velGrid:" + velGridFile;
  try
  {
    velocityAtSrc =
        _velGrids.get(velGId)->getVel(eventLat, eventLon, eventDepth);
  }
  catch (std::range_error &e)
  {
    // Check if we have already excluded the grid because we couldn't load it
    if (_unloadableGrids.find(velGId) != _unloadableGrids.end())
    {
      string msg = strf("Vel grid (%s) not avaliable", velGId.c_str());
      throw Exception(msg.c_str());
    }

    // Load the grid
    try
    {
      _velGrids.put(velGId,
                    shared_ptr<VelGrid>(new VelGrid(velGridFile, _swapBytes)));
    }
    catch (exception &e)
    {
      _unloadableGrids.insert(velGId);
      throw Exception(e.what());
    }

    // Again, get the value from the grid now that it is loaded
    velocityAtSrc =
        _velGrids.get(velGId)->getVel(eventLat, eventLon, eventDepth);
  }

  // Get takeOffAngles
  azimuth      = std::numeric_limits<double>::quiet_NaN();
  takeOffAngle = std::numeric_limits<double>::quiet_NaN();

  string angleGridFile = filePath(_angleGridPath, station, phaseType);
  string angleGId      = "angleGrid:" + angleGridFile;

  try
  {
    _angleGrids.get(angleGId)->getAngles(eventLat, eventLon, eventDepth,
                                         azimuth, takeOffAngle);
  }
  catch (std::range_error &e)
  {
    // Check if we have already excluded the grid because we couldn't load it
    if (_unloadableGrids.find(angleGId) != _unloadableGrids.end())
    {
      string msg = strf("Time grid (%s) not avaliable", angleGId.c_str());
      throw Exception(msg.c_str());
    }

    // Load the grid
    try
    {
      _angleGrids.put(angleGId, shared_ptr<AngleGrid>(
                                    new AngleGrid(angleGridFile, _swapBytes)));
    }
    catch (exception &e)
    {
      _unloadableGrids.insert(angleGId);
      throw Exception(e.what());
    }
    _angleGrids.get(angleGId)->getAngles(eventLat, eventLon, eventDepth,
                                         azimuth, takeOffAngle);
  }

  // approximate angles if not already provided by the grid
  computeApproximatedTakeOffAngles(
      eventLat, eventLon, eventDepth, station, phaseType,
      std::isfinite(azimuth) ? nullptr : &azimuth,
      std::isfinite(takeOffAngle) ? nullptr : &takeOffAngle);
}

} // namespace NLL
} // namespace HDD

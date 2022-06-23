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

void TravelTimeTable::compute(double eventLat,
                              double eventLon,
                              double eventDepth,
                              const Catalog::Station &station,
                              const std::string &phaseType,
                              double &travelTime)
{
  string timeGId =
      "timeGrid:" + Grid::filePath(_timeGridPath, station, phaseType);
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
      _timeGrids.put(timeGId,
                     shared_ptr<TimeGrid>(new TimeGrid(_timeGridPath, station,
                                                       phaseType, _swapBytes)));
    }
    catch (exception &e)
    {
      _unloadableGrids.insert(timeGId);
      throw Exception(e.what());
    }

    travelTime =
        _timeGrids.get(timeGId)->getTime(eventLat, eventLon, eventDepth);
  }
}

void TravelTimeTable::compute(double eventLat,
                              double eventLon,
                              double eventDepth,
                              const Catalog::Station &station,
                              const std::string &phaseType,
                              double &travelTime,
                              double &takeOffAngleAzim,
                              double &takeOffAngleDip,
                              double &velocityAtSrc)
{
  // get travelTime
  compute(eventLat, eventLon, eventDepth, station, phaseType, travelTime);

  // Get velocityAtSrc
  string velGId = "velGrid:" + Grid::filePath(_velGridPath, station, phaseType);
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
      _velGrids.put(velGId, shared_ptr<VelGrid>(new VelGrid(
                                _velGridPath, station, phaseType, _swapBytes)));
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
  takeOffAngleAzim = std::numeric_limits<double>::quiet_NaN();
  takeOffAngleDip  = std::numeric_limits<double>::quiet_NaN();

  string angleGId =
      "angleGrid:" + Grid::filePath(_angleGridPath, station, phaseType);

  try
  {
    _angleGrids.get(angleGId)->getAngles(eventLat, eventLon, eventDepth,
                                         takeOffAngleAzim, takeOffAngleDip);
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
      _angleGrids.put(angleGId,
                      shared_ptr<AngleGrid>(new AngleGrid(
                          _angleGridPath, station, phaseType, _swapBytes)));
    }
    catch (exception &e)
    {
      _unloadableGrids.insert(angleGId);
      throw Exception(e.what());
    }
    _angleGrids.get(angleGId)->getAngles(eventLat, eventLon, eventDepth,
                                         takeOffAngleAzim, takeOffAngleDip);
  }

  // approximate angles if not already provided by the grid
  computeApproximatedTakeOfAngles(
      eventLat, eventLon, eventDepth, station, phaseType,
      std::isfinite(takeOffAngleAzim) ? nullptr : &takeOffAngleAzim,
      std::isfinite(takeOffAngleDip) ? nullptr : &takeOffAngleDip);
}

} // namespace NLL
} // namespace HDD

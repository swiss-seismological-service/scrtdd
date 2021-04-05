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
#include "utils.h"

#include <array>
#include <cstring>
#include <seiscomp3/core/strings.h>
#include <seiscomp3/math/math.h>
#include <seiscomp3/utils/files.h>
#include <stdexcept>

#define SEISCOMP_COMPONENT HDD
#include <seiscomp3/logging/log.h>

using namespace std;
using Seiscomp::Core::stringify;
using TakeOffAngles = Seiscomp::HDD::NLL::AngleGrid::TakeOffAngles;

namespace {

/*
 * function to find value inside a square using Lagrange interpolation
 * 0.0 <= xdiff/zdiff <= 1.0
 * returns interp value at(xdiff, zdiff)
 */
double interpolateSquareLagrange(double xdiff,
                                 double zdiff,
                                 double vval00,
                                 double vval01,
                                 double vval10,
                                 double vval11)
{
  return vval00 * (1.0 - xdiff) * (1.0 - zdiff) +
         vval01 * (1.0 - xdiff) * zdiff + vval10 * xdiff * (1.0 - zdiff) +
         vval11 * xdiff * zdiff;
}

/*
 * function to find value inside a cube using Lagrange interpolation
 * 0.0 <= xdiff/ydiff/zdiff <= 1.0
 * returns interp value at(xdiff, ydiff, zdiff)
 */
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
                               double vval111)
{
  double oneMinusXdiff = 1.0 - xdiff;
  double oneMinusYdiff = 1.0 - ydiff;
  double oneMinusZdiff = 1.0 - zdiff;

  return vval000 * (oneMinusXdiff) * (oneMinusYdiff) * (oneMinusZdiff) +
         vval001 * (oneMinusXdiff) * (oneMinusYdiff)*zdiff +
         vval010 * (oneMinusXdiff)*ydiff * (oneMinusZdiff) +
         vval011 * (oneMinusXdiff)*ydiff * zdiff +
         vval100 * xdiff * (oneMinusYdiff) * (oneMinusZdiff) +
         vval101 * xdiff * (oneMinusYdiff)*zdiff +
         vval110 * xdiff * ydiff * (oneMinusZdiff) +
         vval111 * xdiff * ydiff * zdiff;
}

/*
 * function to find angles inside a cube
 * angle_quality_cutoff value to determine "bad" angles
 * 0.0 <= xdiff/ydiff/zdiff <= 1.0
 */
TakeOffAngles interpolateCubeAngles(double xdiff,
                                    double ydiff,
                                    double zdiff,
                                    TakeOffAngles vval000,
                                    TakeOffAngles vval001,
                                    TakeOffAngles vval010,
                                    TakeOffAngles vval011,
                                    TakeOffAngles vval100,
                                    TakeOffAngles vval101,
                                    TakeOffAngles vval110,
                                    TakeOffAngles vval111,
                                    unsigned angle_quality_cutoff)
{
  // check for lowest quality angles
  unsigned short lowest_qual = std::min(std::initializer_list<unsigned short>(
      {vval000.quality, vval001.quality, vval010.quality, vval011.quality,
       vval100.quality, vval101.quality, vval110.quality, vval111.quality}));

  // if lowest quality is too low, use nearest node */
  if (lowest_qual < angle_quality_cutoff)
  {
    TakeOffAngles nearest;
    if (xdiff < 0.5)
      nearest = (ydiff < 0.5) ? ((zdiff < 0.5) ? vval000 : vval001)
                              : ((zdiff < 0.5) ? vval010 : vval011);
    else
      nearest = (ydiff < 0.5) ? ((zdiff < 0.5) ? vval100 : vval101)
                              : ((zdiff < 0.5) ? vval110 : vval111);

    if (nearest.quality > lowest_qual)
    {
      return nearest;
    }
  }

  /* otherwise interpolate */
  unsigned short azim_interp = interpolateCubeLagrange(
      xdiff, ydiff, zdiff, vval000.azimuth, vval001.azimuth, vval010.azimuth,
      vval011.azimuth, vval100.azimuth, vval101.azimuth, vval110.azimuth,
      vval111.azimuth);
  unsigned short dip_interp = interpolateCubeLagrange(
      xdiff, ydiff, zdiff, vval000.dip, vval001.dip, vval010.dip, vval011.dip,
      vval100.dip, vval101.dip, vval110.dip, vval111.dip);

  return {lowest_qual, dip_interp, azim_interp};
}

} // namespace

namespace Seiscomp {
namespace HDD {
namespace NLL {

NllTravelTimeTable::NllTravelTimeTable(const std::string &type,
                                       const std::string &model)
    : TravelTimeTable(type, model)
{
  static const std::regex regex(";", std::regex::optimize);
  std::vector<std::string> tokens(splitString(model, regex));
  if (tokens.size() != 3 && tokens.size() != 4)
  {
    string msg =
        stringify("Error while initialzing NLL grids: invalid table model (%s)",
                  model.c_str());
    throw runtime_error(msg.c_str());
  }
  _velGridPath   = tokens[0];
  _timeGridPath  = tokens[1];
  _angleGridPath = tokens[2];
  _swapBytes     = false;
  if (tokens.size() > 3 && tokens[3] == "swapBytes")
  {
    _swapBytes = true;
  }
}

void NllTravelTimeTable::compute(double eventLat,
                                 double eventLon,
                                 double eventDepth,
                                 const Catalog::Station &station,
                                 const std::string &phaseType,
                                 double &travelTime)
{
  string timeGId =
      "timeGrid:" + Grid::filePath(_timeGridPath, station, phaseType);
  auto timeIt = _timeGrids.find(timeGId);
  if (timeIt == _timeGrids.end())
  {
    // Check if we have already excluded the grid because we couldn't load it
    if (_unloadableGrids.find(timeGId) != _unloadableGrids.end())
    {
      string msg = stringify("Time grid (%s) not avaliable", timeGId.c_str());
      throw runtime_error(msg.c_str());
    }

    try
    {
      TimeGridPtr tg =
          new TimeGrid(_timeGridPath, station, phaseType, _swapBytes);
      timeIt =
          _timeGrids.insert(std::pair<string, TimeGridPtr>(timeGId, tg)).first;
    }
    catch (exception &e)
    {
      _unloadableGrids.insert(timeGId);
      throw runtime_error(e.what());
    }
  }

  TimeGridPtr timeGrid = timeIt->second;
  travelTime           = timeGrid->getTime(eventLat, eventLon, eventDepth);
}

void NllTravelTimeTable::compute(double eventLat,
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

  string velGId = "velGrid:" + Grid::filePath(_velGridPath, station, phaseType);
  auto velIt    = _velGrids.find(velGId);
  if (velIt == _velGrids.end())
  {
    // Check if we have already excluded the grid because we couldn't load it
    if (_unloadableGrids.find(velGId) != _unloadableGrids.end())
    {
      string msg = stringify("Vel grid (%s) not avaliable", velGId.c_str());
      throw runtime_error(msg.c_str());
    }

    try
    {
      VelGridPtr vg = new VelGrid(_velGridPath, station, phaseType, _swapBytes);
      velIt = _velGrids.insert(std::pair<string, VelGridPtr>(velGId, vg)).first;
    }
    catch (exception &e)
    {
      _unloadableGrids.insert(velGId);
      throw runtime_error(e.what());
    }
  }

  // set velocityAtSrc
  VelGridPtr velGrid = velIt->second;
  velocityAtSrc      = velGrid->getVel(eventLat, eventLon, eventDepth);

  string angleGId =
      "angleGrid:" + Grid::filePath(_angleGridPath, station, phaseType);
  auto angleIt = _angleGrids.find(angleGId);
  if (angleIt == _angleGrids.end())
  {
    // Check if we have already excluded the grid because we couldn't load it
    if (_unloadableGrids.find(angleGId) == _unloadableGrids.end())
    {
      try
      {
        AngleGridPtr ag =
            new AngleGrid(_angleGridPath, station, phaseType, _swapBytes);
        angleIt =
            _angleGrids.insert(std::pair<string, AngleGridPtr>(angleGId, ag))
                .first;
      }
      catch (exception &e)
      {
        _unloadableGrids.insert(angleGId);
        SEISCOMP_WARNING(
            "Cannot load angle grid file: using approximated angles (%s)",
            e.what());
      }
    }
  }

  // set takeOffAngles
  takeOffAngleAzim = std::nan("");
  takeOffAngleDip  = std::nan("");
  if (angleIt != _angleGrids.end())
  {
    AngleGridPtr angleGrid = angleIt->second;
    try
    {
      angleGrid->getAngles(eventLat, eventLon, eventDepth, takeOffAngleAzim,
                           takeOffAngleDip);
    }
    catch (exception &e)
    {
      SEISCOMP_WARNING(
          "Error reading angle grid file: using approximated angles (%s)",
          e.what());
    }
  }
  // approximate angles if not already provided by the grid
  computeApproximatedTakeOfAngles(
      eventLat, eventLon, eventDepth, station, phaseType,
      std::isfinite(takeOffAngleAzim) ? nullptr : &takeOffAngleAzim,
      std::isfinite(takeOffAngleDip) ? nullptr : &takeOffAngleDip);
}

std::string Grid::filePath(const std::string &basePath,
                           const Catalog::Station &station,
                           const std::string &phaseType)
{
  static const std::regex reNet("NETWORK", std::regex::optimize);
  static const std::regex reSta("STATION", std::regex::optimize);
  static const std::regex reLoc("LOCATION", std::regex::optimize);
  static const std::regex rePha("PHASE", std::regex::optimize);

  string out = std::regex_replace(basePath, reNet, station.networkCode);
  out        = std::regex_replace(out, reSta, station.stationCode);
  out        = std::regex_replace(out, reLoc, station.locationCode);
  return std::regex_replace(out, rePha, phaseType);
}

Grid::Grid(Type gridType,
           const std::string &basePath,
           const Catalog::Station &station,
           const std::string &phaseType,
           bool swapBytes)
    : info(parse(filePath(basePath, station, phaseType), gridType, swapBytes))
{

  if (!Util::fileExists(info.bufFilePath))
  {
    string msg =
        stringify("Cannot find grid data file %s", info.bufFilePath.c_str());
    throw runtime_error(msg.c_str());
  }
}

Grid::Info
Grid::parse(const std::string &baseFilePath, Type gridType, bool swapBytes)
{
  Info info;

  info.hdrFilePath = baseFilePath + ".hdr";
  info.bufFilePath = baseFilePath + ".buf";
  info.gridType    = gridType;
  info.swapBytes   = swapBytes; // is there a way to auto-detect this setting?

  // read file one line a time
  int parsedLines = 0;
  ifstream in(info.hdrFilePath);
  while (!in.eof())
  {
    string line;
    std::getline(in, line);
    if (in.bad() || in.fail()) break;

    // split line on blanks
    static const std::regex regex(R"([\s]+)", std::regex::optimize);
    std::vector<std::string> tokens(splitString(line, regex));

    // remove the first empty element if the line start with spaces
    if (!tokens.empty() && tokens[0].empty()) tokens.erase(tokens.begin());

    // skip empty lines
    if (tokens.empty()) continue;

    if (tokens.size() == 10 || tokens.size() == 11) // should be line 1
    {
      info.numx      = std::stoull(tokens[0]);
      info.numy      = std::stoull(tokens[1]);
      info.numz      = std::stoull(tokens[2]);
      info.origx     = std::stod(tokens[3]);
      info.origy     = std::stod(tokens[4]);
      info.origz     = std::stod(tokens[5]);
      info.dx        = std::stod(tokens[6]);
      info.dy        = std::stod(tokens[7]);
      info.dz        = std::stod(tokens[8]);
      info.type      = tokens[9];
      info.useDouble = (tokens.size() == 11) ? (tokens[10] == "DOUBLE") : false;
      ++parsedLines;
    }
    else if (tokens.size() == 4 &&
             (gridType == Type::time || gridType == Type::angle))
    {
      // this should be line 2
      info.label = tokens[0];
      info.srcex = std::stod(tokens[1]);
      info.srcey = std::stod(tokens[2]);
      info.srcez = std::stod(tokens[3]);
      ++parsedLines;
    }
    else if (tokens[0] == "TRANSFORM" || tokens[0] == "TRANS")
    {
      info.transform = std::unique_ptr<Transform>(new Transform(tokens));
      ++parsedLines;
    }
  }

  if (((gridType == Type::time || gridType == Type::angle) &&
       parsedLines != 3) ||
      (gridType == Type::velocity && parsedLines != 2))
  {
    string msg =
        stringify("Cannot load grid header file %s", info.hdrFilePath.c_str());
    throw runtime_error(msg.c_str());
  }

  // make sure that dx for 2D grids is non-zero
  if (info.numx == 1) info.dx = 1.0;

  if (info.useDouble && info.swapBytes)
  {
    throw runtime_error(
        "Grid files with DOUBLE values and byte swapping are not supported");
  }

  return info;
}

bool Grid::isLocationInside(double xloc, double yloc, double zloc) const
{
  if (xloc < info.origx || xloc > info.origx + (info.numx - 1) * info.dx)
    return false;
  if (yloc < info.origy || yloc > info.origy + (info.numy - 1) * info.dy)
    return false;
  if (zloc < info.origz || zloc > info.origz + (info.numz - 1) * info.dz)
    return false;
  return true;
}

bool Grid::isIndexInside(unsigned long long ix,
                         unsigned long long iy,
                         unsigned long long iz) const
{
  return !(ix < 0 || ix >= info.numx || iy < 0 || iy >= info.numy || iz < 0 ||
           iz >= info.numz);
}

/*
 * This is the function that reads the data from the grid files
 * Since the grid files can vary between hundres of MB to hundreds
 * of GB (even TB) we do not load all the grid files in memory, but
 * we read the data directly from the files, for each value!
 * That might raise concerns on the performance, due to the high I/O
 * involved. However the OS + Filesystem put lots of efforts on
 * caching mechanisms. For exampe when a simple float is read from
 * disk a full page (e.g. 4kB) will be read from disk and stored
 * in memory, so that subsequents reads at close offsets will be
 * fulfilled without additional I/O being involved.
 * In general it is hard to do better than OS+filesystem caching
 * and it might not be worth trying so until real-world scenarios
 * (e.g. grids stored on a network folder) prove an additional
 * caching layer to be necessary
 */
template <class GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE Grid::getValueAtIndex(unsigned long long ix,
                                      unsigned long long iy,
                                      unsigned long long iz)
{
  if (!isIndexInside(ix, iy, iz))
  {
    throw runtime_error("Requested index is out of grid boundaries");
  }

  if (!_bufReader.is_open())
  {
    _bufReader.open(info.bufFilePath, std::ios::binary | std::ios::in);
    _bufReader.exceptions(std::ios::badbit | std::ios::failbit);
  }

  unsigned long pos = sizeof(GRID_FLOAT_TYPE) *
                      (ix * info.numy * info.numz + iy * info.numz + iz);

  GRID_FLOAT_TYPE value;
  try
  {
    _bufReader.seekg(pos);
    _bufReader.read(reinterpret_cast<char *>(&value), sizeof(value));
  }
  catch (exception &e)
  {
    _bufReader.close();
    string msg = stringify("Error while reading grid file %s (%s)",
                           info.bufFilePath.c_str(), e.what());
    throw runtime_error(msg.c_str());
  }

  if (info.swapBytes && sizeof(GRID_FLOAT_TYPE) == 4)
  {
    char *unswappedFloat = (char *)&value;
    char swappedFloat[4];
    swappedFloat[0] = unswappedFloat[3];
    swappedFloat[1] = unswappedFloat[2];
    swappedFloat[2] = unswappedFloat[1];
    swappedFloat[3] = unswappedFloat[0];
    std::memcpy(&value, swappedFloat, 4);
  }
  return value;
}

template <class GRID_FLOAT_TYPE>
void Grid::getValuesAt3DLocation(double xloc,
                                 double yloc,
                                 double zloc,
                                 double &xdiff,
                                 double &ydiff,
                                 double &zdiff,
                                 GRID_FLOAT_TYPE &vval000,
                                 GRID_FLOAT_TYPE &vval001,
                                 GRID_FLOAT_TYPE &vval010,
                                 GRID_FLOAT_TYPE &vval011,
                                 GRID_FLOAT_TYPE &vval100,
                                 GRID_FLOAT_TYPE &vval101,
                                 GRID_FLOAT_TYPE &vval110,
                                 GRID_FLOAT_TYPE &vval111)
{
  if (!isLocationInside(xloc, yloc, zloc))
  {
    string msg = stringify("Requested location is out of grid boundaries "
                           "(xloc %.2f yloc %.2f zloc %.2f - grid %s "
                           "origx %.3f origy %.3f origz %.3f "
                           "dx %.2f dy %.2f dz %.2f "
                           "numx %u numy %u numz %u)",
                           xloc, yloc, zloc, info.hdrFilePath.c_str(),
                           info.origx, info.origy, info.origz, info.dx, info.dy,
                           info.dz, info.numx, info.numy, info.numz);
    throw runtime_error(msg.c_str());
  }

  /* calculate grid locations at the vertex of the cube containing the point */

  double xoff = (xloc - info.origx) / info.dx;
  double yoff = (yloc - info.origy) / info.dy;
  double zoff = (zloc - info.origz) / info.dz;

  unsigned long long ix0 = (unsigned long long)xoff;
  unsigned long long iy0 = (unsigned long long)yoff;
  unsigned long long iz0 = (unsigned long long)zoff;

  if (ix0 == info.numx - 1) ix0--;
  if (iy0 == info.numy - 1) iy0--;
  if (iz0 == info.numz - 1) iz0--;

  unsigned long long ix1 = ix0 + 1;
  unsigned long long iy1 = iy0 + 1;
  unsigned long long iz1 = iz0 + 1;

  xdiff = xoff - ix0;
  ydiff = yoff - iy0;
  zdiff = zoff - iz0;

  /* read vertex values from grid file */

  vval000 = getValueAtIndex<GRID_FLOAT_TYPE>(ix0, iy0, iz0);
  vval001 = getValueAtIndex<GRID_FLOAT_TYPE>(ix0, iy0, iz1);
  vval010 = getValueAtIndex<GRID_FLOAT_TYPE>(ix0, iy1, iz0);
  vval011 = getValueAtIndex<GRID_FLOAT_TYPE>(ix0, iy1, iz1);
  vval100 = getValueAtIndex<GRID_FLOAT_TYPE>(ix1, iy0, iz0);
  vval101 = getValueAtIndex<GRID_FLOAT_TYPE>(ix1, iy0, iz1);
  vval110 = getValueAtIndex<GRID_FLOAT_TYPE>(ix1, iy1, iz0);
  vval111 = getValueAtIndex<GRID_FLOAT_TYPE>(ix1, iy1, iz1);
}

template <class GRID_FLOAT_TYPE>
void Grid::getValuesAt2DLocation(double yloc,
                                 double zloc,
                                 double &ydiff,
                                 double &zdiff,
                                 GRID_FLOAT_TYPE &vval00,
                                 GRID_FLOAT_TYPE &vval01,
                                 GRID_FLOAT_TYPE &vval10,
                                 GRID_FLOAT_TYPE &vval11)

{
  double xloc = info.origx;

  if (!isLocationInside(xloc, yloc, zloc))
  {
    string msg = stringify("Requested location is out of grid boundaries "
                           "(xloc %.2f yloc %.2f zloc %.2f - grid %s "
                           "origx %.3f origy %.3f origz %.3f "
                           "dx %.2f dy %.2f dz %.2f "
                           "numx %u numy %u numz %u)",
                           xloc, yloc, zloc, info.hdrFilePath.c_str(),
                           info.origx, info.origy, info.origz, info.dx, info.dy,
                           info.dz, info.numx, info.numy, info.numz);
    throw runtime_error(msg.c_str());
  }

  /* calculate grid locations at the face of the cube containing the point */

  double yoff = (yloc - info.origy) / info.dy;
  double zoff = (zloc - info.origz) / info.dz;

  unsigned long long ix0 = 0;
  unsigned long long iy0 = (unsigned long long)yoff;
  unsigned long long iz0 = (unsigned long long)zoff;

  if (iy0 == info.numy - 1) iy0--;
  if (iz0 == info.numz - 1) iz0--;

  unsigned long long iy1 = iy0 + 1;
  unsigned long long iz1 = iz0 + 1;

  ydiff = yoff - iy0;
  zdiff = zoff - iz0;

  /* read vertex values from grid file */

  vval00 = getValueAtIndex<GRID_FLOAT_TYPE>(ix0, iy0, iz0);
  vval01 = getValueAtIndex<GRID_FLOAT_TYPE>(ix0, iy0, iz1);
  vval10 = getValueAtIndex<GRID_FLOAT_TYPE>(ix0, iy1, iz0);
  vval11 = getValueAtIndex<GRID_FLOAT_TYPE>(ix0, iy1, iz1);
}

template <class GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE Grid::getValue(
    double lat,
    double lon,
    double depth,
    const typename Interpolate3D<GRID_FLOAT_TYPE>::Type &interpolateValues3D,
    const typename Interpolate2D<GRID_FLOAT_TYPE>::Type &interpolateValues2D)
{
  if (is3D()) // 3D grid
  {
    return getValue3D<GRID_FLOAT_TYPE>(lat, lon, depth, interpolateValues3D);
  }
  else // 2D grid (1D model)
  {
    return getValue2D<GRID_FLOAT_TYPE>(lat, lon, depth, interpolateValues2D);
  }
}

template <class GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE Grid::getValue3D(
    double lat,
    double lon,
    double depth,
    const typename Interpolate3D<GRID_FLOAT_TYPE>::Type &interpolateValues3D)
{
  double xLoc, yLoc;
  info.transform->fromLatLon(lat, lon, xLoc, yLoc);

  double xdiff, ydiff, zdiff;
  GRID_FLOAT_TYPE vval000, vval001, vval010, vval011, vval100, vval101, vval110,
      vval111;

  getValuesAt3DLocation<GRID_FLOAT_TYPE>(xLoc, yLoc, depth, xdiff, ydiff, zdiff,
                                         vval000, vval001, vval010, vval011,
                                         vval100, vval101, vval110, vval111);

  return interpolateValues3D(xdiff, ydiff, zdiff, vval000, vval001, vval010,
                             vval011, vval100, vval101, vval110, vval111);
}

template <class GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE Grid::getValue2D(
    double lat,
    double lon,
    double depth,
    const typename Interpolate2D<GRID_FLOAT_TYPE>::Type &interpolateValues2D)
{
  double xLoc, yLoc;
  info.transform->fromLatLon(lat, lon, xLoc, yLoc);

  double ydiff, zdiff;
  GRID_FLOAT_TYPE vval00, vval01, vval10, vval11;

  double dist;

  // This check means the code shuould move to the sub-classes, but I don't
  // like all the code duplication that the move would require
  if (info.gridType == Type::velocity)
  {
    // VEL grid headers don't have _srce informationm, the velocity depends on
    // the depth only, not on the distance to the station
    dist = info.origy;
  }
  else
  {
    dist = info.transform->distance(xLoc, yLoc, info.srcex, info.srcey);
  }

  getValuesAt2DLocation<GRID_FLOAT_TYPE>(dist, depth, ydiff, zdiff, vval00,
                                         vval01, vval10, vval11);

  return interpolateValues2D(ydiff, zdiff, vval00, vval01, vval10, vval11);
}

TimeGrid::TimeGrid(const std::string &basePath,
                   const Catalog::Station &station,
                   const std::string &phaseType,
                   bool swapBytes)
    : Grid(Type::time, basePath, station, phaseType, swapBytes)
{
  if (info.type != "TIME" && info.type != "TIME2D")
  {
    string msg = stringify("Unrecognized time grid type %s (%s)",
                           info.type.c_str(), info.hdrFilePath.c_str());
    throw runtime_error(msg.c_str());
  }
}

double TimeGrid::getTime(double lat, double lon, double depth)
{
  if (info.useDouble)
  {
    return getValue<double>(lat, lon, depth, interpolateValues3D<double>,
                            interpolateValues2D<double>);
  }
  else
  {
    return getValue<float>(lat, lon, depth, interpolateValues3D<float>,
                           interpolateValues2D<float>);
  }
}

template <class GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE TimeGrid::interpolateValues3D(double xdiff,
                                              double ydiff,
                                              double zdiff,
                                              GRID_FLOAT_TYPE vval000,
                                              GRID_FLOAT_TYPE vval001,
                                              GRID_FLOAT_TYPE vval010,
                                              GRID_FLOAT_TYPE vval011,
                                              GRID_FLOAT_TYPE vval100,
                                              GRID_FLOAT_TYPE vval101,
                                              GRID_FLOAT_TYPE vval110,
                                              GRID_FLOAT_TYPE vval111)
{

  if (vval000 < 0.0 || vval010 < 0.0 || vval100 < 0.0 || vval110 < 0.0 ||
      vval001 < 0.0 || vval011 < 0.0 || vval101 < 0.0 || vval111 < 0.0)
  {
    throw runtime_error("Negative times found in the grid file");
  }
  return interpolateCubeLagrange(xdiff, ydiff, zdiff, vval000, vval001, vval010,
                                 vval011, vval100, vval101, vval110, vval111);
}

template <class GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE TimeGrid::interpolateValues2D(double xdiff,
                                              double zdiff,
                                              GRID_FLOAT_TYPE vval00,
                                              GRID_FLOAT_TYPE vval01,
                                              GRID_FLOAT_TYPE vval10,
                                              GRID_FLOAT_TYPE vval11)
{
  if (vval00 < 0.0 || vval01 < 0.0 || vval10 < 0.0 || vval11 < 0.0)
  {
    throw runtime_error("Negative times found in the grid file");
  }
  return interpolateSquareLagrange(xdiff, zdiff, vval00, vval01, vval10,
                                   vval11);
}

AngleGrid::AngleGrid(const std::string &basePath,
                     const Catalog::Station &station,
                     const std::string &phaseType,
                     bool swapBytes)
    : Grid(Type::angle, basePath, station, phaseType, swapBytes)
{
  if (info.type != "ANGLE" && info.type != "ANGLE2D")
  {
    string msg =
        stringify("Unrecognized angle grid type %s", info.type.c_str());
    throw runtime_error(msg.c_str());
  }

  //
  // TakeOffAngles uses bit-fields, but their memory arrangment
  // (how bits are packed together) is not standard. So let's make
  // sure the compiler did a good job
  //
  if (sizeof(TakeOffAngles) != sizeof(float))
  {
    throw runtime_error(
        "Internal error: sizeof(TakeOffAngles) != sizeof(float)");
  }

  unsigned short quality  = 0x06;
  unsigned short dip      = 0xF9;
  unsigned short azimmuth = 0xCC33;

  TakeOffAngles test{quality, dip, azimmuth};

  // This is how we want the bit-fields mapped in memory
  // bits 0-3   : quality
  // bits 4-15  : dip
  // bits 16-31 : azimuth
  unsigned char reference[4];
  reference[0] = ((dip << 4) | (quality & 0xF)) & 0x00FF;
  reference[1] = (((dip << 4) | (quality & 0xF)) & 0xFF00) >> 8;
  reference[2] = azimmuth & 0xFF;
  reference[3] = (azimmuth & 0xFF00) >> 8;

  if (std::memcmp(&test, reference, 4) != 0)
  {
    throw runtime_error(
        "Internal error: TakeOffAngles memory mapping is not ok");
  }

  // assumes angle files store floats only (never double)
  if (info.useDouble)
  {
    throw runtime_error("Angle grid files with DOUBLE values are not "
                        "supported, only FLOAT allowed");
  }
}

void AngleGrid::getAngles(
    double lat, double lon, double depth, double &azim, double &dip)
{
  TakeOffAngles angles = getValue<TakeOffAngles>(
      lat, lon, depth, interpolateValues3D<TakeOffAngles>,
      interpolateValues2D<TakeOffAngles>);

  if (angles.quality < QUALITY_CUTOFF)
  {
    // this is not an error, the code handles nans
    azim = std::nan("");
    dip  = std::nan("");
    return;
  }

  if (is3D())
  {
    azim = angles.azimuth / 10.0; // tenths of degree -> degree
    azim = info.transform->toLatLonAngle(azim);
    azim = deg2rad(azim);
  }
  else
  {
    azim = std::nan("");
  }
  dip = (angles.dip / 10.0);
  dip = deg2rad(dip);
}

template <class GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE AngleGrid::interpolateValues3D(double xdiff,
                                               double ydiff,
                                               double zdiff,
                                               GRID_FLOAT_TYPE vval000,
                                               GRID_FLOAT_TYPE vval001,
                                               GRID_FLOAT_TYPE vval010,
                                               GRID_FLOAT_TYPE vval011,
                                               GRID_FLOAT_TYPE vval100,
                                               GRID_FLOAT_TYPE vval101,
                                               GRID_FLOAT_TYPE vval110,
                                               GRID_FLOAT_TYPE vval111)
{
  return interpolateCubeAngles(xdiff, ydiff, zdiff, vval000, vval001, vval010,
                               vval011, vval100, vval101, vval110, vval111,
                               QUALITY_CUTOFF);
}

template <class GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE AngleGrid::interpolateValues2D(double xdiff,
                                               double zdiff,
                                               GRID_FLOAT_TYPE vval00,
                                               GRID_FLOAT_TYPE vval01,
                                               GRID_FLOAT_TYPE vval10,
                                               GRID_FLOAT_TYPE vval11)
{
  return interpolateCubeAngles(0, xdiff, zdiff, vval00, vval01, vval10, vval11,
                               vval00, vval01, vval10, vval11, QUALITY_CUTOFF);
}

VelGrid::VelGrid(const std::string &basePath,
                 const Catalog::Station &station,
                 const std::string &phaseType,
                 bool swapBytes)
    : Grid(Type::velocity, basePath, station, phaseType, swapBytes)
{
  if (info.numx < 2)
  {
    string msg = stringify(
        "Velocity grid must have xNum greater than 2, found %llu (%s)",
        info.numx, info.hdrFilePath.c_str());
    throw runtime_error(msg.c_str());
  }

  if (info.type == "VELOCITY_METERS")
    convertUnits = [](double vel) -> double { return vel / 1000.0; };
  else if (info.type == "SLOWNESS")
    convertUnits = [](double vel) -> double { return 1.0 / vel; };
  else if (info.type == "SLOW_LEN")
  {
    convertUnits = [this](double vel) -> double {
      return 1.0 / (vel / info.dy);
    };
  }
  else if (info.type == "VEL2")
    convertUnits = [](double vel) -> double { return std::sqrt(vel); };
  else if (info.type == "SLOW2")
    convertUnits = [](double vel) -> double { return std::sqrt(1.0 / vel); };
  else if (info.type == "SLOW2_METERS")
    convertUnits = [](double vel) -> double {
      return std::sqrt(1.0 / vel) / 1000.0;
    };
  else if (info.type == "VELOCITY")
    convertUnits = [](double vel) -> double { return vel; };
  else
  {
    string msg =
        stringify("Unrecognized velocity grid type %s", info.type.c_str());
    throw runtime_error(msg.c_str());
  }
}

double VelGrid::getVel(double lat, double lon, double depth)
{
  double velAtSrc;
  if (info.useDouble)
  {
    velAtSrc = getValue<double>(lat, lon, depth, interpolateValues3D<double>,
                                interpolateValues2D<double>);
  }
  else
  {
    velAtSrc = getValue<float>(lat, lon, depth, interpolateValues3D<float>,
                               interpolateValues2D<float>);
  }
  // velocity -> [km/sec]
  return convertUnits(velAtSrc);
}

template <class GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE VelGrid::interpolateValues3D(double xdiff,
                                             double ydiff,
                                             double zdiff,
                                             GRID_FLOAT_TYPE vval000,
                                             GRID_FLOAT_TYPE vval001,
                                             GRID_FLOAT_TYPE vval010,
                                             GRID_FLOAT_TYPE vval011,
                                             GRID_FLOAT_TYPE vval100,
                                             GRID_FLOAT_TYPE vval101,
                                             GRID_FLOAT_TYPE vval110,
                                             GRID_FLOAT_TYPE vval111)
{
  if (vval000 < 0.0 || vval010 < 0.0 || vval100 < 0.0 || vval110 < 0.0 ||
      vval001 < 0.0 || vval011 < 0.0 || vval101 < 0.0 || vval111 < 0.0)
  {
    throw runtime_error("Negative velocities found in the grid file");
  }
  return interpolateCubeLagrange(xdiff, ydiff, zdiff, vval000, vval001, vval010,
                                 vval011, vval100, vval101, vval110, vval111);
}

template <class GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE VelGrid::interpolateValues2D(double xdiff,
                                             double zdiff,
                                             GRID_FLOAT_TYPE vval00,
                                             GRID_FLOAT_TYPE vval01,
                                             GRID_FLOAT_TYPE vval10,
                                             GRID_FLOAT_TYPE vval11)
{
  if (vval00 < 0.0 || vval01 < 0.0 || vval10 < 0.0 || vval11 < 0.0)
  {
    throw runtime_error("Negative velocities found in the grid file");
  }
  return interpolateSquareLagrange(xdiff, zdiff, vval00, vval01, vval10,
                                   vval11);
}

Transform::Transform(const std::vector<std::string> &tokens)
    : info(parse(tokens))
{}

Transform::Info Transform::parse(const std::vector<string> &tokens)
{
  Info info;

  if (tokens[0] != "TRANSFORM" && tokens[0] != "TRANS")
  {
    throw runtime_error("Malformed transform line");
  }

  info.type = tokens[1];

  if (info.type == "GLOBAL" || info.type == "NONE")
  {
    info.angle  = 0.0;
    info.cosang = std::cos(info.angle);
    info.sinang = std::sin(info.angle);
  }
  else if (info.type == "SIMPLE" || info.type == "SDC")
  {
    info.orig_lat  = std::stod(tokens[3]);
    info.orig_long = std::stod(tokens[5]);
    info.rot       = std::stod(tokens[7]);
    info.angle     = -deg2rad(info.rot);
    info.cosang    = std::cos(info.angle);
    info.sinang    = std::sin(info.angle);

    if (info.orig_lat > 90 || info.orig_lat < -90)
    {
      throw runtime_error("Origin latitude must be in range -90,90");
    }
    if (info.orig_long > 180 || info.orig_long < -180)
    {
      throw runtime_error("Origin longitude must be in range -180,180");
    }
    if (info.rot > 360 || info.rot < -360)
    {
      throw runtime_error("Rotation must be in range -360,360");
    }

    if (info.type == "SDC")
    {
      //  conversion factor for latitude
      double dlt1 =
          std::atan(MAP_TRANS_SDC_DRLT * std::tan(deg2rad(info.orig_lat)));
      double dlt2    = std::atan(MAP_TRANS_SDC_DRLT *
                              std::tan(deg2rad(info.orig_lat + 1.0)));
      double del     = dlt2 - dlt1;
      double r       = ERAD * (1.0 - square(std::sin(dlt1)) * FLATTENING);
      info.sdc_xltkm = r * del;
      //  conversion factor for longitude
      del            = std::acos(1.0 -
                      (1.0 - std::cos(deg2rad(1))) * square(std::cos(dlt1)));
      double bc      = r * del;
      info.sdc_xlnkm = bc / std::cos(dlt1);
    }
  }
  else
  {
    string msg = stringify("Unsupported transform %s", info.type.c_str());
    throw runtime_error(msg.c_str());
  }
  return info;
}

void Transform::fromLatLon(double lat,
                           double lon,
                           double &xLoc,
                           double &yLoc) const
{

  if (info.type == "GLOBAL" || info.type == "NONE")
  {
    xLoc = lon;
    yLoc = lat;
  }
  else if (info.type == "SIMPLE")
  {
    double xtemp = lon - info.orig_long;
    if (xtemp > 180.0) xtemp -= 360.0;
    if (xtemp < -180.0) xtemp += 360.0;
    xtemp        = xtemp * c111 * std::cos(deg2rad(lat));
    double ytemp = (lat - info.orig_lat) * c111;
    xLoc         = xtemp * info.cosang - ytemp * info.sinang;
    yLoc         = ytemp * info.cosang + xtemp * info.sinang;
  }
  else if (info.type == "SDC")
  {
    double xtemp = lon - info.orig_long;
    if (xtemp > 180.0) xtemp -= 360.0;
    if (xtemp < -180.0) xtemp += 360.0;
    double ytemp = lat - info.orig_lat;

    double xlt1 = std::atan(MAP_TRANS_SDC_DRLT *
                            std::tan(deg2rad(lat + info.orig_lat) / 2.0));
    xtemp       = xtemp * info.sdc_xlnkm * std::cos(xlt1);
    ytemp       = ytemp * info.sdc_xltkm;

    xLoc = xtemp * info.cosang - ytemp * info.sinang;
    yLoc = ytemp * info.cosang + xtemp * info.sinang;
  }
  else // this never happens
  {
    string msg = stringify("Unsupported transform %s", info.type.c_str());
    throw runtime_error(msg.c_str());
  }
}

void Transform::toLatLon(double xLoc,
                         double yLoc,
                         double &lat,
                         double &lon) const
{
  if (info.type == "GLOBAL" || info.type == "NONE")
  {
    lat = yLoc;
    lon = xLoc;
  }
  else if (info.type == "SIMPLE")
  {
    double xtemp = xLoc * info.cosang + yLoc * info.sinang;
    double ytemp = yLoc * info.cosang - xLoc * info.sinang;
    lat          = info.orig_lat + ytemp / c111;
    lon          = info.orig_long + xtemp / (c111 * std::cos(deg2rad(lat)));
    if (lon < -180.0)
      lon += 360.0;
    else if (lon > 180.0)
      lon -= 360.0;
  }
  else if (info.type == "SDC")
  {
    double xtemp = xLoc * info.cosang + yLoc * info.sinang;
    double ytemp = yLoc * info.cosang - xLoc * info.sinang;
    ytemp        = ytemp / info.sdc_xltkm;
    lat          = info.orig_lat + ytemp;
    double xlt1  = std::atan(MAP_TRANS_SDC_DRLT *
                            std::tan(deg2rad(lat + info.orig_lat) / 2.0));
    xtemp        = xtemp / (info.sdc_xlnkm * std::cos(xlt1));
    lon          = info.orig_long + xtemp;
    if (lon < -180.0)
      lon += 360.0;
    else if (lon > 180.0)
      lon -= 360.0;
  }
  else // this never happens
  {
    string msg = stringify("Unsupported transform %s", info.type.c_str());
    throw runtime_error(msg.c_str());
  }
}

double Transform::fromLatLonAngle(double latlonAngle) const
{
  if (info.type == "SIMPLE" || info.type == "SDC")
  {
    double angle = latlonAngle + info.rot;
    if (angle < 0.0)
      angle += 360.0;
    else if (angle > 360.0)
      angle -= 360.0;
    return angle;
  }
  else if (info.type == "GLOBAL" || info.type == "NONE")
  {
    return latlonAngle;
  }
  else // this never happens
  {
    string msg = stringify("Unsupported transform %s", info.type.c_str());
    throw runtime_error(msg.c_str());
  }
}

double Transform::toLatLonAngle(double rectAngle) const
{
  if (info.type == "SIMPLE" || info.type == "SDC")
  {
    double angle = rectAngle - info.rot;
    if (angle < 0.0)
      angle += 360.0;
    else if (angle > 360.0)
      angle -= 360.0;
    return (angle);
  }
  else if (info.type == "GLOBAL" || info.type == "NONE")
  {
    return rectAngle;
  }
  else // this never happens
  {
    string msg = stringify("Unsupported transform %s", info.type.c_str());
    throw runtime_error(msg.c_str());
  }
}

double Transform::distance(double xLoc1,
                           double yLoc1,
                           double xLoc2,
                           double yLoc2) const
{
  if (info.type == "GLOBAL")
  {
    return computeDistance(yLoc1, xLoc1, yLoc2, xLoc2);
  }
  else
  {
    return std::sqrt(square(xLoc2 - xLoc1) + square(yLoc2 - yLoc1));
  }
}

double Transform::distance(double xLoc1,
                           double yLoc1,
                           double zLoc1,
                           double xLoc2,
                           double yLoc2,
                           double zLoc2) const
{
  if (info.type == "GLOBAL")
  {
    return computeDistance(yLoc1, xLoc1, zLoc1, yLoc2, xLoc2, zLoc2);
  }
  else
  {
    return std::sqrt(square(xLoc2 - xLoc1) + square(yLoc2 - yLoc1) +
                     square(zLoc2 - zLoc1));
  }
}

} // namespace NLL
} // namespace HDD
} // namespace Seiscomp

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

#ifndef __HDD_NLLGRID_H__
#define __HDD_NLLGRID_H__

#include "catalog.h"

#include <cmath>
#include <fstream>
#include <functional>
#include <memory>
#include <unordered_map>
#include <vector>

namespace HDD {
namespace NLL {

class Transform
{
public:
  static std::unique_ptr<Transform>
  parse(const std::vector<std::string> &tokens);

  Transform(const std::string &type,
            double orgLat,
            double _orgLong,
            double rot);
  virtual ~Transform() = default;

  Transform(const Transform &other)            = delete;
  Transform &operator=(const Transform &other) = delete;

  virtual void
  fromLatLon(double lat, double lon, double &xLoc, double &yLoc) const = 0;
  virtual void
  toLatLon(double xLoc, double yLoc, double &lat, double &lon) const = 0;

  virtual double fromLatLonAngle(double latlonAngle) const;
  virtual double toLatLonAngle(double rectAngle) const;

  virtual double
  distance(double xLoc1, double yLoc1, double xLoc2, double yLoc2) const;
  virtual double distance(double xLoc1,
                          double yLoc1,
                          double zLoc1,
                          double xLoc2,
                          double yLoc2,
                          double zLoc2) const;

protected:
  void rotate(double &xLoc, double &yLoc) const
  {
    double xtemp = xLoc;
    double ytemp = yLoc;
    xLoc         = xtemp * _cosang - ytemp * _sinang;
    yLoc         = ytemp * _cosang + xtemp * _sinang;
  }

  const std::string _type;
  const double _orgLat;
  const double _orgLong;
  const double _rot;
  const double _angle;
  const double _cosang;
  const double _sinang;
};

class Grid
{
public:
  enum class Type
  {
    time,
    angle,
    velocity
  };

  static std::string filePath(const std::string &basePath,
                              const Catalog::Station &station,
                              const std::string &phaseType);
  Grid(Type gridType,
       const std::string &basePath,
       const Catalog::Station &station,
       const std::string &phaseType,
       bool swapBytes);

  ~Grid() = default;

  bool isLocationInside(double xloc, double yloc, double zloc) const;
  bool isIndexInside(unsigned long long ix,
                     unsigned long long iy,
                     unsigned long long iz) const;

  template <typename GRID_FLOAT_TYPE> struct Interpolate2D
  {
    typedef std::function<GRID_FLOAT_TYPE(double,
                                          double,
                                          GRID_FLOAT_TYPE,
                                          GRID_FLOAT_TYPE,
                                          GRID_FLOAT_TYPE,
                                          GRID_FLOAT_TYPE)>
        Type;
  };

  template <typename GRID_FLOAT_TYPE> struct Interpolate3D
  {
    typedef std::function<GRID_FLOAT_TYPE(double,
                                          double,
                                          double,
                                          GRID_FLOAT_TYPE,
                                          GRID_FLOAT_TYPE,
                                          GRID_FLOAT_TYPE,
                                          GRID_FLOAT_TYPE,
                                          GRID_FLOAT_TYPE,
                                          GRID_FLOAT_TYPE,
                                          GRID_FLOAT_TYPE,
                                          GRID_FLOAT_TYPE)>
        Type;
  };

  template <typename GRID_FLOAT_TYPE>
  GRID_FLOAT_TYPE
  getValue3D(double lat,
             double lon,
             double depth,
             const typename Interpolate3D<GRID_FLOAT_TYPE>::Type &);

  template <typename GRID_FLOAT_TYPE>
  GRID_FLOAT_TYPE
  getValue2D(double lat,
             double lon,
             double depth,
             const typename Interpolate2D<GRID_FLOAT_TYPE>::Type &);

  template <typename GRID_FLOAT_TYPE>
  void getValuesAt3DLocation(double xloc,
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
                             GRID_FLOAT_TYPE &vval111);

  template <typename GRID_FLOAT_TYPE>
  void getValuesAt2DLocation(double yloc,
                             double zloc,
                             double &ydiff,
                             double &zdiff,
                             GRID_FLOAT_TYPE &vval00,
                             GRID_FLOAT_TYPE &vval01,
                             GRID_FLOAT_TYPE &vval10,
                             GRID_FLOAT_TYPE &vval11);

  template <typename GRID_FLOAT_TYPE>
  GRID_FLOAT_TYPE getValueAtIndex(unsigned long long ix,
                                  unsigned long long iy,
                                  unsigned long long iz);

  struct Info
  {
    std::string hdrFilePath;
    std::string bufFilePath;
    Type gridType;
    bool swapBytes; // should disk values bytes be swapped?

    unsigned long long numx, numy, numz;
    double origx, origy, origz; // km
    double dx, dy, dz;          // km
    std::string type;
    bool useDouble; // grid values stored as double instead of float
    std::string label;
    double srcex, srcey, srcez;
    std::unique_ptr<Transform> transform;
  };
  const Info info;

  static Info
  parse(const std::string &baseFilePath, Type gridType, bool swapBytes);

private:
  std::ifstream _bufReader;
};

class TimeGrid
{
public:
  TimeGrid(const std::string &basePath,
           const Catalog::Station &station,
           const std::string &phaseType,
           bool swapBytes);
  ~TimeGrid() = default;

  double getTime(double lat, double lon, double depth);

  bool is3D() const { return _grid.info.numx > 1; }

private:
  template <typename GRID_FLOAT_TYPE>
  static GRID_FLOAT_TYPE interpolateValues2D(double xdiff,
                                             double zdiff,
                                             GRID_FLOAT_TYPE vval00,
                                             GRID_FLOAT_TYPE vval01,
                                             GRID_FLOAT_TYPE vval10,
                                             GRID_FLOAT_TYPE vval11);
  template <typename GRID_FLOAT_TYPE>
  static GRID_FLOAT_TYPE interpolateValues3D(double xdiff,
                                             double ydiff,
                                             double zdiff,
                                             GRID_FLOAT_TYPE vval000,
                                             GRID_FLOAT_TYPE vval001,
                                             GRID_FLOAT_TYPE vval010,
                                             GRID_FLOAT_TYPE vval011,
                                             GRID_FLOAT_TYPE vval100,
                                             GRID_FLOAT_TYPE vval101,
                                             GRID_FLOAT_TYPE vval110,
                                             GRID_FLOAT_TYPE vval111);
  Grid _grid;
};

class AngleGrid
{
public:
  AngleGrid(const std::string &basePath,
            const Catalog::Station &station,
            const std::string &phaseType,
            bool swapBytes);
  ~AngleGrid() = default;

  void
  getAngles(double lat, double lon, double depth, double &azim, double &dip);

  bool is3D() const { return _grid.info.numx > 1; }

  struct TakeOffAngles
  {
    unsigned short quality : 4;  // 0 to 10
    unsigned short dip : 12;     // 0 (down) to 1800 (up) in tenths deg
    unsigned short azimuth : 16; // 0 to 3600 in tenths deg
  };

  static const unsigned QUALITY_CUTOFF = 5;

private:
  template <typename GRID_FLOAT_TYPE>
  static GRID_FLOAT_TYPE interpolateValues2D(double xdiff,
                                             double zdiff,
                                             GRID_FLOAT_TYPE vval00,
                                             GRID_FLOAT_TYPE vval01,
                                             GRID_FLOAT_TYPE vval10,
                                             GRID_FLOAT_TYPE vval11);
  template <typename GRID_FLOAT_TYPE>
  static GRID_FLOAT_TYPE interpolateValues3D(double xdiff,
                                             double ydiff,
                                             double zdiff,
                                             GRID_FLOAT_TYPE vval000,
                                             GRID_FLOAT_TYPE vval001,
                                             GRID_FLOAT_TYPE vval010,
                                             GRID_FLOAT_TYPE vval011,
                                             GRID_FLOAT_TYPE vval100,
                                             GRID_FLOAT_TYPE vval101,
                                             GRID_FLOAT_TYPE vval110,
                                             GRID_FLOAT_TYPE vval111);
  Grid _grid;
};

class VelGrid
{
public:
  VelGrid(const std::string &basePath,
          const Catalog::Station &station,
          const std::string &phaseType,
          bool swapBytes);

  ~VelGrid() = default;

  double getVel(double lat, double lon, double depth);

  bool is3D() const { return _grid.info.numx > 2; }

private:
  template <typename GRID_FLOAT_TYPE>
  static GRID_FLOAT_TYPE interpolateValues2D(double xdiff,
                                             double zdiff,
                                             GRID_FLOAT_TYPE vval00,
                                             GRID_FLOAT_TYPE vval01,
                                             GRID_FLOAT_TYPE vval10,
                                             GRID_FLOAT_TYPE vval11);
  template <typename GRID_FLOAT_TYPE>
  static GRID_FLOAT_TYPE interpolateValues3D(double xdiff,
                                             double ydiff,
                                             double zdiff,
                                             GRID_FLOAT_TYPE vval000,
                                             GRID_FLOAT_TYPE vval001,
                                             GRID_FLOAT_TYPE vval010,
                                             GRID_FLOAT_TYPE vval011,
                                             GRID_FLOAT_TYPE vval100,
                                             GRID_FLOAT_TYPE vval101,
                                             GRID_FLOAT_TYPE vval110,
                                             GRID_FLOAT_TYPE vval111);

  std::function<double(double)> convertUnits; // velocity -> km/sec

  Grid _grid;
};

} // namespace NLL
} // namespace HDD

#endif

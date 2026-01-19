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

#ifndef __HDD_NLL_GRID_H__
#define __HDD_NLL_GRID_H__

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

  Transform(const Transform &)            = delete;
  Transform &operator=(const Transform &) = delete;
  Transform(Transform &&)                 = delete;
  Transform &operator=(Transform &&)      = delete;

  std::string getType() const { return _type; }

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
  void rotate(double &xLoc, double &yLoc) const;
  void inverseRotate(double &xLoc, double &yLoc) const;

  const std::string _type;
  const double _orgLat;
  const double _orgLong;
  const double _rot;
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

  enum class OpenMode
  {
    LoadIntoMemory,
    Mmap,
    Pread
  };

  Grid(Type gridType, const std::string &filePath, bool swapBytes);
  ~Grid();

  // mmap: memory mapping offers better performance, but cerain file systems
  // might not work well, so it is optional
  void open(OpenMode mode = OpenMode::Mmap);
  bool isOpen();
  void close();

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
    bool swapBytes; // should disk values be byte swapped?

    unsigned long long numx, numy, numz;
    double origx, origy, origz; // km
    double dx, dy, dz;          // km
    std::string type;
    bool useDouble; // grid values stored as double instead of float
    std::string label;
    double srcex, srcey, srcez;
    std::unique_ptr<Transform> transform;

    double minX() const { return origx; }
    double maxX() const { return origx + (numx - 1) * dx; }
    double minY() const { return origy; }
    double maxY() const { return origy + (numy - 1) * dy; }
    double minZ() const { return origz; }
    double maxZ() const { return origz + (numz - 1) * dz; }
    bool isLocationInside(double xloc, double yloc, double zloc) const
    {
      return !(xloc < minX() || xloc > maxX() || yloc < minY() ||
               yloc > maxY() || zloc < minZ() || zloc > maxZ());
    }
    bool isIndexInside(unsigned long long ix,
                       unsigned long long iy,
                       unsigned long long iz) const
    {
      return ix < numx && iy < numy && iz < numz;
    }
  };
  const Info info;

  static Info
  parse(const std::string &baseFilePath, Type gridType, bool swapBytes);

private:
  OpenMode _mode;
  int _fd                  = -1;
  void *_mapped            = nullptr;
  std::size_t _mappedSize  = 0;
  std::vector<char> _bytes = {};
};

class TimeGrid
{
public:
  TimeGrid(const std::string &filePath, bool swapBytes);

  void open(Grid::OpenMode mode);
  bool isOpen();
  void close();

  double getTime(double lat, double lon, double depth);

  double getTimeAtIndex(unsigned long long ix,
                        unsigned long long iy,
                        unsigned long long iz);

  bool is3D() const { return _grid.info.numx > 1; }

  bool isGlobal() const { return _grid.info.transform->getType() == "GLOBAL"; }

  const Grid::Info &getInfo() const { return _grid.info; }

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
  AngleGrid(const std::string &filePath,
            bool swapBytes,
            unsigned quality_cutoff = 5);

  void open(Grid::OpenMode mode);
  bool isOpen();
  void close();

  void
  getAngles(double lat, double lon, double depth, double &azim, double &dip);

  void getAnglesAtIndex(unsigned long long ix,
                        unsigned long long iy,
                        unsigned long long iz,
                        double &azim,
                        double &dip);

  bool is3D() const { return _grid.info.numx > 1; }

  bool isGlobal() const { return _grid.info.transform->getType() == "GLOBAL"; }

  const Grid::Info &getInfo() const { return _grid.info; }

  struct TakeOffAngles
  {
    unsigned short quality : 4;  // 0 to 10
    unsigned short dip : 12;     // 0 (down) to 1800 (up) in tenths deg
    unsigned short azimuth : 16; // 0 to 3600 in tenths deg, clockwise
  };

private:
  template <typename GRID_FLOAT_TYPE>
  GRID_FLOAT_TYPE interpolateValues2D(double xdiff,
                                      double zdiff,
                                      GRID_FLOAT_TYPE vval00,
                                      GRID_FLOAT_TYPE vval01,
                                      GRID_FLOAT_TYPE vval10,
                                      GRID_FLOAT_TYPE vval11);
  template <typename GRID_FLOAT_TYPE>
  GRID_FLOAT_TYPE interpolateValues3D(double xdiff,
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

  void convertAngles(const TakeOffAngles &in, double &azim, double &dip);

  Grid _grid;
  const unsigned _quality_cutoff;
};

class VelGrid
{
public:
  VelGrid(const std::string &filePath, bool swapBytes);

  void open(Grid::OpenMode mode);
  bool isOpen();
  void close();

  double getVel(double lat, double lon, double depth);

  double getVelAtIndex(unsigned long long ix,
                       unsigned long long iy,
                       unsigned long long iz);

  bool is3D() const { return _grid.info.numx > 2; }

  bool isGlobal() const { return _grid.info.transform->getType() == "GLOBAL"; }

  const Grid::Info &getInfo() const { return _grid.info; }

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

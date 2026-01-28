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

#include "grid.h"
#include "../log.h"
#include "../utils.h"
#include "map_project.h"

#include <array>
#include <cstring>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;
using namespace HDD::Logger;
using TakeOffAngles = HDD::NLL::AngleGrid::TakeOffAngles;

namespace {

using namespace HDD;
using namespace HDD::NLL;

class GlobalTransform : public Transform
{
public:
  GlobalTransform(const std::string &type);
  virtual ~GlobalTransform() = default;
  void
  fromLatLon(double lat, double lon, double &xLoc, double &yLoc) const override;
  void
  toLatLon(double xLoc, double yLoc, double &lat, double &lon) const override;
  double distance(double xLoc1,
                  double yLoc1,
                  double xLoc2,
                  double yLoc2) const override;
  double distance(double xLoc1,
                  double yLoc1,
                  double zLoc1,
                  double xLoc2,
                  double yLoc2,
                  double zLoc2) const override;
};

class SimpleTransform : public Transform
{
public:
  /*
   * Adopting NLL constants to improve compatibility
   */
  static constexpr double AVG_ERAD = 6371.0;
  static constexpr double c111 =
      M_PI * AVG_ERAD / 180.; // kilometers per degree

public:
  SimpleTransform(const std::string &type,
                  double orgLat,
                  double orgLong,
                  double rot);
  void
  fromLatLon(double lat, double lon, double &xLoc, double &yLoc) const override;
  void
  toLatLon(double xLoc, double yLoc, double &lat, double &lon) const override;
};

class SDCTransform : public Transform
{
public:
  /*
   * Adopting NLL constants to improve compatibility
   */
  static constexpr double FLATTENING =
      1.0 / 298.26; // Earth flattening (WGS-72)
  static constexpr double ERAD               = 6378.135; // WGS-72
  static constexpr double MAP_TRANS_SDC_DRLT = 0.99330647;

public:
  SDCTransform(const std::string &type,
               double orgLat,
               double orgLong,
               double rot,
               double xltkm,
               double xlnkm);
  void
  fromLatLon(double lat, double lon, double &xLoc, double &yLoc) const override;
  void
  toLatLon(double xLoc, double yLoc, double &lat, double &lon) const override;

private:
  const double _xltkm;
  const double _xlnkm;
};

class LambertTransform : public Transform
{
public:
  LambertTransform(const std::string &type,
                   double orgLat,
                   double orgLong,
                   double rot,
                   const std::string &refEllip,
                   double pha,
                   double phb);
  void
  fromLatLon(double lat, double lon, double &xLoc, double &yLoc) const override;
  void
  toLatLon(double xLoc, double yLoc, double &lat, double &lon) const override;

private:
  std::string _refEllip;
  GMT::LAMBERT _proj;
};

class TransMercTransform : public Transform
{
public:
  TransMercTransform(const std::string &type,
                     double orgLat,
                     double orgLong,
                     double rot,
                     const std::string &refEllip,
                     bool useFalseEasting,
                     long falseEasting,
                     double mapScaleFactor);
  void
  fromLatLon(double lat, double lon, double &xLoc, double &yLoc) const override;
  void
  toLatLon(double xLoc, double yLoc, double &lat, double &lon) const override;

private:
  const std::string _refEllip;
  GMT::TRANS_MERCATOR _proj;
};

class AziEqdistTransform : public Transform
{
public:
  AziEqdistTransform(const std::string &type,
                     double orgLat,
                     double orgLong,
                     double rot,
                     const std::string &refEllip);
  void
  fromLatLon(double lat, double lon, double &xLoc, double &yLoc) const override;
  void
  toLatLon(double xLoc, double yLoc, double &lat, double &lon) const override;

private:
  const std::string _refEllip;
  GMT::AZIMUTHAL_EQUIDIST _proj;
};

GlobalTransform::GlobalTransform(const std::string &type)
    : Transform(type, 0, 0, 0)
{}

void GlobalTransform::fromLatLon(double lat,
                                 double lon,
                                 double &xLoc,
                                 double &yLoc) const
{
  xLoc = lon;
  yLoc = lat;
}

void GlobalTransform::toLatLon(double xLoc,
                               double yLoc,
                               double &lat,
                               double &lon) const
{
  lat = yLoc;
  lon = xLoc;
}

double GlobalTransform::distance(double xLoc1,
                                 double yLoc1,
                                 double xLoc2,
                                 double yLoc2) const
{
  return radToDeg(
      computeDistance(yLoc1, xLoc1, yLoc2, xLoc2, nullptr, nullptr, 0, true));
}

double GlobalTransform::distance(double xLoc1,
                                 double yLoc1,
                                 double zLoc1,
                                 double xLoc2,
                                 double yLoc2,
                                 double zLoc2) const
{
  throw Exception("GlobalTransform doesn't support 3D grids");
}

SimpleTransform::SimpleTransform(const std::string &type,
                                 double orgLat,
                                 double orgLong,
                                 double rot)
    : Transform(type, orgLat, orgLong, rot)
{}

void SimpleTransform::fromLatLon(double lat,
                                 double lon,
                                 double &xLoc,
                                 double &yLoc) const
{
  xLoc = normalizeLon(lon - _orgLong);
  xLoc *= c111 * std::cos(degToRad(lat));
  yLoc = (lat - _orgLat) * c111;
  rotate(xLoc, yLoc);
}

void SimpleTransform::toLatLon(double xLoc,
                               double yLoc,
                               double &lat,
                               double &lon) const
{
  inverseRotate(xLoc, yLoc);
  lat = _orgLat + yLoc / c111;
  lon = _orgLong + xLoc / (c111 * std::cos(degToRad(lat)));
  lon = normalizeLon(lon);
}

SDCTransform::SDCTransform(const std::string &type,
                           double orgLat,
                           double orgLong,
                           double rot,
                           double xltkm,
                           double xlnkm)
    : Transform(type, orgLat, orgLong, rot), _xltkm(xltkm), _xlnkm(xlnkm)
{}

void SDCTransform::fromLatLon(double lat,
                              double lon,
                              double &xLoc,
                              double &yLoc) const
{
  xLoc = normalizeLon(lon - _orgLong);
  yLoc = lat - _orgLat;
  double xlt1 =
      std::atan(MAP_TRANS_SDC_DRLT * std::tan(degToRad(lat + _orgLat) / 2.0));
  xLoc *= _xlnkm * std::cos(xlt1);
  yLoc *= _xltkm;
  rotate(xLoc, yLoc);
}

void SDCTransform::toLatLon(double xLoc,
                            double yLoc,
                            double &lat,
                            double &lon) const
{
  inverseRotate(xLoc, yLoc);
  yLoc /= _xltkm;
  lat = _orgLat + yLoc;
  double xlt1 =
      std::atan(MAP_TRANS_SDC_DRLT * std::tan(degToRad(lat + _orgLat) / 2.0));
  xLoc /= (_xlnkm * std::cos(xlt1));
  lon = _orgLong + xLoc;
  lon = normalizeLon(lon);
}

LambertTransform::LambertTransform(const std::string &type,
                                   double orgLat,
                                   double orgLong,
                                   double rot,
                                   const std::string &refEllip,
                                   double pha,
                                   double phb)
    : Transform(type, orgLat, orgLong, rot), _refEllip(refEllip)
{
  if (pha > 90 || pha < -90)
  {
    throw Exception("FirstStdParal must be in range -90,90");
  }
  if (phb > 90 || phb < -90)
  {
    throw Exception("SecondStdParal must be in range -90,90");
  }
  _proj = GMT::vlamb(_refEllip.c_str(), _orgLong, _orgLat, pha, phb);
}

void LambertTransform::fromLatLon(double lat,
                                  double lon,
                                  double &xLoc,
                                  double &yLoc) const
{
  GMT::lamb(_proj, lon, lat, &xLoc, &yLoc);
  xLoc /= 1000.0; /* m -> km */
  yLoc /= 1000.0; /* m -> km */
  rotate(xLoc, yLoc);
}

void LambertTransform::toLatLon(double xLoc,
                                double yLoc,
                                double &lat,
                                double &lon) const
{
  inverseRotate(xLoc, yLoc);
  GMT::ilamb(_proj, &lon, &lat, xLoc * 1000.0, yLoc * 1000.0);
  lon = normalizeLon(lon);
}

TransMercTransform::TransMercTransform(const std::string &type,
                                       double orgLat,
                                       double orgLong,
                                       double rot,
                                       const std::string &refEllip,
                                       bool useFalseEasting,
                                       long falseEasting,
                                       double mapScaleFactor)
    : Transform(type, orgLat, orgLong, rot), _refEllip(refEllip)
{
  _proj = GMT::vtm(_refEllip.c_str(), _orgLong, _orgLat, useFalseEasting,
                   falseEasting, mapScaleFactor);
}

void TransMercTransform::fromLatLon(double lat,
                                    double lon,
                                    double &xLoc,
                                    double &yLoc) const
{
  GMT::tm(_proj, lon, lat, &xLoc, &yLoc);
  xLoc /= 1000.0; /* m -> km */
  yLoc /= 1000.0; /* m -> km */
  rotate(xLoc, yLoc);
}

void TransMercTransform::toLatLon(double xLoc,
                                  double yLoc,
                                  double &lat,
                                  double &lon) const
{
  inverseRotate(xLoc, yLoc);
  GMT::itm(_proj, &lon, &lat, xLoc * 1000.0, yLoc * 1000.0);
  lon = normalizeLon(lon);
}

AziEqdistTransform::AziEqdistTransform(const std::string &type,
                                       double orgLat,
                                       double orgLong,
                                       double rot,
                                       const std::string &refEllip)
    : Transform(type, orgLat, orgLong, rot), _refEllip(refEllip)
{
  _proj = GMT::vazeqdist(_refEllip.c_str(), _orgLong, _orgLat);
}

void AziEqdistTransform::fromLatLon(double lat,
                                    double lon,
                                    double &xLoc,
                                    double &yLoc) const
{
  GMT::azeqdist(_proj, lon, lat, &xLoc, &yLoc);
  xLoc /= 1000.0; /* m -> km */
  yLoc /= 1000.0; /* m -> km */
  rotate(xLoc, yLoc);
}

void AziEqdistTransform::toLatLon(double xLoc,
                                  double yLoc,
                                  double &lat,
                                  double &lon) const
{
  inverseRotate(xLoc, yLoc);
  GMT::iazeqdist(_proj, &lon, &lat, xLoc * 1000.0, yLoc * 1000.0);
  lon = normalizeLon(lon);
}

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
  // Lambda for linear interpolation
  auto lerp = [](double a, double b, double t) { return a + t * (b - a); };

  // Step 1: Interpolate along Z-axis (reducing 4 points to 2)
  double z0 = lerp(vval00, vval01, zdiff);
  double z1 = lerp(vval10, vval11, zdiff);

  // Step 2: Interpolate along X-axis (final value)
  return lerp(z0, z1, xdiff);
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
  // Define the linear interpolation lambda
  auto lerp = [](double a, double b, double t) { return a + t * (b - a); };

  // Step 1: Interpolate along Z-axis (reducing 8 points to 4)
  double z00 = lerp(vval000, vval001, zdiff);
  double z01 = lerp(vval010, vval011, zdiff);
  double z10 = lerp(vval100, vval101, zdiff);
  double z11 = lerp(vval110, vval111, zdiff);

  // Step 2: Interpolate along Y-axis (reducing 4 points to 2)
  double y0 = lerp(z00, z01, ydiff);
  double y1 = lerp(z10, z11, ydiff);

  // Step 3: Interpolate along X-axis (final value)
  return lerp(y0, y1, xdiff);
}

/*
 * Similar to interpolateSquareLagrange but for angles. Since they can wrap
 * at 360 degree, they need a special function
 */
double interpolateSquareAngles(
    double xdiff, double zdiff, double v00, double v01, double v10, double v11)
{
  // 1. Interpolate X-components
  double x_comp = interpolateSquareLagrange(
      xdiff, zdiff, cos(degToRad(v00)), cos(degToRad(v01)), cos(degToRad(v10)),
      cos(degToRad(v11)));

  // 2. Interpolate Z-components
  double z_comp = interpolateSquareLagrange(
      xdiff, zdiff, sin(degToRad(v00)), sin(degToRad(v01)), sin(degToRad(v10)),
      sin(degToRad(v11)));

  // 3. Convert back to degrees
  double result = atan2(z_comp, x_comp);

  // Normalize to [0, 360)
  result = radToDeg(result);
  if (result < 0) result += 360.0;
  return result;
}

/*
 * Similar to interpolateCubeLagrange but for angles. Since they can wrap
 * at 360 degree, they need a special function
 */
double interpolateCubeAngles(double xdiff,
                             double ydiff,
                             double zdiff,
                             double v000,
                             double v001,
                             double v010,
                             double v011,
                             double v100,
                             double v101,
                             double v110,
                             double v111)
{
  // 1. Interpolate the X-components (cosines)
  double cos_comp = interpolateCubeLagrange(
      xdiff, ydiff, zdiff, cos(degToRad(v000)), cos(degToRad(v001)),
      cos(degToRad(v010)), cos(degToRad(v011)), cos(degToRad(v100)),
      cos(degToRad(v101)), cos(degToRad(v110)), cos(degToRad(v111)));

  // 2. Interpolate the Y-components (sines)
  double sin_comp = interpolateCubeLagrange(
      xdiff, ydiff, zdiff, sin(degToRad(v000)), sin(degToRad(v001)),
      sin(degToRad(v010)), sin(degToRad(v011)), sin(degToRad(v100)),
      sin(degToRad(v101)), sin(degToRad(v110)), sin(degToRad(v111)));

  // 3. Reconstruct the angle from the average vector components
  double result = atan2(sin_comp, cos_comp);

  // Normalize result to [0, 360)
  result = radToDeg(result);
  if (result < 0) result += 360.0;
  return result;
}

} // namespace

namespace HDD {
namespace NLL {

void Transform::rotate(double &xLoc, double &yLoc) const
{
  const double xtemp = xLoc;
  const double ytemp = yLoc;
  xLoc               = xtemp * _cosang - ytemp * _sinang;
  yLoc               = ytemp * _cosang + xtemp * _sinang;
}

void Transform::inverseRotate(double &xLoc, double &yLoc) const
{
  const double xtemp = xLoc;
  const double ytemp = yLoc;
  xLoc               = xtemp * _cosang + ytemp * _sinang;
  yLoc               = ytemp * _cosang - xtemp * _sinang;
}

unique_ptr<Transform> Transform::parse(const std::vector<string> &tokens)
{
  if (tokens.at(0) != "TRANSFORM" && tokens.at(0) != "TRANS")
  {
    throw Exception("Malformed transform line");
  }

  string type = tokens.at(1);

  if (type == "GLOBAL")
  {
    return unique_ptr<Transform>(new GlobalTransform(type));
  }
  else if (type == "SDC" || type == "SIMPLE")
  {
    if (tokens.at(2) != "LatOrig" || tokens.at(4) != "LongOrig" ||
        tokens.at(6) != "RotCW")
    {
      throw Exception("Cannot parse grid header");
    }

    double orgLat  = std::stod(tokens.at(3));
    double orgLong = std::stod(tokens.at(5));
    double rot     = std::stod(tokens.at(7));

    if (type == "SIMPLE")
    {
      return unique_ptr<Transform>(
          new SimpleTransform(type, orgLat, orgLong, rot));
    }
    else
    {
      //  conversion factor for latitude
      double dlt1 = std::atan(SDCTransform::MAP_TRANS_SDC_DRLT *
                              std::tan(degToRad(orgLat)));
      double dlt2 = std::atan(SDCTransform::MAP_TRANS_SDC_DRLT *
                              std::tan(degToRad(orgLat + 1.0)));
      double del  = dlt2 - dlt1;
      double r    = SDCTransform::ERAD *
                 (1.0 - square(std::sin(dlt1)) * SDCTransform::FLATTENING);
      double sdc_xltkm = r * del;
      //  conversion factor for longitude
      del              = std::acos(1.0 -
                                   (1.0 - std::cos(degToRad(1))) * square(std::cos(dlt1)));
      double bc        = r * del;
      double sdc_xlnkm = bc / std::cos(dlt1);

      return unique_ptr<Transform>(
          new SDCTransform(type, orgLat, orgLong, rot, sdc_xltkm, sdc_xlnkm));
    }
  }
  else if (type == "LAMBERT")
  {
    if (tokens.at(2) != "RefEllipsoid" || tokens.at(4) != "LatOrig" ||
        tokens.at(6) != "LongOrig" || tokens.at(8) != "FirstStdParal" ||
        tokens.at(10) != "SecondStdParal" || tokens.at(12) != "RotCW")
    {
      throw Exception("Cannot parse grid header");
    }

    string refEllip = tokens.at(3);
    double orgLat   = std::stod(tokens.at(5));
    double orgLong  = std::stod(tokens.at(7));
    double pha      = std::stod(tokens.at(9));
    double phb      = std::stod(tokens.at(11));
    double rot      = std::stod(tokens.at(13));

    return unique_ptr<Transform>(
        new LambertTransform(type, orgLat, orgLong, rot, refEllip, pha, phb));
  }
  else if (type == "TRANS_MERC")
  {
    if (tokens.at(2) != "RefEllipsoid" || tokens.at(4) != "LatOrig" ||
        tokens.at(6) != "LongOrig" || tokens.at(8) != "RotCW")
    {
      throw Exception("Cannot parse grid header");
    }

    string refEllip = tokens.at(3);
    double orgLat   = std::stod(tokens.at(5));
    double orgLong  = std::stod(tokens.at(7));
    double rot      = std::stod(tokens.at(9));

    bool useFalseEasting  = false;
    long falseEasting     = 500000;
    double mapScaleFactor = 1.0;

    if (tokens.size() > 10)
    {
      if (tokens.at(10) != "UseFalseEasting")
        throw Exception("Cannot parse grid header");
      useFalseEasting = std::stod(tokens.at(11)) != 0;
    }

    if (tokens.size() > 12)
    {
      if (tokens.at(12) != "FalseEasting")
        throw Exception("Cannot parse grid header");
      falseEasting = std::stol(tokens.at(13));
    }

    if (tokens.size() > 14)
    {
      if (tokens.at(14) != "ScaleFactor")
        throw Exception("Cannot parse grid header");
      mapScaleFactor = std::stod(tokens.at(15));
    }

    return unique_ptr<Transform>(
        new TransMercTransform(type, orgLat, orgLong, rot, refEllip,
                               useFalseEasting, falseEasting, mapScaleFactor));
  }
  else if (type == "AZIMUTHAL_EQUIDIST")
  {
    if (tokens.at(2) != "RefEllipsoid" || tokens.at(4) != "LatOrig" ||
        tokens.at(6) != "LongOrig" || tokens.at(8) != "RotCW")
    {
      throw Exception("Cannot parse grid header");
    }
    string refEllip = tokens.at(3);
    double orgLat   = std::stod(tokens.at(5));
    double orgLong  = std::stod(tokens.at(7));
    double rot      = std::stod(tokens.at(9));

    return unique_ptr<Transform>(
        new AziEqdistTransform(type, orgLat, orgLong, rot, refEllip));
  }
  else
  {
    string msg = strf("Unsupported transform %s", type.c_str());
    throw Exception(msg.c_str());
  }
}

Transform::Transform(const std::string &type,
                     double orgLat,
                     double orgLong,
                     double rot)
    : _type(type), _orgLat(orgLat), _orgLong(orgLong), _rot(rot),
      _cosang(std::cos(-degToRad(rot))), _sinang(std::sin(-degToRad(rot)))
{
  if (_orgLat > 90 || _orgLat < -90)
  {
    throw Exception("Origin latitude must be in range -90,90");
  }

  if (_orgLong > 180 || _orgLong < -180)
  {
    throw Exception("Origin longitude must be in range -180,180");
  }

  if (_rot > 360 || _rot < -360)
  {
    throw Exception("Rotation must be in range -360,360");
  }
}

double Transform::fromLatLonAngle(double latlonAngle) const
{
  double angle = latlonAngle + _rot;
  if (angle < 0.0)
    angle += 360.0;
  else if (angle > 360.0)
    angle -= 360.0;
  return angle;
}

double Transform::toLatLonAngle(double rectAngle) const
{
  double angle = rectAngle - _rot;
  if (angle < 0.0)
    angle += 360.0;
  else if (angle > 360.0)
    angle -= 360.0;
  return (angle);
}

double Transform::distance(double xLoc1,
                           double yLoc1,
                           double xLoc2,
                           double yLoc2) const
{
  return std::sqrt(square(xLoc2 - xLoc1) + square(yLoc2 - yLoc1));
}

double Transform::distance(double xLoc1,
                           double yLoc1,
                           double zLoc1,
                           double xLoc2,
                           double yLoc2,
                           double zLoc2) const
{
  return std::sqrt(square(xLoc2 - xLoc1) + square(yLoc2 - yLoc1) +
                   square(zLoc2 - zLoc1));
}

Grid::Grid(Type gridType, const std::string &filePath, bool swapBytes)
    : info(parse(filePath, gridType, swapBytes))
{

  if (!pathExists(info.bufFilePath))
  {
    string msg =
        strf("Cannot find grid data file %s", info.bufFilePath.c_str());
    throw Exception(msg.c_str());
  }
}

Grid::~Grid() { close(); }

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
    if (!tokens.empty() && tokens.at(0).empty()) tokens.erase(tokens.begin());

    // skip empty lines
    if (tokens.empty()) continue;

    if (parsedLines == 0 &&
        (tokens.size() == 10 || tokens.size() == 11)) // should be line 1
    {
      info.numx  = std::stoull(tokens.at(0));
      info.numy  = std::stoull(tokens.at(1));
      info.numz  = std::stoull(tokens.at(2));
      info.origx = std::stod(tokens.at(3));
      info.origy = std::stod(tokens.at(4));
      info.origz = std::stod(tokens.at(5));
      info.dx    = std::stod(tokens.at(6));
      info.dy    = std::stod(tokens.at(7));
      info.dz    = std::stod(tokens.at(8));
      info.type  = tokens.at(9);
      info.useDouble =
          (tokens.size() == 11) ? (tokens.at(10) == "DOUBLE") : false;
      ++parsedLines;
    }
    else if (parsedLines == 1 && tokens.size() == 4 &&
             (gridType == Type::time || gridType == Type::angle))
    {
      // this should be line 2 (if present)
      info.label = tokens.at(0);
      info.srcex = std::stod(tokens.at(1));
      info.srcey = std::stod(tokens.at(2));
      info.srcez = std::stod(tokens.at(3));
      ++parsedLines;
    }
    else if (tokens.at(0) == "TRANSFORM" || tokens.at(0) == "TRANS")
    {
      info.transform = Transform::parse(tokens);
      ++parsedLines;
    }
  }

  if (((gridType == Type::time || gridType == Type::angle) &&
       parsedLines != 3) ||
      (gridType == Type::velocity && parsedLines != 2))
  {
    string msg =
        strf("Cannot load grid header file %s", info.hdrFilePath.c_str());
    throw Exception(msg.c_str());
  }

  // make sure that dx for 2D grids is non-zero
  if (info.numx == 1) info.dx = 1.0;

  if (info.useDouble && info.swapBytes)
  {
    throw Exception(
        "Grid files with DOUBLE values and byte swapping are not supported");
  }

  return info;
}

void Grid::open(OpenMode mode)
{
  if (!isOpen())
  {
    _mode = mode;

    if (_fd < 0)
    {
      _fd = ::open(info.bufFilePath.c_str(), O_RDONLY);
    }

    if (_fd < 0)
    {
      throw Exception("Failed to open grid file " + info.bufFilePath);
    }

    struct stat st;
    if (::fstat(_fd, &st) != 0)
    {
      throw Exception(
          strf("Failed to stat grid file %s", info.bufFilePath.c_str()));
    }

    _mappedSize = static_cast<std::size_t>(st.st_size);

    if (_mode == OpenMode::Mmap)
    {
      if (!_mapped || _mapped == MAP_FAILED)
      {
        _mapped = ::mmap(nullptr, _mappedSize, PROT_READ, MAP_SHARED, _fd, 0);
        if (_mapped == MAP_FAILED)
        {
          throw Exception(
              strf("Failed to mmap grid file %s", info.bufFilePath.c_str()));
        }
      }

      ::close(_fd);
      _fd = -1;
    }
    else if (_mode == OpenMode::LoadIntoMemory)
    {
      _bytes.resize(_mappedSize);
      ssize_t bytesRead = ::pread(_fd, _bytes.data(), _mappedSize, 0);

      if (bytesRead != (ssize_t)_mappedSize)
      {
        throw Exception("pread failed or reached EOF");
      }

      ::close(_fd);
      _fd = -1;
    }
  }
}

bool Grid::isOpen()
{
  if (_mode == OpenMode::Mmap)
  {
    return _mapped && _mapped != MAP_FAILED;
  }
  else if (_mode == OpenMode::Pread)
  {
    return _fd >= 0;
  }
  else
  {
    return _bytes.size() > 0;
  }
}

void Grid::close()
{
  if (isOpen())
  {
    logDebugF("Closing grid file %s", info.bufFilePath.c_str());
  }

  if (_mapped && _mapped != MAP_FAILED)
  {
    ::munmap(_mapped, _mappedSize);
    _mapped = nullptr;
  }

  if (_fd >= 0)
  {
    ::close(_fd);
    _fd = -1;
  }

  _bytes.clear();
  _bytes.shrink_to_fit();
  _mappedSize = 0;
}

template <typename GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE Grid::getValueAtIndex(unsigned long long ix,
                                      unsigned long long iy,
                                      unsigned long long iz)
{
  if (!isOpen())
  {
    throw Exception("Grid file not open: " + info.bufFilePath);
  }

  if (!info.isIndexInside(ix, iy, iz))
  {
    throw Exception("Requested index is out of grid boundaries");
  }

  unsigned long long index = ix * info.numy * info.numz + iy * info.numz + iz;
  unsigned long long byteOffset = index * sizeof(GRID_FLOAT_TYPE);

  GRID_FLOAT_TYPE value;
  if (_mode == OpenMode::Mmap)
  {
    if (byteOffset + sizeof(GRID_FLOAT_TYPE) > _mappedSize)
    {
      throw Exception("Requested data past end of grid file " +
                      info.bufFilePath);
    }

    std::memcpy(&value, static_cast<const char *>(_mapped) + byteOffset,
                sizeof(GRID_FLOAT_TYPE));
  }
  else if (_mode == OpenMode::Pread)
  {
    ssize_t bytesRead = ::pread(_fd, &value, sizeof(value), byteOffset);

    if (bytesRead != sizeof(value))
    {
      throw Exception("pread failed or reached EOF");
    }
  }
  else
  {
    if (byteOffset + sizeof(GRID_FLOAT_TYPE) > _bytes.size())
    {
      throw Exception("Requested data past end of grid file " +
                      info.bufFilePath);
    }
    std::memcpy(&value, _bytes.data() + byteOffset, sizeof(GRID_FLOAT_TYPE));
  }

  if (info.swapBytes)
  {
    if constexpr (sizeof(GRID_FLOAT_TYPE) == 4)
    {
      auto *b = reinterpret_cast<unsigned char *>(&value);
      std::swap(b[0], b[3]);
      std::swap(b[1], b[2]);
    }
    else if constexpr (sizeof(GRID_FLOAT_TYPE) == 8)
    {
      auto *b = reinterpret_cast<unsigned char *>(&value);
      std::reverse(b, b + 8);
    }
  }
  return value;
}

template <typename GRID_FLOAT_TYPE>
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
  if (!info.isLocationInside(xloc, yloc, zloc))
  {
    string msg = strf("Requested location is out of grid boundaries "
                      "(xloc %.2f yloc %.2f zloc %.2f - grid %s "
                      "origx %.3f origy %.3f origz %.3f "
                      "dx %.2f dy %.2f dz %.2f "
                      "numx %llu numy %llu numz %llu)",
                      xloc, yloc, zloc, info.hdrFilePath.c_str(), info.origx,
                      info.origy, info.origz, info.dx, info.dy, info.dz,
                      info.numx, info.numy, info.numz);
    throw Exception(msg.c_str());
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

template <typename GRID_FLOAT_TYPE>
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

  if (!info.isLocationInside(xloc, yloc, zloc))
  {
    string msg = strf("Requested location is out of grid boundaries "
                      "(xloc %.2f yloc %.2f zloc %.2f - grid %s "
                      "origx %.3f origy %.3f origz %.3f "
                      "dx %.2f dy %.2f dz %.2f "
                      "numx %llu numy %llu numz %llu)",
                      xloc, yloc, zloc, info.hdrFilePath.c_str(), info.origx,
                      info.origy, info.origz, info.dx, info.dy, info.dz,
                      info.numx, info.numy, info.numz);
    throw Exception(msg.c_str());
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

template <typename GRID_FLOAT_TYPE>
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

template <typename GRID_FLOAT_TYPE>
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

TimeGrid::TimeGrid(const std::string &filePath, bool swapBytes)
    : _grid(Grid::Type::time, filePath, swapBytes)
{
  if (_grid.info.type != "TIME" && _grid.info.type != "TIME2D")
  {
    string msg = strf("Unrecognized time grid type %s (%s)",
                      _grid.info.type.c_str(), _grid.info.hdrFilePath.c_str());
    throw Exception(msg.c_str());
  }
}

void TimeGrid::open(Grid::OpenMode mode) { _grid.open(mode); }
bool TimeGrid::isOpen() { return _grid.isOpen(); }
void TimeGrid::close() { _grid.close(); }

double TimeGrid::getTimeAtIndex(unsigned long long ix,
                                unsigned long long iy,
                                unsigned long long iz)
{
  if (_grid.info.useDouble)
  {
    return _grid.getValueAtIndex<double>(ix, iy, iz);
  }
  else
  {
    return _grid.getValueAtIndex<float>(ix, iy, iz);
  }
}

double TimeGrid::getTime(double lat, double lon, double depth)
{
  if (_grid.info.useDouble)
  {
    if (is3D())
    {
      return _grid.getValue3D<double>(lat, lon, depth,
                                      interpolateValues3D<double>);
    }
    else
    {
      return _grid.getValue2D<double>(lat, lon, depth,
                                      interpolateValues2D<double>);
    }
  }
  else
  {
    if (is3D())
    {
      return _grid.getValue3D<float>(lat, lon, depth,
                                     interpolateValues3D<float>);
    }
    else
    {
      return _grid.getValue2D<float>(lat, lon, depth,
                                     interpolateValues2D<float>);
    }
  }
}

template <typename GRID_FLOAT_TYPE>
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
    throw Exception("Negative times found in the grid file");
  }
  return interpolateCubeLagrange(xdiff, ydiff, zdiff, vval000, vval001, vval010,
                                 vval011, vval100, vval101, vval110, vval111);
}

template <typename GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE TimeGrid::interpolateValues2D(double xdiff,
                                              double zdiff,
                                              GRID_FLOAT_TYPE vval00,
                                              GRID_FLOAT_TYPE vval01,
                                              GRID_FLOAT_TYPE vval10,
                                              GRID_FLOAT_TYPE vval11)
{
  if (vval00 < 0.0 || vval01 < 0.0 || vval10 < 0.0 || vval11 < 0.0)
  {
    throw Exception("Negative times found in the grid file");
  }
  return interpolateSquareLagrange(xdiff, zdiff, vval00, vval01, vval10,
                                   vval11);
}

AngleGrid::AngleGrid(const std::string &filePath,
                     bool swapBytes,
                     unsigned quality_cutoff)
    : _grid(Grid::Type::angle, filePath, swapBytes),
      _quality_cutoff(quality_cutoff)
{
  if (_grid.info.type != "ANGLE" && _grid.info.type != "ANGLE2D")
  {
    string msg =
        strf("Unrecognized angle grid type %s", _grid.info.type.c_str());
    throw Exception(msg.c_str());
  }

  //
  // TakeOffAngles uses bit-fields, but their memory arrangment
  // (how bits are packed together) is not standard. So let's make
  // sure the compiler did a good job
  //
  if (sizeof(TakeOffAngles) != sizeof(float))
  {
    throw Exception("Internal error: sizeof(TakeOffAngles) != sizeof(float)");
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
    throw Exception("Internal error: TakeOffAngles memory mapping is not ok");
  }

  // assumes angle files store floats only (never double)
  if (_grid.info.useDouble)
  {
    throw Exception("Angle grid files with DOUBLE values are not "
                    "supported, only FLOAT allowed");
  }
}

void AngleGrid::open(Grid::OpenMode mode) { _grid.open(mode); }
bool AngleGrid::isOpen() { return _grid.isOpen(); }
void AngleGrid::close() { _grid.close(); }

void AngleGrid::getAnglesAtIndex(unsigned long long ix,
                                 unsigned long long iy,
                                 unsigned long long iz,
                                 double &azim,
                                 double &dip)
{
  TakeOffAngles angles = _grid.getValueAtIndex<TakeOffAngles>(ix, iy, iz);
  convertAngles(angles, azim, dip);
}

void AngleGrid::getAngles(
    double lat, double lon, double depth, double &azim, double &dip)
{
  TakeOffAngles angles;

  if (is3D())
  {
    angles = _grid.getValue3D<TakeOffAngles>(
        lat, lon, depth,
        std::bind(
            &AngleGrid::interpolateValues3D<TakeOffAngles>, this,
            std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
            std::placeholders::_4, std::placeholders::_5, std::placeholders::_6,
            std::placeholders::_7, std::placeholders::_8, std::placeholders::_9,
            std::placeholders::_10, std::placeholders::_11));
  }
  else
  {
    angles = _grid.getValue2D<TakeOffAngles>(
        lat, lon, depth,
        std::bind(&AngleGrid::interpolateValues2D<TakeOffAngles>, this,
                  std::placeholders::_1, std::placeholders::_2,
                  std::placeholders::_3, std::placeholders::_4,
                  std::placeholders::_5, std::placeholders::_6));
  }
  convertAngles(angles, azim, dip);
}

void AngleGrid::convertAngles(const TakeOffAngles &angles,
                              double &azim,
                              double &dip)
{
  if (angles.quality < _quality_cutoff)
  {
    string msg = strf("Angle quality too low (%u) for selected cut off (%u)",
                      angles.quality, _quality_cutoff);
    throw Exception(msg.c_str());
  }

  if (is3D())
  {
    azim = angles.azimuth / 10.0; // tenths of degree -> degree
    azim = _grid.info.transform->toLatLonAngle(azim);
  }
  else
  {
    azim = std::numeric_limits<double>::quiet_NaN();
  }
  dip = angles.dip / 10.0; // tenths of degree -> degree
}

template <typename GRID_FLOAT_TYPE>
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
  // get lowest quality angles
  unsigned short lowest_qual = std::min(std::initializer_list<unsigned short>(
      {vval000.quality, vval001.quality, vval010.quality, vval011.quality,
       vval100.quality, vval101.quality, vval110.quality, vval111.quality}));

  // due to the azimuth angle range that span 0/360, we must use a special
  // interpolation function that handle the wrapping at 360
  unsigned short azim_interp =
      interpolateCubeAngles(
          xdiff, ydiff, zdiff, vval000.azimuth / 10., vval001.azimuth / 10.,
          vval010.azimuth / 10., vval011.azimuth / 10., vval100.azimuth / 10.,
          vval101.azimuth / 10., vval110.azimuth / 10., vval111.azimuth / 10.) *
      10.0;

  // dip angles are in the range 0-180, so simple interpolation
  unsigned short dip_interp = interpolateCubeLagrange(
      xdiff, ydiff, zdiff, vval000.dip, vval001.dip, vval010.dip, vval011.dip,
      vval100.dip, vval101.dip, vval110.dip, vval111.dip);

  return {lowest_qual, dip_interp, azim_interp};
}

template <typename GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE AngleGrid::interpolateValues2D(double xdiff,
                                               double zdiff,
                                               GRID_FLOAT_TYPE vval00,
                                               GRID_FLOAT_TYPE vval01,
                                               GRID_FLOAT_TYPE vval10,
                                               GRID_FLOAT_TYPE vval11)
{
  // get lowest quality angles
  unsigned short lowest_qual = std::min(std::initializer_list<unsigned short>(
      {vval00.quality, vval01.quality, vval10.quality, vval11.quality}));

  // There is no azimuth in 2D Angle Grids: skip the value computation
  unsigned short azim_interp = 0;
  //    interpolateSquareAngles(xdiff, zdiff, vval00.azimuth / 10.,
  //                            vval01.azimuth / 10., vval10.azimuth / 10.,
  //                            vval11.azimuth / 10.) * 10.0;

  // dip angles are in the range 0-180, so simple interpolation
  unsigned short dip_interp = interpolateSquareLagrange(
      xdiff, zdiff, vval00.dip, vval01.dip, vval10.dip, vval11.dip);

  return {lowest_qual, dip_interp, azim_interp};
}

VelGrid::VelGrid(const std::string &filePath, bool swapBytes)
    : _grid(Grid::Type::velocity, filePath, swapBytes)
{
  if (_grid.info.numx < 2)
  {
    string msg =
        strf("Velocity grid must have xNum greater than 2, found %llu (%s)",
             _grid.info.numx, _grid.info.hdrFilePath.c_str());
    throw Exception(msg.c_str());
  }

  if (_grid.info.type == "VELOCITY_METERS")
    convertUnits = [](double vel) -> double { return vel / 1000.0; };
  else if (_grid.info.type == "SLOWNESS")
    convertUnits = [](double vel) -> double { return 1.0 / vel; };
  else if (_grid.info.type == "SLOW_LEN")
  {
    convertUnits = [this](double vel) -> double {
      return 1.0 / (vel / _grid.info.dy);
    };
  }
  else if (_grid.info.type == "VEL2")
    convertUnits = [](double vel) -> double { return std::sqrt(vel); };
  else if (_grid.info.type == "SLOW2")
    convertUnits = [](double vel) -> double { return std::sqrt(1.0 / vel); };
  else if (_grid.info.type == "SLOW2_METERS")
    convertUnits = [](double vel) -> double {
      return std::sqrt(1.0 / vel) / 1000.0;
    };
  else if (_grid.info.type == "VELOCITY")
    convertUnits = [](double vel) -> double { return vel; };
  else
  {
    string msg =
        strf("Unrecognized velocity grid type %s", _grid.info.type.c_str());
    throw Exception(msg.c_str());
  }
}

void VelGrid::open(Grid::OpenMode mode) { _grid.open(mode); }
bool VelGrid::isOpen() { return _grid.isOpen(); }
void VelGrid::close() { _grid.close(); }

double VelGrid::getVelAtIndex(unsigned long long ix,
                              unsigned long long iy,
                              unsigned long long iz)
{
  double velAtSrc;
  if (_grid.info.useDouble)
  {
    velAtSrc = _grid.getValueAtIndex<double>(ix, iy, iz);
  }
  else
  {
    velAtSrc = _grid.getValueAtIndex<float>(ix, iy, iz);
  }
  // velocity -> [km/sec]
  return convertUnits(velAtSrc);
}

double VelGrid::getVel(double lat, double lon, double depth)
{
  double velAtSrc;

  if (_grid.info.useDouble)
  {
    if (is3D())
    {
      velAtSrc = _grid.getValue3D<double>(lat, lon, depth,
                                          interpolateValues3D<double>);
    }
    else
    {
      velAtSrc = _grid.getValue2D<double>(lat, lon, depth,
                                          interpolateValues2D<double>);
    }
  }
  else
  {
    if (is3D())
    {
      velAtSrc =
          _grid.getValue3D<float>(lat, lon, depth, interpolateValues3D<float>);
    }
    else
    {
      velAtSrc =
          _grid.getValue2D<float>(lat, lon, depth, interpolateValues2D<float>);
    }
  }

  // velocity -> [km/sec]
  return convertUnits(velAtSrc);
}

template <typename GRID_FLOAT_TYPE>
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
    throw Exception("Negative velocities found in the grid file");
  }
  return interpolateCubeLagrange(xdiff, ydiff, zdiff, vval000, vval001, vval010,
                                 vval011, vval100, vval101, vval110, vval111);
}

template <typename GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE VelGrid::interpolateValues2D(double xdiff,
                                             double zdiff,
                                             GRID_FLOAT_TYPE vval00,
                                             GRID_FLOAT_TYPE vval01,
                                             GRID_FLOAT_TYPE vval10,
                                             GRID_FLOAT_TYPE vval11)
{
  if (vval00 < 0.0 || vval01 < 0.0 || vval10 < 0.0 || vval11 < 0.0)
  {
    throw Exception("Negative velocities found in the grid file");
  }
  return interpolateSquareLagrange(xdiff, zdiff, vval00, vval01, vval10,
                                   vval11);
}

} // namespace NLL
} // namespace HDD

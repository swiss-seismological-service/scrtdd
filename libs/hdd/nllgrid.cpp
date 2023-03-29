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
#include "log.h"
#include "map_project.h"
#include "utils.h"

#include <array>
#include <cstring>

using namespace std;
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
  virtual ~SimpleTransform() = default;
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
  virtual ~SDCTransform() = default;
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
  virtual ~LambertTransform() = default;
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
  virtual ~TransMercTransform() = default;
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
  virtual ~AziEqdistTransform() = default;
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
  return computeDistance(yLoc1, xLoc1, yLoc2, xLoc2);
}

double GlobalTransform::distance(double xLoc1,
                                 double yLoc1,
                                 double zLoc1,
                                 double xLoc2,
                                 double yLoc2,
                                 double zLoc2) const
{
  return computeDistance(yLoc1, xLoc1, zLoc1, yLoc2, xLoc2, zLoc2);
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
  rotate(xLoc, yLoc);
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
  rotate(xLoc, yLoc);
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
  rotate(xLoc, yLoc);
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
  rotate(xLoc, yLoc);
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
  rotate(xLoc, yLoc);
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

namespace HDD {
namespace NLL {

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
      _angle(-degToRad(_rot)), _cosang(std::cos(_angle)),
      _sinang(std::sin(_angle))
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

  if (!pathExists(info.bufFilePath))
  {
    string msg =
        strf("Cannot find grid data file %s", info.bufFilePath.c_str());
    throw Exception(msg.c_str());
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
template <typename GRID_FLOAT_TYPE>
GRID_FLOAT_TYPE Grid::getValueAtIndex(unsigned long long ix,
                                      unsigned long long iy,
                                      unsigned long long iz)
{
  if (!isIndexInside(ix, iy, iz))
  {
    throw Exception("Requested index is out of grid boundaries");
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
    string msg = strf("Error while reading grid file %s (%s)",
                      info.bufFilePath.c_str(), e.what());
    throw Exception(msg.c_str());
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
  if (!isLocationInside(xloc, yloc, zloc))
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

  if (!isLocationInside(xloc, yloc, zloc))
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
    : _grid(Grid::Type::time, basePath, station, phaseType, swapBytes)
{
  if (_grid.info.type != "TIME" && _grid.info.type != "TIME2D")
  {
    string msg = strf("Unrecognized time grid type %s (%s)",
                      _grid.info.type.c_str(), _grid.info.hdrFilePath.c_str());
    throw Exception(msg.c_str());
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

AngleGrid::AngleGrid(const std::string &basePath,
                     const Catalog::Station &station,
                     const std::string &phaseType,
                     bool swapBytes)
    : _grid(Grid::Type::angle, basePath, station, phaseType, swapBytes)
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

void AngleGrid::getAngles(
    double lat, double lon, double depth, double &azim, double &dip)
{
  TakeOffAngles angles;

  if (is3D())
  {
    angles = _grid.getValue3D<TakeOffAngles>(
        lat, lon, depth, interpolateValues3D<TakeOffAngles>);
  }
  else
  {
    angles = _grid.getValue2D<TakeOffAngles>(
        lat, lon, depth, interpolateValues2D<TakeOffAngles>);
  }

  if (angles.quality < QUALITY_CUTOFF)
  {
    // this is not an error, the code handles nans
    azim = std::numeric_limits<double>::quiet_NaN();
    dip  = std::numeric_limits<double>::quiet_NaN();
    return;
  }

  if (is3D())
  {
    azim = angles.azimuth / 10.0; // tenths of degree -> degree
    azim = _grid.info.transform->toLatLonAngle(azim);
    azim = degToRad(azim);
  }
  else
  {
    azim = std::numeric_limits<double>::quiet_NaN();
  }
  dip = (angles.dip / 10.0);
  dip = degToRad(dip);
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
  return interpolateCubeAngles(xdiff, ydiff, zdiff, vval000, vval001, vval010,
                               vval011, vval100, vval101, vval110, vval111,
                               QUALITY_CUTOFF);
}

template <typename GRID_FLOAT_TYPE>
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
    : _grid(Grid::Type::velocity, basePath, station, phaseType, swapBytes)
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

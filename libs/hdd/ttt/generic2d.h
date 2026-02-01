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

#ifndef __HDD_TTT_GENERIC2D_H__
#define __HDD_TTT_GENERIC2D_H__

#include "../ttt.h"
#include "../utils.h"
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace HDD {
namespace TTT {

class Generic2D : public HDD::TravelTimeTable
{
public:
  class OutOfRange : public HDD::Exception
  {
  public:
    using Exception::Exception;
  };

public:
  Generic2D(const std::string &gridPath, const std::string &gridModel);

  Generic2D(const Generic2D &other)            = default;
  Generic2D &operator=(const Generic2D &other) = default;

  Generic2D(Generic2D &&other)            = default;
  Generic2D &operator=(Generic2D &&other) = default;

  Generic2D()  = default;
  ~Generic2D() = default;

  double compute(double eventLat,
                 double eventLon,
                 double eventDepth,
                 double stationLat,
                 double stationLon,
                 double stationElevation,
                 const std::string &phaseType) override;

  void compute(double eventLat,
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
               double &dtdh) override;

  std::unordered_set<std::string> availablePhases() const;

public:
  static const std::string COL_ELEVATION;
  static const std::string COL_DISTDEG;
  static const std::string COL_DISTKM;
  static const std::string COL_DEPTH;
  static const std::string COL_TTIME;
  static const std::string COL_TAKEOFF;
  static const std::string COL_DTDD;
  static const std::string COL_DTDH;
  static const std::string COL_VELOCITY;

private:
  struct GridData
  {
    std::vector<double> distances; // unique sorted distances in deg
    std::vector<double> depths;    // unique sorted depths in km

    // 2D grids of values, indexed by [distance_idx][depth_idx]
    std::vector<std::vector<double>> travelTimes;
    std::vector<std::vector<double>> takeoffDips;
    std::vector<std::vector<double>> dtdds;
    std::vector<std::vector<double>> dtdhs;
  };

  void findIndexes(const std::string &phaseType,
                   double stationElevation,
                   double distance_deg,
                   double eventDepth,
                   const GridData *&grid,
                   size_t &i_dist1,
                   double &dist_frac,
                   size_t &i_depth1,
                   double &depth_frac) const;

private:
  std::string _gridPath;
  std::string _gridModel;
  // keys phase, elevation
  std::unordered_map<std::string, std::map<double, GridData>> _phaseData;
};

} // namespace TTT
} // namespace HDD

#endif

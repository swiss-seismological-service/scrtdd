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

#ifndef __HDD_TTT_NLLGRID_H__
#define __HDD_TTT_NLLGRID_H__

#include "../3rd-party/lrucache.h"
#include "../nll/grid.h"
#include "../ttt.h"
#include "kdtree.h"

#include <unordered_set>

namespace HDD {

namespace TTT {

class NLLGrid : public HDD::TravelTimeTable
{
public:
  NLLGrid(const std::string &gridPath,
          const std::string &gridModel,
          double maxSearchDistance = 0.1,
          bool swapBytes           = false,
          unsigned maxOpenFiles    = 512);

  void compute(double eventLat,
               double eventLon,
               double eventDepth,
               double stationLat,
               double stationLon,
               double stationElevation,
               const std::string &phaseType,
               double &travelTime,
               double &azimuth,
               double &takeOffAngle,
               double &velocityAtSrc) override;

  double compute(double eventLat,
                 double eventLon,
                 double eventDepth,
                 double stationLat,
                 double stationLon,
                 double stationElevation,
                 const std::string &phaseType) override;

  void closeOpenFiles();

private:
  std::string _gridPath;
  std::string _gridModel;
  bool _swapBytes;
  double _maxSearchDistance; // meters

  std::shared_ptr<NLL::VelGrid> _PVelGrid;
  std::shared_ptr<NLL::VelGrid> _SVelGrid;

  using TimeKDTree = KDTree<std::shared_ptr<NLL::TimeGrid>>;
  TimeKDTree _PTimeKDTree;
  TimeKDTree _STimeKDTree;

  using AngleKDTree = KDTree<std::shared_ptr<NLL::AngleGrid>>;
  AngleKDTree _PAngleKDTree;
  AngleKDTree _SAngleKDTree;

  struct GridHolder
  {
    GridHolder(const std::shared_ptr<NLL::VelGrid> &g) : velGrid(g) {}
    GridHolder(const std::shared_ptr<NLL::TimeGrid> &g) : timeGrid(g) {}
    GridHolder(const std::shared_ptr<NLL::AngleGrid> &g) : angleGrid(g) {}
    std::shared_ptr<NLL::VelGrid> velGrid;
    std::shared_ptr<NLL::TimeGrid> timeGrid;
    std::shared_ptr<NLL::AngleGrid> angleGrid;
  };
  lru_cache<std::string, GridHolder> _openGrids;
};

} // namespace TTT

} // namespace HDD

#endif

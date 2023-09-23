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

#include <unordered_set>

namespace HDD {

namespace TTT {

class NLLGrid : public HDD::TravelTimeTable
{
public:
  NLLGrid(const std::string &velGridPath,
          const std::string &timeGridPath,
          const std::string &angleGridPath,
          bool swapBytes,
          unsigned cacheSize = 255);

  virtual ~NLLGrid() = default;

  void freeResources() override;

  void compute(double eventLat,
               double eventLon,
               double eventDepth,
               const Catalog::Station &station,
               const std::string &phaseType,
               double &travelTime,
               double &azimuth,
               double &takeOffAngle,
               double &velocityAtSrc) override;

  double compute(double eventLat,
                 double eventLon,
                 double eventDepth,
                 const Catalog::Station &station,
                 const std::string &phaseType) override;

  static std::string filePath(const std::string &basePath,
                              const Catalog::Station &station,
                              const std::string &phaseType);

private:
  std::string _velGridPath;
  std::string _timeGridPath;
  std::string _angleGridPath;
  bool _swapBytes;
  lru_cache<std::string, std::shared_ptr<NLL::VelGrid>> _velGrids;
  lru_cache<std::string, std::shared_ptr<NLL::TimeGrid>> _timeGrids;
  lru_cache<std::string, std::shared_ptr<NLL::AngleGrid>> _angleGrids;
  std::unordered_set<std::string> _unloadableGrids;
};

} // namespace TTT

} // namespace HDD

#endif

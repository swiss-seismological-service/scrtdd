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

#ifndef __HDD_TTT_H__
#define __HDD_TTT_H__

#include "catalog.h"

namespace HDD {

class TravelTimeTable
{
public:
  TravelTimeTable()          = default;
  virtual ~TravelTimeTable() = default;

  TravelTimeTable(const TravelTimeTable &other)            = delete;
  TravelTimeTable &operator=(const TravelTimeTable &other) = delete;

  virtual void freeResources() = 0;

  /*
   * The implementation of this interface MUST compute:
   * - travelTime    [seconds]
   * - azimuth       [rad]
   * - backAzimuth   [rad]
   * - takeOffAngle  [rad] 0(down):180(up)
   * - velocityAtSrc [km/sec]
   */
  virtual void compute(double eventLat,
                       double eventLon,
                       double eventDepth,
                       const Catalog::Station &station,
                       const std::string &phaseType,
                       double &travelTime,
                       double &azimuth,
                       double &takeOffAngle,
                       double &velocityAtSrc) = 0;

  void compute(const Catalog::Event &event,
               const Catalog::Station &station,
               const std::string &phaseType,
               double &travelTime,
               double &azimuth,
               double &takeOffAngle,
               double &velocityAtSrc)
  {
    return compute(event.latitude, event.longitude, event.depth, station,
                   phaseType, travelTime, azimuth, takeOffAngle, velocityAtSrc);
  }

  /*
   * The implementation of this interface MUST return:
   * - travelTime [seconds]
   * This method is provided for the common case where only the travel time is
   * required (no angles and no velocity).
   * Depending on the specific TTT implementation this method might be faster
   * since less information is required and we want to take adavantage of that.
   * If that is not the case, this method can be implemented as a call to the
   * similar method with additional information and returning only the
   * travelTime value
   */
  virtual double compute(double eventLat,
                         double eventLon,
                         double eventDepth,
                         const Catalog::Station &station,
                         const std::string &phaseType) = 0;

  double compute(const Catalog::Event &event,
                 const Catalog::Station &station,
                 const std::string &phaseType)
  {
    return compute(event.latitude, event.longitude, event.depth, station,
                   phaseType);
  }

protected:
  /*
   * Utility function to compute straight ray path approximation
   * for takeOffAngles angles when those are not available in the
   * travel time tables
   */
  static void computeApproximatedTakeOffAngles(double eventLat,
                                               double eventLon,
                                               double eventDepth,
                                               const Catalog::Station &station,
                                               const std::string &phaseType,
                                               double *azimuth      = nullptr,
                                               double *takeOffAngle = nullptr);
};

} // namespace HDD

#endif

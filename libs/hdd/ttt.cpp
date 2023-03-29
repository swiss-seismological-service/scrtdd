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

#include "ttt.h"
#include "log.h"
#include "utils.h"

#include <sstream>

using namespace std;

namespace HDD {

void TravelTimeTable::computeApproximatedTakeOffAngles(
    double eventLat,
    double eventLon,
    double eventDepth,
    const Catalog::Station &station,
    const std::string &phaseType,
    double *azimuth,
    double *takeOffAngle)
{

  if (!azimuth && !takeOffAngle)
  {
    return;
  }

  double distance =
      computeDistance(eventLat, eventLon, eventDepth, station.latitude,
                      station.longitude, -(station.elevation / 1000.), azimuth);

  if (takeOffAngle)
  {
    double VertDist = eventDepth + station.elevation / 1000.;
    *takeOffAngle   = std::asin(VertDist / distance);
    *takeOffAngle += degToRad(90); // -90(down):+90(up) -> 0(down):180(up)
  }
}

} // namespace HDD

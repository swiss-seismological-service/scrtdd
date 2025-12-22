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

#include "hdd/utils.h"
#include "hdd/catalog.h"
#include "hdd/ttt/nllgrid.h"
#include "hdd/ttt/homogeneous.h"

#include <vector>
#include <string>

namespace {

struct TTTParams
{
  std::string type;
  std::vector<std::string> args;
};

const std::vector<TTTParams> tttList = {
    {"Native:homogeneous",  {"5.8", "3.36"}},
    {"Native:NonLinLoc", {"./data/nll/iasp91_2D_simple/time", "iasp91"}},
    {"Native:NonLinLoc", {"./data/nll/iasp91_2D_sdc/time", "iasp91"}},
    {"Native:NonLinLoc", {"./data/nll/iasp91_2D_global/time", "iasp91"}},
    {"Native:NonLinLoc", {"./data/nll/iasp91_2D_azimuthal_equidist/time", "iasp91"}},
    {"Native:NonLinLoc", {"./data/nll/iasp91_2D_merc/time", "iasp91"}},
    {"Native:NonLinLoc", {"./data/nll/iasp91_2D_lambert/time", "iasp91"}},
    {"Native:NonLinLoc", {"./data/nll/iasp91_3D_simple/time", "iasp91"}},
    {"Native:NonLinLoc", {"./data/nll/iasp91_3D_sdc/time", "iasp91"}},
    {"Native:NonLinLoc", {"./data/nll/iasp91_3D_azimuthal_equidist/time", "iasp91"}},
    {"Native:NonLinLoc", {"./data/nll/iasp91_3D_merc/time", "iasp91"}},
    {"Native:NonLinLoc", {"./data/nll/iasp91_3D_lambert/time", "iasp91"}}
};

std::unique_ptr<HDD::TravelTimeTable> createTTT(const TTTParams &prms)
{
  if (prms.type == "Native:NonLinLoc")
  {
    return std::unique_ptr<HDD::TravelTimeTable>(
        new HDD::TTT::NLLGrid(prms.args[0], prms.args[1], 1.0, false, 256));
  }
  else if (prms.type == "Native:homogeneous")
  {
    double pVel = std::stod(prms.args[0]);
    double sVel = std::stod(prms.args[1]);
    return std::unique_ptr<HDD::TravelTimeTable>(
        new HDD::TTT::Homogeneous(pVel, sVel));
  }
  return nullptr;
}

// Those station parameters must be consistent with nonlinloc
// generated grids files
const std::vector<HDD::Catalog::Station> stationList = {
    {"N1.ST01A", 47.1, 8.6,   500, "N1", "ST01A", ""},
    {"N1.ST02A", 47.1, 8.4,  1000, "N1", "ST02A", ""},
    {"N1.ST03A", 46.9, 8.4,  1500, "N1", "ST03A", ""},
    {"N1.ST04A", 46.9, 8.6,  2000, "N1", "ST04A", ""},
    {"N1.ST01B", 47.0, 8.7, -2000, "N1", "ST01B", ""},
    {"N1.ST02B", 47.0, 8.3, -1500, "N1", "ST02B", ""},
    {"N1.ST03B", 47.2, 8.5, -1000, "N1", "ST03B", ""},
    {"N1.ST04B", 46.8, 8.5,  -500, "N1", "ST04B", ""},
};

const std::vector<HDD::Catalog::Station> stationListZeroElevation = {
    {"N2.ST01A", 47.1, 8.6,     0, "N2", "ST01A", ""},
    {"N2.ST02A", 47.1, 8.4,     0, "N2", "ST02A", ""},
    {"N2.ST03A", 46.9, 8.4,     0, "N2", "ST03A", ""},
    {"N2.ST04A", 46.9, 8.6,     0, "N2", "ST04A", ""},
    {"N2.ST01B", 47.0, 8.7,     0, "N2", "ST01B", ""},
    {"N2.ST02B", 47.0, 8.3,     0, "N2", "ST02B", ""},
    {"N2.ST03B", 47.2, 8.5,     0, "N2", "ST03B", ""},
    {"N2.ST04B", 46.8, 8.5,     0, "N2", "ST04B", ""}
};

// This is in the origin of the nll grid files (TRANS statement)
const double nllGridCentroidLat = 47.0;
const double nllGridCentroidLon = 8.5;

} // namespace

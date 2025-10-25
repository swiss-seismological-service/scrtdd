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

#include "xcorrcache.h"
#include "csvreader.h"
#include "utils.h"
#include <fstream>
#include <vector>

using namespace std;

namespace HDD {

void XCorrCache::writeToFile(const Catalog &cat, const std::string &file) const
{
  ofstream os(file);
  os << "eventId1,eventId2,networkCode,stationCode,locationCode,"
        "phaseType,valid,coefficient,lag"
     << endl;

  auto callback = [&os, &cat](unsigned ev1, unsigned ev2,
                              const std::string &stationId,
                              const Catalog::Phase::Type &type,
                              const XCorrCache::Entry &e) {
    const Catalog::Station &sta = cat.getStations().at(stationId);

    os << strf("%u,%u,%s,%s,%s,%c,%s,%f,%f", ev1, ev2, sta.networkCode.c_str(),
               sta.stationCode.c_str(), sta.locationCode.c_str(),
               static_cast<char>(type), e.valid ? "true" : "false", e.coeff,
               e.lag)
       << endl;
  };

  forEach(callback);
}

XCorrCache XCorrCache::readFromFile(const Catalog &cat, const std::string &file)
{
  auto strToBool = [](const std::string &s) -> bool {
    return s == "1" || s == "true" || s == "True" || s == "TRUE";
  };
  auto strToPhaseType = [](const std::string &s) -> Catalog::Phase::Type {
    return (s == "P" || s == "p") ? Catalog::Phase::Type::P
                                  : Catalog::Phase::Type::S;
  };

  XCorrCache xcorr;
  int row_count = 0;
  try
  {
    vector<unordered_map<string, string>> lines = CSV::readWithHeader(file);
    for (const auto &row : lines)
    {
      unsigned ev1              = std::stoul(row.at("eventId1"));
      unsigned ev2              = std::stoul(row.at("eventId2"));
      string networkCode        = row.at("networkCode");
      string stationCode        = row.at("stationCode");
      string locationCode       = row.at("locationCode");
      Catalog::Phase::Type type = strToPhaseType(row.at("phaseType"));
      bool valid                = strToBool(row.at("valid"));
      double coeff              = std::stod(row.at("coefficient"));
      double lag                = std::stod(row.at("lag"));

      std::string stationId =
          cat.searchStation(networkCode, stationCode, locationCode)->second.id;

      if (!xcorr.has(ev1, ev2, stationId, type))
      {
        xcorr.add(ev1, ev2, stationId, type, valid, coeff, lag);
      }
    }
  }
  catch (std::exception &e)
  {
    string msg = strf("Error while parsing file '%s' at row %d: %s",
                      file.c_str(), row_count, e.what());
    throw Exception(msg);
  }

  return xcorr;
}

} // namespace HDD

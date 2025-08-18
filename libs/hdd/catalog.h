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

#ifndef __HDD_CATALOG_H__
#define __HDD_CATALOG_H__

#include <map>
#include <memory>
#include <unordered_map>
#include <vector>

#include "utctime.h"

namespace HDD {

// DD background catalog
class Catalog
{

public:
  struct Station
  {
    std::string id;
    double latitude;
    double longitude;
    double elevation; // meter
    std::string networkCode;
    std::string stationCode;
    std::string locationCode;

    // Compare attributes when the id is not known (works between multiple
    // catalogs).
    bool operator==(const Station &other) const
    {
      return (networkCode == other.networkCode) &&
             (stationCode == other.stationCode) &&
             (locationCode == other.locationCode) &&
             (latitude == other.latitude) && (longitude == other.longitude) &&
             (elevation == other.elevation);
    }
    bool operator!=(const Station &other) const { return !operator==(other); }
    operator std::string() const { return id; }
  };

  struct Event
  {
    unsigned id; // makes it unique in the catalog
    UTCTime time;
    double latitude;
    double longitude;
    double depth; // km
    double magnitude;

    struct
    {
      bool isRelocated = false;
      double startRms;
      double finalRms;
      unsigned numNeighbours;
      unsigned usedP;
      unsigned usedS;

      struct
      {
        unsigned numTTp;
        unsigned numTTs;
        unsigned numCCp;
        unsigned numCCs;
        double startResidualMedian;
        double startResidualMAD;
        double finalResidualMedian;
        double finalResidualMAD;
      } dd; // double-difference

    } relocInfo;

    // Compare attributes when the id is not known (works between multiple
    // catalogs).
    bool operator==(const Event &other) const
    {
      return (time == other.time) && (latitude == other.latitude) &&
             (longitude == other.longitude) && (depth == other.depth) &&
             (magnitude == other.magnitude);
    }
    bool operator!=(const Event &other) const { return !operator==(other); }
    operator std::string() const { return std::to_string(id); }
  };

  struct Phase
  {
    unsigned eventId;
    std::string stationId;
    UTCTime time;
    double lowerUncertainty;
    double upperUncertainty;
    std::string type;
    std::string networkCode;
    std::string stationCode;
    std::string locationCode;
    std::string channelCode;

    enum class Type : char
    {
      P = 'P',
      S = 'S'
    };

    enum class Source
    {
      CATALOG,
      RT_EVENT_MANUAL,
      RT_EVENT_AUTOMATIC,
    };

    struct
    {
      Type type;
      double classWeight; // 0-1 interval
      Source source;
    } procInfo;

    struct
    {
      bool isRelocated = false;
      double weight;
      double residual;
    } relocInfo;

    bool trusted() const
    {
      return procInfo.source == Source::CATALOG ||
             procInfo.source == Source::RT_EVENT_MANUAL;
    }

    // Compare attributes when the id is not known (works between multiple
    // catalogs).
    bool operator==(const Phase &other) const
    {
      return (time == other.time) &&
             (lowerUncertainty == other.lowerUncertainty) &&
             (upperUncertainty == other.upperUncertainty) &&
             (type == other.type) && (networkCode == other.networkCode) &&
             (stationCode == other.stationCode) &&
             (locationCode == other.locationCode) &&
             (channelCode == other.channelCode);
    }
    bool operator!=(const Phase &other) const { return !operator==(other); }
    operator std::string() const
    {
      return type + "@" + networkCode + "." + stationCode + "." + locationCode +
             "." + channelCode + ":" + UTCClock::toString(time) + ":evId-" +
             std::to_string(eventId);
    }
  };

  Catalog() = default;

  Catalog(std::unordered_map<std::string, Station> &&stations,
          std::map<unsigned, Event> &&events,
          std::unordered_multimap<unsigned, Phase> &&phases)
      : _stations(stations), _events(events), _phases(phases)
  {}

  Catalog(const std::unordered_map<std::string, Station> &stations,
          const std::map<unsigned, Event> &events,
          const std::unordered_multimap<unsigned, Phase> &phases)
      : _stations(stations), _events(events), _phases(phases)
  {}

  Catalog(const std::string &stationFile,
          const std::string &eventFile,
          const std::string &phaseFile,
          bool loadRelocationInfo = false);

  void add(const Catalog &other, bool keepEvId);
  unsigned add(unsigned evId, const Catalog &eventCatalog, bool keepEvId);
  std::unique_ptr<Catalog> extractEvent(unsigned eventId, bool keepEvId) const;

  void removeEvent(unsigned eventId);
  void removePhase(unsigned eventId,
                   const std::string &stationId,
                   const Phase::Type &type);

  std::string addStation(const Station &);
  unsigned addEvent(const Event &);
  void addPhase(const Phase &);

  bool updateStation(const Station &newStation, bool addIfMissing = false);
  bool updateEvent(const Event &newEv, bool addIfMissing = false);
  bool updatePhase(const Phase &newPh, bool addIfMissing = false);

  const std::unordered_map<std::string, Station> &getStations() const
  {
    return _stations;
  }
  const std::map<unsigned, Event> &getEvents() const { return _events; }
  const std::unordered_multimap<unsigned, Phase> &getPhases() const
  {
    return _phases;
  }

  std::unordered_map<std::string, Station>::const_iterator
  searchStation(const std::string &networkCode,
                const std::string &stationCode,
                const std::string &locationCode) const;
  std::map<unsigned, Event>::const_iterator searchEvent(const Event &) const;
  std::unordered_map<unsigned, Phase>::const_iterator
  searchPhase(unsigned eventId,
              const std::string &stationId,
              const Phase::Type &type) const;

  void writeToFile(const std::string &eventFile,
                   const std::string &phaseFile,
                   const std::string &stationFile) const;

  //
  //  static
  //
  static double
  computePickWeight(const Catalog::Phase &phase,
                    const std::vector<double> &uncertaintyClasses);
  static Catalog
  filterPhasesAndSetWeights(const Catalog &catalog,
                            const Catalog::Phase::Source &source,
                            const std::vector<std::string> &PphaseToKeep,
                            const std::vector<std::string> &SphaseToKeep,
                            const std::vector<double> &uncertaintyClasses);

private:
  std::unordered_map<std::string, Station> _stations; // indexed by station id
  std::map<unsigned, Event> _events;                  // indexed by event id
  std::unordered_multimap<unsigned, Phase> _phases;   // indexed by event id
};

} // namespace HDD

#endif

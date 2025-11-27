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

#include "catalog.h"
#include "csvreader.h"
#include "log.h"
#include "utils.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace HDD::Logger;

namespace {

template <typename Map, class Val>
typename Map::iterator searchByValue(Map &SearchMap, const Val &SearchVal)
{
  typename Map::iterator iRet = SearchMap.end();
  for (typename Map::iterator iTer = SearchMap.begin(); iTer != SearchMap.end();
       iTer++)
  {
    if (iTer->second == SearchVal)
    {
      iRet = iTer;
      break;
    }
  }
  return iRet;
}

template <typename Map, class Val>
typename Map::const_iterator searchByValue(const Map &SearchMap,
                                           const Val &SearchVal)
{
  typename Map::const_iterator iRet = SearchMap.end();
  for (typename Map::const_iterator iTer = SearchMap.begin();
       iTer != SearchMap.end(); iTer++)
  {
    if (iTer->second == SearchVal)
    {
      iRet = iTer;
      break;
    }
  }
  return iRet;
}

bool strToBool(const std::string &s)
{
  return s == "1" || s == "true" || s == "True" || s == "TRUE";
}

} // namespace

namespace HDD {

Catalog::Catalog(const string &stationFile,
                 const string &eventFile,
                 const string &phaFile,
                 bool loadRelocationInfo)
{
  if (!pathExists(stationFile))
  {
    string msg = "File " + stationFile + " does not exist";
    throw Exception(msg);
  }

  if (!pathExists(eventFile))
  {
    string msg = "File " + eventFile + " does not exist";
    throw Exception(msg);
  }

  if (!pathExists(phaFile))
  {
    string msg = "File " + phaFile + " does not exist";
    throw Exception(msg);
  }

  unsigned row_count = 0;
  try
  {
    vector<unordered_map<string, string>> stations =
        CSV::readWithHeader(stationFile);
    for (const auto &row : stations)
    {
      row_count++;
      Station sta;
      sta.latitude     = std::stod(row.at("latitude"));
      sta.longitude    = std::stod(row.at("longitude"));
      sta.elevation    = std::stod(row.at("elevation"));
      sta.networkCode  = row.at("networkCode");
      sta.stationCode  = row.at("stationCode");
      sta.locationCode = row.at("locationCode");
      sta.id = sta.networkCode + "." + sta.stationCode + "." + sta.locationCode;
      _stations[sta.id] = sta;
    }
  }
  catch (std::exception &e)
  {
    string msg = strf("Error while parsing file '%s' at row %d: %s",
                      stationFile.c_str(), row_count, e.what());
    throw Exception(msg);
  }

  row_count = 0;
  try
  {
    vector<unordered_map<string, string>> events =
        CSV::readWithHeader(eventFile);
    for (const auto &row : events)
    {
      row_count++;
      Event ev;
      ev.id        = std::stoul(row.at("id"));
      ev.time      = UTCClock::fromString(row.at("isotime"));
      ev.latitude  = std::stod(row.at("latitude"));
      ev.longitude = std::stod(row.at("longitude"));
      ev.depth     = std::stod(row.at("depth"));
      try
      {
        ev.magnitude = std::stod(row.at("magnitude"));
      }
      catch (std::exception &e)
      {
        ev.magnitude = numeric_limits<double>::quiet_NaN();
      }
      ev.relocInfo.isRelocated = false;
      if (loadRelocationInfo && (row.count("relocated") != 0) &&
          strToBool(row.at("relocated")))
      {
        ev.relocInfo.isRelocated = true;
        ev.relocInfo.startRms    = std::stod(row.at("startRms"));
        ev.relocInfo.finalRms    = std::stod(row.at("finalRms"));
      }
      _events[ev.id] = ev;
    }
  }
  catch (std::exception &e)
  {
    string msg = strf("Error while parsing file '%s' at row %d: %s",
                      eventFile.c_str(), row_count, e.what());
    throw Exception(msg);
  }

  row_count = 0;
  try
  {
    vector<unordered_map<string, string>> phases = CSV::readWithHeader(phaFile);
    for (const auto &row : phases)
    {
      row_count++;
      Phase ph;
      ph.eventId = std::stoul(row.at("eventId"));
      ph.time    = UTCClock::fromString(row.at("isotime"));
      try
      {
        ph.lowerUncertainty = std::stod(row.at("lowerUncertainty"));
        ph.upperUncertainty = std::stod(row.at("upperUncertainty"));
      }
      catch (std::exception &e)
      {
        ph.lowerUncertainty = numeric_limits<double>::quiet_NaN();
        ph.upperUncertainty = numeric_limits<double>::quiet_NaN();
      }
      ph.type         = row.at("type");
      ph.networkCode  = row.at("networkCode");
      ph.stationCode  = row.at("stationCode");
      ph.locationCode = row.at("locationCode");
      ph.channelCode  = row.at("channelCode");
      ph.stationId =
          ph.networkCode + "." + ph.stationCode + "." + ph.locationCode;
      ph.relocInfo.isRelocated = false;
      if (loadRelocationInfo && (row.count("usedInReloc") != 0) &&
          strToBool(row.at("usedInReloc")))
      {
        ph.relocInfo.isRelocated   = true;
        ph.relocInfo.startResidual = std::stod(row.at("startResidual"));
        ph.relocInfo.finalResidual = std::stod(row.at("finalResidual"));
        ph.relocInfo.weight        = std::stod(row.at("weight"));
      }

      if (_events.find(ph.eventId) == _events.end())
      {
        throw Exception(strf("unknown event %u", ph.eventId));
      }

      if (_stations.find(ph.stationId) == _stations.end())
      {
        throw Exception(strf("unknown station %s", ph.stationId.c_str()));
      }

      _phases.emplace(ph.eventId, ph);
    }
  }
  catch (std::exception &e)
  {
    string msg = strf("Error while parsing file '%s' at row %d: %s",
                      phaFile.c_str(), row_count, e.what());
    throw Exception(msg);
  }
}

bool Catalog::empty() const
{
  return _stations.empty() && _events.empty() && _phases.empty();
}

void Catalog::add(const Catalog &other, bool keepEvId)
{
  for (const auto &kv : other.getEvents())
  {
    const Catalog::Event &event = kv.second;
    if (keepEvId && _events.find(event.id) != _events.end())
    {
      logDebugF("Skipping duplicated event id %u", event.id);
      continue;
    }
    this->add(event.id, other, keepEvId);
  }
}

Catalog Catalog::extractEvent(unsigned eventId, bool keepEvId) const
{
  Catalog eventToExtract{};

  auto search = this->getEvents().find(eventId);
  if (search == this->getEvents().end())
  {
    string msg = strf("Cannot find event id %u in the catalog.", eventId);
    throw Exception(msg);
  }

  const Catalog::Event &event = search->second;
  unsigned newEventId;

  if (keepEvId)
  {
    eventToExtract._events[event.id] = event;
    newEventId                       = event.id;
  }
  else
  {
    newEventId = eventToExtract.addEvent(event);
  }

  auto eqlrng = this->getPhases().equal_range(event.id);
  for (auto it = eqlrng.first; it != eqlrng.second; ++it)
  {
    Catalog::Phase phase = it->second;

    const Catalog::Station &station = _stations.at(phase.stationId);
    eventToExtract.addStation(station);

    phase.eventId = newEventId;
    eventToExtract.addPhase(phase);
  }

  return eventToExtract;
}

unsigned Catalog::add(unsigned evId, const Catalog &evCat, bool keepEvId)
{
  unsigned newEventId;

  const Catalog::Event &event = evCat._events.find(evId)->second;

  if (keepEvId)
  {
    if (_events.find(event.id) != _events.end())
    {
      throw Exception("Cannot add event, the same id exists already");
    }
    _events[event.id] = event;
    newEventId        = event.id;
  }
  else
  {
    // don't add an event with same values, but keep merging phases
    newEventId = addEvent(event);
  }

  auto eqlrng = evCat._phases.equal_range(event.id);
  for (auto it = eqlrng.first; it != eqlrng.second; ++it)
  {
    Catalog::Phase phase = it->second;

    const Catalog::Station &station = evCat._stations.at(phase.stationId);
    addStation(station);

    phase.eventId = newEventId;
    addPhase(phase);
  }

  return newEventId;
}

void Catalog::removeEvent(unsigned eventId)
{
  map<unsigned, Catalog::Event>::const_iterator it = _events.find(eventId);
  if (it != _events.end())
  {
    _events.erase(it);
  }
  auto eqlrng = _phases.equal_range(eventId);
  _phases.erase(eqlrng.first, eqlrng.second);
}

void Catalog::removePhase(unsigned eventId,
                          const std::string &stationId,
                          const std::string &type)
{
  unordered_map<unsigned, Phase>::const_iterator it =
      searchPhase(eventId, stationId, type);
  if (it != _phases.end())
  {
    _phases.erase(it);
  }
}

bool Catalog::updateStation(const Station &newStation, bool addIfMissing)
{
  unordered_map<string, Catalog::Station>::iterator it =
      _stations.find(newStation.id);
  if (it != _stations.end())
  {
    it->second = newStation;
    return true;
  }
  else if (addIfMissing)
  {
    addStation(newStation);
  }
  return false;
}

bool Catalog::updateEvent(const Event &newEv, bool addIfMissing)
{
  map<unsigned, Catalog::Event>::iterator it = _events.find(newEv.id);
  if (it != _events.end())
  {
    it->second = newEv;
    return true;
  }
  else if (addIfMissing)
  {
    addEvent(newEv);
  }
  return false;
}

bool Catalog::updatePhase(const Phase &newPh, bool addIfMissing)
{
  auto eqlrng = _phases.equal_range(newPh.eventId);
  for (auto it = eqlrng.first; it != eqlrng.second; ++it)
  {
    Catalog::Phase &ph = it->second;
    if (ph.stationId == newPh.stationId &&
        ph.procInfo.type == newPh.procInfo.type)
    {
      ph = newPh;
      return true;
    }
  }

  if (addIfMissing)
  {
    addPhase(newPh);
  }
  return false;
}

map<unsigned, Catalog::Event>::const_iterator
Catalog::searchEvent(const Event &event) const
{
  return searchByValue(_events, event);
}

unordered_map<std::string, Catalog::Station>::const_iterator
Catalog::searchStation(const std::string &networkCode,
                       const std::string &stationCode,
                       const std::string &locationCode) const
{
  string stationId = networkCode + "." + stationCode + "." + locationCode;
  return _stations.find(stationId);
}

unordered_map<unsigned, Catalog::Phase>::const_iterator
Catalog::searchPhase(unsigned eventId,
                     const std::string &stationId,
                     const std::string &type) const
{
  auto eqlrng = _phases.equal_range(eventId);
  for (auto it = eqlrng.first; it != eqlrng.second; ++it)
  {
    const Catalog::Phase &ph = it->second;
    if (ph.stationId == stationId && ph.type == type) return it;
  }
  return _phases.end();
}

string Catalog::addStation(const Station &sta)
{
  string stationId =
      sta.networkCode + "." + sta.stationCode + "." + sta.locationCode;
  if (_stations.find(stationId) == _stations.end())
  {
    Station newSta       = sta;
    newSta.id            = stationId;
    _stations[newSta.id] = newSta;
  }
  return stationId;
}

unsigned Catalog::addEvent(const Event &event)
{
  decltype(_events)::key_type maxKey =
      _events.empty() ? 0 : _events.rbegin()->first;
  Event newEvent       = event;
  newEvent.id          = maxKey + 1;
  _events[newEvent.id] = newEvent;
  return newEvent.id;
}

void Catalog::addPhase(const Phase &phase)
{
  _phases.emplace(phase.eventId, phase);
}

void Catalog::writeToFile(const string &eventFile,
                          const string &phaseFile,
                          const string &stationFile) const
{
  /*
   * write events
   */
  stringstream evStreamNoReloc;
  stringstream evStreamReloc;

  evStreamNoReloc << "id,isotime,latitude,longitude,depth,magnitude";
  evStreamReloc << evStreamNoReloc.str() << ",relocated,startRms,finalRms"
                << endl;
  evStreamNoReloc << endl;

  bool relocInfo = false;
  for (const auto &kv : _events)
  {
    const Catalog::Event &ev = kv.second;

    stringstream evStream;
    evStream << strf("%u,%s,%.12f,%.12f,%g,%.2f", ev.id,
                     UTCClock::toString(ev.time).c_str(), ev.latitude,
                     ev.longitude, ev.depth, ev.magnitude);

    evStreamNoReloc << evStream.str() << endl;
    evStreamReloc << evStream.str();

    if (!ev.relocInfo.isRelocated)
    {
      evStreamReloc << ",false,,";
    }
    else
    {
      relocInfo = true;
      evStreamReloc << strf(",true,%g,%g", ev.relocInfo.startRms,
                            ev.relocInfo.finalRms);
    }
    evStreamReloc << endl;
  }

  ofstream evStream(eventFile);
  evStream << (relocInfo ? evStreamReloc.str() : evStreamNoReloc.str());

  /*
   * write phases
   */
  if (!phaseFile.empty())
  {
    ofstream phStream(phaseFile);

    phStream << "eventId,isotime,lowerUncertainty,upperUncertainty,"
                "type,networkCode,stationCode,locationCode,channelCode";
    if (relocInfo)
    {
      phStream << ",usedInReloc,weight,startResidual,finalResidual";
    }
    phStream << endl;

    const multimap<unsigned, Catalog::Phase> orderedPhases(_phases.begin(),
                                                           _phases.end());
    for (const auto &kv : orderedPhases)
    {
      const Catalog::Phase &ph = kv.second;
      phStream << strf("%u,%s,%g,%g,%s,%s,%s,%s,%s", ph.eventId,
                       UTCClock::toString(ph.time).c_str(), ph.lowerUncertainty,
                       ph.upperUncertainty, ph.type.c_str(),
                       ph.networkCode.c_str(), ph.stationCode.c_str(),
                       ph.locationCode.c_str(), ph.channelCode.c_str());

      if (relocInfo)
      {
        if (!ph.relocInfo.isRelocated)
        {
          phStream << ",false,,,";
        }
        else
        {
          phStream << strf(",true,%g,%g,%g", ph.relocInfo.weight,
                           ph.relocInfo.startResidual,
                           ph.relocInfo.finalResidual);
        }
      }
      phStream << endl;
    }
  }

  /*
   * write stations
   */
  if (!stationFile.empty())
  {
    ofstream staStream(stationFile);
    staStream
        << "latitude,longitude,elevation,networkCode,stationCode,locationCode"
        << endl;

    const map<string, Catalog::Station> orderedStations(_stations.begin(),
                                                        _stations.end());
    for (const auto &kv : orderedStations)
    {
      const Catalog::Station &sta = kv.second;
      staStream << strf("%.12f,%.12f,%g,%s,%s,%s", sta.latitude, sta.longitude,
                        sta.elevation, sta.networkCode.c_str(),
                        sta.stationCode.c_str(), sta.locationCode.c_str())
                << endl;
    }
  }
}

/*
 * Build a catalog with requested phases, only. Besides, make sure, that for an
 * event/station pair there is only one P and one S phase. If multiple phases
 * are found, keep the one with the highest priority.
 */
Catalog
Catalog::filterPhasesAndSetWeights(const Catalog &catalog,
                                   const Phase::Source &source,
                                   const std::vector<std::string> &PphaseToKeep,
                                   const std::vector<std::string> &SphaseToKeep,
                                   const vector<double> &uncertaintyClasses)
{
  unordered_multimap<unsigned, Phase> filteredS;
  unordered_multimap<unsigned, Phase> filteredP;

  // loop through each event
  for (const auto &kv : catalog.getEvents())
  {
    const Event &event = kv.second;

    // loop through all phases of current event
    const auto &eqlrng = catalog.getPhases().equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
      // keep the phase only if it has a higher priority than an existing one
      // or if this is the only one for a given station
      const Phase &phase = it->second;

      // P phase
      auto itpp = find(PphaseToKeep.begin(), PphaseToKeep.end(), phase.type);
      if (itpp != PphaseToKeep.end())
      {
        auto priority = std::distance(PphaseToKeep.begin(), itpp);

        // fetch already selected P phases for current event, and check if
        // there is already a P phase for the same station
        bool inserted = false;
        auto eqlrng2  = filteredP.equal_range(event.id);
        for (auto it2 = eqlrng2.first; it2 != eqlrng2.second; ++it2)
        {
          Phase &existingPhase = it2->second;
          auto existingPriority =
              std::distance(PphaseToKeep.begin(),
                            find(PphaseToKeep.begin(), PphaseToKeep.end(),
                                 existingPhase.type));
          if (existingPhase.type == phase.type &&
              existingPhase.stationId == phase.stationId &&
              existingPriority < priority)
          {
            existingPhase = phase;
            inserted      = true;
            break;
          }
        }
        if (!inserted) filteredP.emplace(phase.eventId, phase);
        continue;
      }

      // S phase
      auto itsp = find(SphaseToKeep.begin(), SphaseToKeep.end(), phase.type);
      if (itsp != SphaseToKeep.end())
      {
        auto priority = std::distance(SphaseToKeep.begin(), itsp);

        // fetch already selected S phases for current event, and check if
        // there is already a S phase for the same station
        bool inserted = false;
        auto eqlrng2  = filteredS.equal_range(event.id);
        for (auto it2 = eqlrng2.first; it2 != eqlrng2.second; ++it2)
        {
          Phase &existingPhase = it2->second;
          auto existingPriority =
              std::distance(SphaseToKeep.begin(),
                            find(SphaseToKeep.begin(), SphaseToKeep.end(),
                                 existingPhase.type));
          if (existingPhase.type == phase.type &&
              existingPhase.stationId == phase.stationId &&
              existingPriority < priority)
          {
            existingPhase = phase;
            inserted      = true;
            break;
          }
        }
        if (!inserted) filteredS.emplace(phase.eventId, phase);
        continue;
      }

      logDebugF("Discard phase (%s), the type is not among the selected ones",
                string(phase).c_str());
    }
  }

  // loop through selected phases and replace the actual phase name with a
  // generic P or S
  unordered_multimap<unsigned, Phase> filteredPhases;
  for (auto &it : filteredP)
  {
    Phase &phase               = it.second;
    phase.procInfo.classWeight = computePickWeight(phase, uncertaintyClasses);
    phase.procInfo.type        = Phase::Type::P;
    phase.procInfo.source      = source;

    filteredPhases.emplace(phase.eventId, phase);
  }
  for (auto &it : filteredS)
  {
    Phase &phase               = it.second;
    phase.procInfo.classWeight = computePickWeight(phase, uncertaintyClasses);
    phase.procInfo.type        = Phase::Type::S;
    phase.procInfo.source      = source;

    filteredPhases.emplace(phase.eventId, phase);
  }

  return Catalog(catalog.getStations(), catalog.getEvents(), filteredPhases);
}

double Catalog::computePickWeight(const Phase &phase,
                                  const vector<double> &uncertaintyClasses)
{
  if (uncertaintyClasses.size() < 2)
  {
    return 1.;
  }

  double uncertainty =
      (phase.lowerUncertainty + phase.upperUncertainty) / 2.; // secs

  // set lowest class as default
  unsigned uncertaintyClass = uncertaintyClasses.size() - 1;

  if (uncertainty >= 0 && std::isfinite(uncertainty) &&
      uncertaintyClasses.size() > 1 && uncertainty < uncertaintyClasses.back())
  {

    for (unsigned curr = 0, next = 1; next < uncertaintyClasses.size();
         curr++, next++)
    {
      if (uncertainty >= uncertaintyClasses.at(curr) &&
          uncertainty <= uncertaintyClasses.at(next))
      {
        uncertaintyClass = curr;
        break;
      }
    }
  }

  return 1. / pow(2, uncertaintyClass);
}

} // namespace HDD

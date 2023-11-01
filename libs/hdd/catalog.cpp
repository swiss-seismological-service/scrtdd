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
    string msg = strf("Error while parsing file '%s' at row %d (%s)",
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
      ev.id                    = std::stoul(row.at("id"));
      ev.time                  = UTCClock::fromString(row.at("isotime"));
      ev.latitude              = std::stod(row.at("latitude"));
      ev.longitude             = std::stod(row.at("longitude"));
      ev.depth                 = std::stod(row.at("depth"));
      ev.magnitude             = std::stod(row.at("magnitude"));
      ev.relocInfo.isRelocated = false;
      if (loadRelocationInfo && (row.count("relocated") != 0) &&
          strToBool(row.at("relocated")))
      {
        ev.relocInfo.isRelocated   = true;
        ev.relocInfo.startRms      = std::stod(row.at("startRms"));
        ev.relocInfo.finalRms      = std::stod(row.at("finalRms"));
        ev.relocInfo.locChange     = std::stod(row.at("locChange"));
        ev.relocInfo.depthChange   = std::stod(row.at("depthChange"));
        ev.relocInfo.timeChange    = std::stod(row.at("timeChange"));
        ev.relocInfo.numNeighbours = std::stoul(row.at("numNeighbours"));
        ev.relocInfo.phases.usedP  = std::stoul(row.at("ph_usedP"));
        ev.relocInfo.phases.usedS  = std::stoul(row.at("ph_usedS"));
        ev.relocInfo.phases.stationDistMin =
            std::stod(row.at("ph_stationDistMin"));
        ev.relocInfo.phases.stationDistMedian =
            std::stod(row.at("ph_stationDistMedian"));
        ev.relocInfo.phases.stationDistMax =
            std::stod(row.at("ph_stationDistMax"));
        ev.relocInfo.dd.numTTp = std::stoul(row.at("dd_numTTp"));
        ev.relocInfo.dd.numTTs = std::stoul(row.at("dd_numTTs"));
        ev.relocInfo.dd.numCCp = std::stoul(row.at("dd_numCCp"));
        ev.relocInfo.dd.numCCs = std::stoul(row.at("dd_numCCs"));
        ev.relocInfo.dd.startResidualMedian =
            std::stod(row.at("dd_startResidualMedian"));
        ev.relocInfo.dd.startResidualMAD =
            std::stod(row.at("dd_startResidualMAD"));
        ev.relocInfo.dd.finalResidualMedian =
            std::stod(row.at("dd_finalResidualMedian"));
        ev.relocInfo.dd.finalResidualMAD =
            std::stod(row.at("dd_finalResidualMAD"));
      }
      _events[ev.id] = ev;
    }
  }
  catch (std::exception &e)
  {
    string msg = strf("Error while parsing file '%s' at row %d (%s)",
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
      ph.eventId          = std::stoul(row.at("eventId"));
      ph.time             = UTCClock::fromString(row.at("isotime"));
      ph.lowerUncertainty = std::stod(row.at("lowerUncertainty"));
      ph.upperUncertainty = std::stod(row.at("upperUncertainty"));
      ph.type             = row.at("type");
      ph.networkCode      = row.at("networkCode");
      ph.stationCode      = row.at("stationCode");
      ph.locationCode     = row.at("locationCode");
      ph.channelCode      = row.at("channelCode");
      ph.stationId =
          ph.networkCode + "." + ph.stationCode + "." + ph.locationCode;
      ph.isManual              = row.at("evalMode") == "manual";
      ph.relocInfo.isRelocated = false;
      if (loadRelocationInfo && (row.count("usedInReloc") != 0) &&
          strToBool(row.at("usedInReloc")))
      {
        ph.relocInfo.isRelocated     = true;
        ph.procInfo.weight           = std::stod(row.at("startWeight"));
        ph.relocInfo.finalWeight     = std::stod(row.at("finalWeight"));
        ph.relocInfo.startTTResidual = std::stod(row.at("startTTResidual"));
        ph.relocInfo.finalTTResidual = std::stod(row.at("finalTTResidual"));
        ph.relocInfo.numTTObs        = std::stoul(row.at("numTTObs"));
        ph.relocInfo.numCCObs        = std::stoul(row.at("numCCObs"));
        ph.relocInfo.startMeanDDResidual =
            std::stod(row.at("startMeanDDResidual"));
        ph.relocInfo.finalMeanDDResidual =
            std::stod(row.at("finalMeanDDResidual"));
      }
      _phases.emplace(ph.eventId, ph);
    }
  }
  catch (std::exception &e)
  {
    string msg = strf("Error while parsing file '%s' at row %d (%s)",
                      phaFile.c_str(), row_count, e.what());
    throw Exception(msg);
  }
}

void Catalog::add(const Catalog &other, bool keepEvId)
{
  for (const auto &kv : other.getEvents())
  {
    const Catalog::Event &event = kv.second;
    if (keepEvId && _events.find(event.id) != _events.end())
    {
      logDebug("Skipping duplicated event id %u", event.id);
      continue;
    }
    this->add(event.id, other, keepEvId);
  }
}

unique_ptr<Catalog> Catalog::extractEvent(unsigned eventId, bool keepEvId) const
{
  unique_ptr<Catalog> eventToExtract(new Catalog());

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
    eventToExtract->_events[event.id] = event;
    newEventId                        = event.id;
  }
  else
  {
    newEventId = eventToExtract->addEvent(event);
  }

  auto eqlrng = this->getPhases().equal_range(event.id);
  for (auto it = eqlrng.first; it != eqlrng.second; ++it)
  {
    Catalog::Phase phase = it->second;

    const Catalog::Station &station = _stations.at(phase.stationId);
    eventToExtract->addStation(station);

    phase.eventId = newEventId;
    eventToExtract->addPhase(phase);
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
      throw Exception("Cannot add event, internal logic error");
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
                          const Phase::Type &type)
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
                     const Phase::Type &type) const
{
  auto eqlrng = _phases.equal_range(eventId);
  for (auto it = eqlrng.first; it != eqlrng.second; ++it)
  {
    const Catalog::Phase &ph = it->second;
    if (ph.stationId == stationId && ph.procInfo.type == type) return it;
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
  evStreamReloc << evStreamNoReloc.str()
                << ",relocated,startRms,finalRms,locChange,depthChange,"
                   "timeChange,numNeighbours,"
                   "ph_usedP,ph_usedS,ph_stationDistMin,ph_stationDistMedian,"
                   "ph_stationDistMax,"
                   "dd_numTTp,dd_numTTs,dd_numCCp,dd_numCCs,"
                   "dd_startResidualMedian,dd_startResidualMAD,"
                   "dd_finalResidualMedian,dd_finalResidualMAD"
                << endl;
  evStreamNoReloc << endl;

  bool relocInfo = false;
  for (const auto &kv : _events)
  {
    const Catalog::Event &ev = kv.second;

    stringstream evStream;
    evStream << strf("%u,%s,%.12f,%.12f,%.10f,%.2f", ev.id,
                     UTCClock::toString(ev.time).c_str(), ev.latitude,
                     ev.longitude, ev.depth, ev.magnitude);

    evStreamNoReloc << evStream.str() << endl;
    evStreamReloc << evStream.str();

    if (!ev.relocInfo.isRelocated)
    {
      evStreamReloc << ",false,,,,,,,,,,,,,,,,,,,";
    }
    else
    {
      relocInfo = true;
      evStreamReloc << strf(
          ",true,%.3f,%.3f,%.3f,%.3f,%.3f,%u,%u,%u,%.3f,%.3f,%.3f,%u,%u,%u,%u,%"
          ".4f,%.4f,%.4f,%.4f",
          ev.relocInfo.startRms, ev.relocInfo.finalRms, ev.relocInfo.locChange,
          ev.relocInfo.depthChange, ev.relocInfo.timeChange,
          ev.relocInfo.numNeighbours, ev.relocInfo.phases.usedP,
          ev.relocInfo.phases.usedS, ev.relocInfo.phases.stationDistMin,
          ev.relocInfo.phases.stationDistMedian,
          ev.relocInfo.phases.stationDistMax, ev.relocInfo.dd.numTTp,
          ev.relocInfo.dd.numTTs, ev.relocInfo.dd.numCCp,
          ev.relocInfo.dd.numCCs, ev.relocInfo.dd.startResidualMedian,
          ev.relocInfo.dd.startResidualMAD, ev.relocInfo.dd.finalResidualMedian,
          ev.relocInfo.dd.finalResidualMAD);
    }
    evStreamReloc << endl;
  }

  ofstream evStream(eventFile);
  evStream << (relocInfo ? evStreamReloc.str() : evStreamNoReloc.str());

  /*
   * write phases
   */
  ofstream phStream(phaseFile);

  phStream << "eventId,isotime,lowerUncertainty,upperUncertainty,"
              "type,networkCode,stationCode,locationCode,channelCode,evalMode";
  if (relocInfo)
  {
    phStream << ",usedInReloc,startTTResidual,finalTTResidual,"
                "startWeight,finalWeight,numTTObs,numCCObs,"
                "startMeanDDResidual,finalMeanDDResidual";
  }
  phStream << endl;

  const multimap<unsigned, Catalog::Phase> orderedPhases(_phases.begin(),
                                                         _phases.end());
  for (const auto &kv : orderedPhases)
  {
    const Catalog::Phase &ph = kv.second;
    phStream << strf("%u,%s,%.3f,%.3f,%s,%s,%s,%s,%s,%s", ph.eventId,
                     UTCClock::toString(ph.time).c_str(), ph.lowerUncertainty,
                     ph.upperUncertainty, ph.type.c_str(),
                     ph.networkCode.c_str(), ph.stationCode.c_str(),
                     ph.locationCode.c_str(), ph.channelCode.c_str(),
                     (ph.isManual ? "manual" : "automatic"));

    if (relocInfo)
    {
      if (!ph.relocInfo.isRelocated)
      {
        phStream << ",false,,,,,,,,";
      }
      else
      {
        phStream << strf(
            ",true,%.3f,%.3f,%.2f,%.2f,%u,%u,%.4f,%.4f",
            ph.relocInfo.startTTResidual, ph.relocInfo.finalTTResidual,
            ph.procInfo.weight, ph.relocInfo.finalWeight, ph.relocInfo.numTTObs,
            ph.relocInfo.numCCObs, ph.relocInfo.startMeanDDResidual,
            ph.relocInfo.finalMeanDDResidual);
      }
    }
    phStream << endl;
  }

  /*
   * write stations
   */
  ofstream staStream(stationFile);
  staStream
      << "latitude,longitude,elevation,networkCode,stationCode,locationCode"
      << endl;

  const map<string, Catalog::Station> orderedStations(_stations.begin(),
                                                      _stations.end());
  for (const auto &kv : orderedStations)
  {
    const Catalog::Station &sta = kv.second;
    staStream << strf("%.12f,%.12f,%.7f,%s,%s,%s", sta.latitude, sta.longitude,
                      sta.elevation, sta.networkCode.c_str(),
                      sta.stationCode.c_str(), sta.locationCode.c_str())
              << endl;
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
                                   const std::vector<std::string> &SphaseToKeep)
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

      logDebug("Discard phase (%s), the type is not among the selected ones",
               string(phase).c_str());
    }
  }

  // loop through selected phases and replace the actual phase name with a
  // generic P or S
  unordered_multimap<unsigned, Phase> filteredPhases;
  for (auto &it : filteredP)
  {
    Phase &phase          = it.second;
    phase.procInfo.weight = computePickWeight(phase);
    phase.procInfo.type   = Phase::Type::P;
    phase.procInfo.source = source;

    filteredPhases.emplace(phase.eventId, phase);
  }
  for (auto &it : filteredS)
  {
    Phase &phase          = it.second;
    phase.procInfo.weight = computePickWeight(phase);
    phase.procInfo.type   = Phase::Type::S;
    phase.procInfo.source = source;

    filteredPhases.emplace(phase.eventId, phase);
  }

  return Catalog(catalog.getStations(), catalog.getEvents(), filteredPhases);
}

/*
 * Fixed weighting scheme based on pick time uncertainties
 * Class 0: 0     - 0.025  sec
 *       1: 0.025 - 0.050  sec
 *       2: 0.050 - 0.100  sec
 *       3: 0.100 - 0.200  sec
 *       4: 0.200 - 0.400  sec
 *       5: 0.400 -        sec
 *  weight = 1 / 2^class
 */
double Catalog::computePickWeight(double uncertainty /* secs */)
{
  unsigned uncertaintyClass;

  if (uncertainty >= 0.000 && uncertainty <= 0.025)
    uncertaintyClass = 0;
  else if (uncertainty > 0.025 && uncertainty <= 0.050)
    uncertaintyClass = 1;
  else if (uncertainty > 0.050 && uncertainty <= 0.100)
    uncertaintyClass = 2;
  else if (uncertainty > 0.100 && uncertainty <= 0.200)
    uncertaintyClass = 3;
  else if (uncertainty > 0.200 && uncertainty <= 0.400)
    uncertaintyClass = 4;
  else
    uncertaintyClass = 5;

  return 1 / std::pow(2, uncertaintyClass);
}

double Catalog::computePickWeight(const Phase &phase)
{
  return computePickWeight((phase.lowerUncertainty + phase.upperUncertainty) /
                           2.);
}

} // namespace HDD

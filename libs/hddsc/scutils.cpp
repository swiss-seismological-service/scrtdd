/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU Affero General Public License as published*
 * by the Free Software Foundation, either version 3 of the License, or    *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This program is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#include "scutils.h"

#include "hdd/csvreader.h"
#include "hdd/log.h"

#include <fstream>
#include <iostream>

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp/logging/log.h>

#include <seiscomp/datamodel/amplitude.h>
#include <seiscomp/datamodel/arrival.h>
#include <seiscomp/datamodel/event.h>
#include <seiscomp/datamodel/magnitude.h>
#include <seiscomp/datamodel/origin.h>
#include <seiscomp/datamodel/originreference.h>
#include <seiscomp/datamodel/stationmagnitude.h>
#include <seiscomp/datamodel/stationmagnitudecontribution.h>
#include <seiscomp/math/geo.h>
#include <seiscomp/utils/files.h>

using namespace std;
using namespace Seiscomp;
using PhaseSrc       = HDD::Catalog::Phase::Source;
using XCorrEvalStats = HDD::DD::XCorrEvalStats;
using Seiscomp::Core::stringify;

namespace {

double normalizeAz(double az)
{
  if (az < 0)
    az += 360.0;
  else if (az >= 360.0)
    az -= 360.0;
  return az;
}

std::pair<double, double> getPickUncertainty(DataModel::Pick *pick)
{
  pair<double, double> uncertainty(-1, -1); // secs
  try
  {
    // symmetric uncertainty
    uncertainty.first = uncertainty.second = pick->time().uncertainty();
  }
  catch (Core::ValueException &)
  {
    // unsymmetric uncertainty
    try
    {
      uncertainty.first  = pick->time().lowerUncertainty();
      uncertainty.second = pick->time().upperUncertainty();
    }
    catch (Core::ValueException &)
    {}
  }

  if (uncertainty.first < 0 && uncertainty.second < 0)
  {
    try
    {
      uncertainty.first = uncertainty.second =
          (pick->evaluationMode() == Seiscomp::DataModel::MANUAL)
              ? HDD::Catalog::DEFAULT_MANUAL_PICK_UNCERTAINTY
              : HDD::Catalog::DEFAULT_AUTOMATIC_PICK_UNCERTAINTY;
    }
    catch (Core::ValueException &)
    {
      uncertainty.first = uncertainty.second =
          HDD::Catalog::DEFAULT_AUTOMATIC_PICK_UNCERTAINTY;
    }
  }

  return uncertainty;
}

DataModel::SensorLocation *findSensorLocation(const std::string &networkCode,
                                              const std::string &stationCode,
                                              const std::string &locationCode,
                                              const Core::Time &atTime)
{
  DataModel::Inventory *inv = Client::Inventory::Instance()->inventory();
  if (!inv)
  {
    SEISCOMP_DEBUG("Inventory not available");
    return nullptr;
  }

  DataModel::InventoryError error;
  DataModel::SensorLocation *loc = DataModel::getSensorLocation(
      inv, networkCode, stationCode, locationCode, atTime, &error);

  if (!loc)
  {
    SEISCOMP_DEBUG(
        "Unable to fetch SensorLocation information (%s.%s.%s at %s): %s",
        networkCode.c_str(), stationCode.c_str(), locationCode.c_str(),
        atTime.iso().c_str(), error.toString());
  }
  return loc;
}

} // namespace

namespace HDD {
namespace SCAdapter {

DataSource::DataSource(DataModel::DatabaseQuery *query,
                       DataModel::PublicObjectTimeSpanBuffer *cache)
    : _query(query), _cache(cache)
{}

DataSource::DataSource(DataModel::EventParameters *eventParameters)
    : _eventParameters(eventParameters)
{}

DataSource::DataSource(DataModel::DatabaseQuery *query,
                       DataModel::PublicObjectTimeSpanBuffer *cache,
                       DataModel::EventParameters *eventParameters)
    : _query(query), _cache(cache), _eventParameters(eventParameters)
{}

DataModel::PublicObject *DataSource::getObject(const Core::RTTI &classType,
                                               const std::string &publicID)
{
  DataModel::PublicObject *ret = nullptr;

  if (_eventParameters && !ret)
  {
    if (classType == DataModel::Pick::TypeInfo())
      ret = _eventParameters->findPick(publicID);
    else if (classType == DataModel::Amplitude::TypeInfo())
      ret = _eventParameters->findAmplitude(publicID);
    else if (classType == DataModel::Origin::TypeInfo())
      ret = _eventParameters->findOrigin(publicID);
    else if (classType == DataModel::Event::TypeInfo())
      ret = _eventParameters->findEvent(publicID);
  }

  if (_cache && !ret)
  {
    ret = _cache->find(classType, publicID);
  }

  return ret;
}

void DataSource::loadArrivals(DataModel::Origin *org)
{
  if (_query)
  {
    if (org->arrivalCount() == 0) _query->loadArrivals(org);
  }
}

void DataSource::loadMagnitudes(DataModel::Origin *org,
                                bool loadStationMagnitudeContributions,
                                bool loadStationMagnitudes)
{
  if (_query)
  {
    if (org->magnitudeCount() == 0) _query->loadMagnitudes(org);

    if (loadStationMagnitudeContributions)
    {
      for (size_t i = 0; i < org->magnitudeCount(); i++)
      {
        DataModel::Magnitude *mag = org->magnitude(i);
        if (mag->stationMagnitudeContributionCount() == 0)
          _query->loadStationMagnitudeContributions(mag);
      }
    }

    if (loadStationMagnitudes && org->stationMagnitudeCount() == 0)
    {
      _query->loadStationMagnitudes(org);
    }
  }
}

DataModel::Event *DataSource::getParentEvent(const std::string &originID)
{
  DataModel::Event *ret = nullptr;

  if (_eventParameters && !ret)
  {
    for (size_t i = 0; i < _eventParameters->eventCount() && !ret; i++)
    {
      DataModel::Event *ev = _eventParameters->event(i);
      for (size_t j = 0; j < ev->originReferenceCount() && !ret; j++)
      {
        DataModel::OriginReference *orgRef = ev->originReference(j);
        if (orgRef->originID() == originID) ret = ev;
      }
    }
  }

  if (_query && !ret)
  {
    ret = _query->getEvent(originID);
  }

  return ret;
}

std::unordered_map<unsigned, DataModel::OriginPtr>
addToCatalog(HDD::Catalog &cat,
             const std::vector<DataModel::OriginPtr> &origins,
             DataSource &dataSrc)
{
  std::unordered_map<unsigned, DataModel::OriginPtr> idmap;
  for (DataModel::OriginPtr org : origins)
  {
    dataSrc.loadArrivals(org.get());

    if (org->arrivalCount() == 0)
    {
      SEISCOMP_WARNING("Origin %s doesn't have any arrival. Skip it.",
                       org->publicID().c_str());
      continue;
    }

    // add event
    HDD::Catalog::Event ev;
    ev.id        = 0;
    ev.time      = fromSC(org->time().value());
    ev.latitude  = org->latitude();
    ev.longitude = org->longitude();
    ev.depth     = org->depth(); // km

    DataModel::MagnitudePtr mag;
    // try to fetch preferred magnitude of the event
    DataModel::EventPtr parentEvent = dataSrc.getParentEvent(org->publicID());
    if (parentEvent)
    {
      mag = dataSrc.get<DataModel::Magnitude>(
          parentEvent->preferredMagnitudeID());
    }
    if (mag)
    {
      ev.magnitude = mag->magnitude();
    }
    else
    {
      SEISCOMP_DEBUG("Origin %s: cannot load preferred magnitude from parent "
                     "event, set it to 0",
                     org->publicID().c_str());
      ev.magnitude = 0.;
    }

    SEISCOMP_DEBUG("Adding origin '%s' to the Catalog",
                   org->publicID().c_str());

    unsigned newEventId = cat.addEvent(ev);

    // add phases
    for (size_t i = 0; i < org->arrivalCount(); ++i)
    {
      DataModel::Arrival *orgArr    = org->arrival(i);
      const DataModel::Phase &orgPh = orgArr->phase();

      DataModel::PickPtr pick = dataSrc.get<DataModel::Pick>(orgArr->pickID());
      if (!pick)
      {
        SEISCOMP_ERROR("Cannot load pick '%s' (origin %s)",
                       orgArr->pickID().c_str(), org->publicID().c_str());
        continue;
      }

      // find the station
      HDD::Catalog::Station sta;
      sta.networkCode  = pick->waveformID().networkCode();
      sta.stationCode  = pick->waveformID().stationCode();
      sta.locationCode = pick->waveformID().locationCode();

      // skip not selected picks/phases or those which have 0 weight, unless
      // manual
      try
      {
        if (pick->evaluationMode() != DataModel::MANUAL &&
            (orgArr->weight() == 0 || !orgArr->timeUsed()))
        {
          SEISCOMP_DEBUG("Discarding not used %s phase %s.%s",
                         orgPh.code().c_str(), sta.networkCode.c_str(),
                         sta.stationCode.c_str());
          continue;
        }
      }
      catch (Core::ValueException &)
      {}

      // add station if not already there
      if (cat.searchStation(sta.networkCode, sta.stationCode,
                            sta.locationCode) == cat.getStations().end())
      {
        DataModel::SensorLocation *loc = findSensorLocation(
            sta.networkCode, sta.stationCode, sta.locationCode, pick->time());

        if (!loc)
        {
          SEISCOMP_ERROR(
              "Cannot load sensor location %s.%s.%s information for arrival "
              "'%s' (origin '%s'). All picks associated with this station will "
              "not be used.",
              sta.networkCode.c_str(), sta.stationCode.c_str(),
              sta.locationCode.c_str(), orgArr->pickID().c_str(),
              org->publicID().c_str());
          continue;
        }

        sta.latitude  = loc->latitude();
        sta.longitude = loc->longitude();
        sta.elevation = loc->elevation(); // meter
        cat.addStation(sta);
      }
      // the station must be available at this point
      sta =
          cat.searchStation(sta.networkCode, sta.stationCode, sta.locationCode)
              ->second;

      // get uncertainty
      pair<double, double> uncertainty = getPickUncertainty(pick.get());

      HDD::Catalog::Phase ph;
      ph.eventId          = newEventId;
      ph.stationId        = sta.id;
      ph.time             = fromSC(pick->time().value());
      ph.lowerUncertainty = uncertainty.first;
      ph.upperUncertainty = uncertainty.second;
      ph.type             = orgPh.code();
      ph.networkCode      = pick->waveformID().networkCode();
      ph.stationCode      = pick->waveformID().stationCode();
      ph.locationCode     = pick->waveformID().locationCode();
      ph.channelCode      = pick->waveformID().channelCode();
      ph.isManual         = (pick->evaluationMode() == DataModel::MANUAL);
      cat.addPhase(ph);
    }
    idmap[newEventId] = org;
  }
  return idmap;
}

std::unordered_map<unsigned, DataModel::OriginPtr> addToCatalog(
    HDD::Catalog &cat, const std::vector<std::string> &ids, DataSource &dataSrc)
{
  vector<DataModel::OriginPtr> origins;

  for (const string &id : ids)
  {
    DataModel::OriginPtr org = dataSrc.get<DataModel::Origin>(id);
    if (!org)
    {
      SEISCOMP_ERROR("Cannot find origin with id %s", id.c_str());
      continue;
    }
    origins.push_back(org);
  }

  return addToCatalog(cat, origins, dataSrc);
}

std::unordered_map<unsigned, DataModel::OriginPtr>
addToCatalog(HDD::Catalog &cat, const std::string &idFile, DataSource &dataSrc)
{
  if (!Util::fileExists(idFile))
  {
    string msg = "File " + idFile + " does not exist";
    throw std::runtime_error(msg);
  }

  SEISCOMP_INFO("Reading file %s which must contain at least a column "
                "with header 'origin' and an origin id per line",
                idFile.c_str());

  vector<string> ids;
  vector<unordered_map<string, string>> rows = HDD::CSV::readWithHeader(idFile);

  for (const auto &row : rows)
  {
    const string &id = row.at("origin");
    ids.push_back(id);
  }

  return addToCatalog(cat, ids, dataSrc);
}

void convertOrigin(DataSource &dataSrc,
                   const HDD::Catalog &relocatedOrg,
                   DataModel::Origin *org, // can be nullptr
                   const string &author,
                   const string &agencyID,
                   const string &methodID,
                   const string &earthModelID,
                   bool includeMagnitude,
                   bool includeExistingPicks,
                   DataModel::OriginPtr &newOrg,                 // return value
                   std::vector<DataModel::PickPtr> &newOrgPicks) // return value
{
  // there must be only one event in the catalog, the relocated origin
  const HDD::Catalog::Event &event = relocatedOrg.getEvents().begin()->second;

  newOrg = DataModel::Origin::Create();

  DataModel::CreationInfo ci;
  ci.setAgencyID(agencyID);
  ci.setAuthor(author);
  ci.setCreationTime(Core::Time::GMT());

  newOrg->setCreationInfo(ci);
  newOrg->setEarthModelID(earthModelID);
  newOrg->setMethodID(methodID);
  newOrg->setEvaluationMode(DataModel::EvaluationMode(DataModel::AUTOMATIC));

  newOrg->setTime(DataModel::TimeQuantity(toSC(event.time)));

  DataModel::RealQuantity latitude = DataModel::RealQuantity(event.latitude);
  newOrg->setLatitude(event.latitude);

  DataModel::RealQuantity longitude =
      DataModel::RealQuantity(normalizeLon(event.longitude));
  newOrg->setLongitude(longitude);

  DataModel::RealQuantity depth = DataModel::RealQuantity(event.depth);
  newOrg->setDepth(depth);

  if (event.relocInfo.isRelocated)
  {
    DataModel::CommentPtr comment = new DataModel::Comment();
    comment->setId("relocation::report");
    comment->setText(HDD::DD::relocationReport(relocatedOrg));
    comment->setCreationInfo(ci);
    newOrg->add(comment.get());
  }

  auto evPhases      = relocatedOrg.getPhases().equal_range(event.id);
  int usedPhaseCount = 0;
  vector<double> azi;
  vector<double> staDistances;
  set<string> associatedStations;
  set<string> usedStations;

  // If we know the origin before relocation fetch some information from it
  if (org)
  {
    //
    // keep track of the triggering origin of this relocation
    //
    DataModel::CommentPtr comment = new DataModel::Comment();
    comment->setId("relocation::sourceOrigin");
    comment->setText(org->publicID());
    comment->setCreationInfo(ci);
    newOrg->add(comment.get());
    //
    // Copy magnitude from org if that is Manual
    //
    if (includeMagnitude)
    {
      dataSrc.loadMagnitudes(org, true, true);

      unordered_map<string, string> staMagIdMap;
      for (size_t i = 0; i < org->stationMagnitudeCount(); i++)
      {
        DataModel::StationMagnitude *staMag = org->stationMagnitude(i);
        DataModel::StationMagnitude *newStaMag =
            DataModel::StationMagnitude::Create();
        *newStaMag = *staMag;
        newOrg->add(newStaMag);
        staMagIdMap[staMag->publicID()] = newStaMag->publicID();
      }

      for (size_t i = 0; i < org->magnitudeCount(); i++)
      {
        DataModel::Magnitude *mag    = org->magnitude(i);
        DataModel::Magnitude *newMag = DataModel::Magnitude::Create();
        *newMag                      = *mag;

        for (size_t j = 0; j < mag->stationMagnitudeContributionCount(); j++)
        {
          DataModel::StationMagnitudeContribution *contrib =
              mag->stationMagnitudeContribution(j);
          DataModel::StationMagnitudeContributionPtr newContrib =
              new DataModel::StationMagnitudeContribution(*contrib);
          try
          {
            newContrib->setStationMagnitudeID(
                staMagIdMap.at(contrib->stationMagnitudeID()));
          }
          catch (...)
          {}
          newMag->add(newContrib.get());
        }
        newOrg->add(newMag);
      }
    }

    //
    // add all arrivals that were in the original Origin (before relocation)
    //
    for (size_t i = 0; i < org->arrivalCount(); i++)
    {
      DataModel::Arrival *orgArr = org->arrival(i);

      DataModel::ArrivalPtr newArr = new DataModel::Arrival();
      newArr->setPickID(orgArr->pickID());
      newArr->setPhase(orgArr->phase());
      newArr->setWeight(0.);
      newArr->setTimeUsed(false);

      newOrg->add(newArr.get());

      DataModel::PickPtr pick = dataSrc.get<DataModel::Pick>(orgArr->pickID());
      if (pick)
      {
        associatedStations.insert(pick->waveformID().networkCode() + "." +
                                  pick->waveformID().stationCode());

        if (includeExistingPicks)
          newOrgPicks.push_back(DataModel::Pick::Cast(pick->clone()));
      }
    }
  }

  // add missing arrivals and fill in all the properties
  for (auto it = evPhases.first; it != evPhases.second; ++it)
  {
    const HDD::Catalog::Phase &phase = it->second;
    bool phaseUsed =
        phase.relocInfo.isRelocated && phase.relocInfo.finalWeight != 0;

    // drop phases discovered via cross-correlation if those phases were not
    // used for the relocations
    if ((phase.procInfo.source == PhaseSrc::THEORETICAL ||
         phase.procInfo.source == PhaseSrc::XCORR) &&
        !phaseUsed)
    {
      continue;
    }

    associatedStations.insert(phase.networkCode + "." + phase.stationCode);

    // check if this phase has been already added
    bool alreadyAdded = false;
    DataModel::Arrival *newArr;

    for (size_t i = 0; i < newOrg->arrivalCount(); i++)
    {
      newArr                  = newOrg->arrival(i);
      DataModel::PickPtr pick = dataSrc.get<DataModel::Pick>(newArr->pickID());

      if (pick && toSC(phase.time) == pick->time().value() &&
          phase.networkCode == pick->waveformID().networkCode() &&
          phase.stationCode == pick->waveformID().stationCode() &&
          phase.locationCode == pick->waveformID().locationCode() &&
          phase.channelCode == pick->waveformID().channelCode())
      {
        alreadyAdded = true;
        break;
      }
    }

    if (!alreadyAdded)
    {
      // prepare the new pick
      DataModel::PickPtr newPick = DataModel::Pick::Create();
      newPick->setCreationInfo(ci);
      newPick->setMethodID(methodID);
      newPick->setEvaluationMode(
          phase.isManual ? DataModel::EvaluationMode(DataModel::MANUAL)
                         : DataModel::EvaluationMode(DataModel::AUTOMATIC));
      DataModel::TimeQuantity pickTime(toSC(phase.time));
      pickTime.setLowerUncertainty(phase.lowerUncertainty);
      pickTime.setUpperUncertainty(phase.upperUncertainty);
      newPick->setTime(pickTime);
      newPick->setPhaseHint(DataModel::Phase(phase.type));
      newPick->setWaveformID(DataModel::WaveformStreamID(
          phase.networkCode, phase.stationCode, phase.locationCode,
          phase.channelCode, ""));
      newOrgPicks.push_back(newPick);

      // prepare the new arrival
      newArr = new DataModel::Arrival();
      newArr->setCreationInfo(ci);
      newArr->setPickID(newPick->publicID());
      newArr->setPhase(phase.type);

      newOrg->add(newArr);
    }

    newArr->setWeight(phase.relocInfo.isRelocated ? phase.relocInfo.finalWeight
                                                  : 0.);
    newArr->setTimeUsed(phaseUsed);
    newArr->setTimeResidual(
        phase.relocInfo.isRelocated ? phase.relocInfo.finalTTResidual : 0.);

    auto search = relocatedOrg.getStations().find(phase.stationId);
    if (search == relocatedOrg.getStations().end())
    {
      SEISCOMP_WARNING("Cannot find station id '%s' referenced by phase '%s'."
                       "Cannot add Arrival to relocated origin",
                       phase.stationId.c_str(), string(phase).c_str());
      continue;
    }
    const HDD::Catalog::Station &station = search->second;

    double distance, az, baz;
    Math::Geo::delazi(event.latitude, event.longitude, station.latitude,
                      station.longitude, &distance, &az, &baz);

    newArr->setAzimuth(normalizeAz(az));
    newArr->setDistance(distance);

    // update stats
    if (newArr->timeUsed())
    {
      usedPhaseCount++;
      staDistances.push_back(distance);
      azi.push_back(az);
      usedStations.insert(phase.stationId);
    }
  }

  // finish computing stats
  double primaryAz = 360., secondaryAz = 360.;
  if (azi.size() >= 2)
  {
    primaryAz = secondaryAz = 0.;
    sort(azi.begin(), azi.end());
    vector<double>::size_type aziCount = azi.size();
    azi.push_back(azi[0] + 360.);
    azi.push_back(azi[1] + 360.);
    for (vector<double>::size_type i = 0; i < aziCount; i++)
    {
      double gap = azi[i + 1] - azi[i];
      if (gap > primaryAz) primaryAz = gap;
      gap = azi[i + 2] - azi[i];
      if (gap > secondaryAz) secondaryAz = gap;
    }
  }

  // add quality
  DataModel::OriginQuality oq;
  oq.setAssociatedPhaseCount(newOrg->arrivalCount());
  oq.setUsedPhaseCount(usedPhaseCount);
  oq.setAssociatedStationCount(associatedStations.size());
  oq.setUsedStationCount(usedStations.size());
  oq.setStandardError(event.relocInfo.isRelocated ? event.relocInfo.finalRms
                                                  : 0);
  oq.setMedianDistance(computeMedian(staDistances));
  oq.setMinimumDistance(*min_element(staDistances.begin(), staDistances.end()));
  oq.setMaximumDistance(*max_element(staDistances.begin(), staDistances.end()));
  oq.setAzimuthalGap(primaryAz);
  oq.setSecondaryAzimuthalGap(secondaryAz);
  newOrg->setQuality(oq);
}

void printEvalXcorrStats(
    const XCorrEvalStats &pTotStats,
    const XCorrEvalStats &sTotStats,
    const map<string, XCorrEvalStats> &pStatsByStation,
    const map<string, XCorrEvalStats> &sStatsByStation,
    const map<unsigned, XCorrEvalStats> &pStatsByStaDistance,
    const map<unsigned, XCorrEvalStats> &sStatsByStaDistance,
    const map<unsigned, XCorrEvalStats> &pStatsByInterEvDistance,
    const map<unsigned, XCorrEvalStats> &sStatsByInterEvDistance,
    double interEvDistStep,
    double staDistStep,
    double completionPercent)
{
  unsigned skipped, performed;
  double meanCoeff, meanCoeffAbsDev, medianCoeff, medianCoeffAbsDev;
  double meanLag, meanLagAbsDev, medianLag, medianLagAbsDev;

  string log;
  if (completionPercent >= 1)
    log += "---FINAL STATS--\n";
  else
    log += stringify("---PARTIAL STATS %.1f%%---\n", completionPercent * 100);

  auto header = [](const string &firstColumn) {
    return stringify(
        "%-20s   #CC   #Skip meanCC (MAD) medianCC (MAD) meanLag (MAD) "
        "medianLag (MAD)\n",
        firstColumn.c_str());
  };

  auto row = [&](const string &firstColumn, const XCorrEvalStats &stats) {
    stats.summarize(skipped, performed, meanCoeff, meanCoeffAbsDev, medianCoeff,
                    medianCoeffAbsDev, meanLag, meanLagAbsDev, medianLag,
                    medianLagAbsDev);
    return stringify(
        "%-15s %10u %7u %5.2f (%3.2f) %7.2f (%3.2f) %7.f (%3.f) %9.f (%3.f)\n",
        firstColumn.c_str(), performed, skipped, meanCoeff, meanCoeffAbsDev,
        medianCoeff, medianCoeffAbsDev, meanLag * 1000, meanLagAbsDev * 1000,
        medianLag * 1000, medianLagAbsDev * 1000);
  };

  log += header("Total xcorr P phases");
  log += row("", pTotStats);
  log += header("Total xcorr S phases");
  log += row("", sTotStats);

  log += stringify("Xcorr P phases by inter-event distance in %.2f km step\n",
                   interEvDistStep);
  log += header(" EvDist [km]");
  for (const auto &kv : pStatsByInterEvDistance)
    log += row(stringify("%5.2f-%-5.2f", kv.first * interEvDistStep,
                         (kv.first + 1) * interEvDistStep),
               kv.second);

  log += stringify("Xcorr S phases by inter-event distance in %.2f km step\n",
                   interEvDistStep);
  log += header(" EvDist [km]");
  for (const auto &kv : sStatsByInterEvDistance)
    log += row(stringify("%5.2f-%-5.2f", kv.first * interEvDistStep,
                         (kv.first + 1) * interEvDistStep),
               kv.second);

  log +=
      stringify("XCorr P phases by event to station distance in %.2f km step\n",
                staDistStep);
  log += header("StaDist [km]");
  for (const auto &kv : pStatsByStaDistance)
    log += row(stringify("%3d-%-3d", int(kv.first * staDistStep),
                         int((kv.first + 1) * staDistStep)),
               kv.second);

  log +=
      stringify("XCorr S phases by event to station distance in %.2f km step\n",
                staDistStep);
  log += header("StaDist [km]");
  for (const auto &kv : sStatsByStaDistance)
    log += row(stringify("%3d-%-3d", int(kv.first * staDistStep),
                         int((kv.first + 1) * staDistStep)),
               kv.second);

  log += stringify("XCorr P phases by station\n");
  log += header("Station");
  for (const auto &kv : pStatsByStation)
    log += row(stringify("%-12s", kv.first.c_str()), kv.second);

  log += stringify("XCorr S phases by station\n");
  log += header("Station");
  for (const auto &kv : sStatsByStation)
    log += row(stringify("%-12s", kv.first.c_str()), kv.second);

  std::cout << log;
}
} // namespace SCAdapter
} // namespace HDD

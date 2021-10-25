/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as          *
 * published by the Free Software Foundation, either version 3 of the      *
 * License, or (at your option) any later version.                         *
 *                                                                         *
 * This software is distributed in the hope that it will be useful,        *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#include "sccatalog.h"
#include "csvreader.h"

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <seiscomp3/client/inventory.h>
#include <seiscomp3/core/datetime.h>
#include <seiscomp3/core/strings.h>
#include <seiscomp3/datamodel/amplitude.h>
#include <seiscomp3/datamodel/event.h>
#include <seiscomp3/datamodel/magnitude.h>
#include <seiscomp3/datamodel/origin.h>
#include <seiscomp3/datamodel/pick.h>
#include <seiscomp3/utils/files.h>
#include <stdexcept>

#define SEISCOMP_COMPONENT HDD
#include <seiscomp3/logging/log.h>

using namespace std;
using namespace Seiscomp;
using Seiscomp::Core::stringify;

namespace {

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

} // namespace

namespace Seiscomp {
namespace HDD {

DataModel::PublicObject *
ScCatalog::DataSource::getObject(const Core::RTTI &classType,
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

void ScCatalog::DataSource::loadArrivals(DataModel::Origin *org)
{
  if (_query)
  {
    if (org->arrivalCount() == 0) _query->loadArrivals(org);
  }
}

void ScCatalog::DataSource::loadMagnitudes(
    DataModel::Origin *org,
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

DataModel::Event *
ScCatalog::DataSource::getParentEvent(const std::string &originID)
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
ScCatalog::add(const std::vector<DataModel::OriginPtr> &origins,
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
    Event ev;
    ev.id        = 0;
    ev.time      = org->time().value();
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

    SEISCOMP_DEBUG("Adding origin '%s' to the ScCatalog",
                   org->publicID().c_str());

    unsigned newEventId = this->addEvent(ev);

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
      Station sta;
      sta.networkCode  = pick->waveformID().networkCode();
      sta.stationCode  = pick->waveformID().stationCode();
      sta.locationCode = pick->waveformID().locationCode();

      // skip not selected picks/phases or those which have 0 weight, unless
      // manual
      try
      {
        if (pick->evaluationMode() != Seiscomp::DataModel::MANUAL &&
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
      if (searchStation(sta.networkCode, sta.stationCode, sta.locationCode) ==
          _stations.end())
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
        this->addStation(sta);
      }
      // the station must be available at this point
      sta = searchStation(sta.networkCode, sta.stationCode, sta.locationCode)
                ->second;

      // get uncertainty
      pair<double, double> uncertainty = getPickUncertainty(pick.get());

      Phase ph;
      ph.eventId          = newEventId;
      ph.stationId        = sta.id;
      ph.time             = pick->time().value();
      ph.lowerUncertainty = uncertainty.first;
      ph.upperUncertainty = uncertainty.second;
      ph.type             = orgPh.code();
      ph.networkCode      = pick->waveformID().networkCode();
      ph.stationCode      = pick->waveformID().stationCode();
      ph.locationCode     = pick->waveformID().locationCode();
      ph.channelCode      = pick->waveformID().channelCode();
      ph.isManual = (pick->evaluationMode() == Seiscomp::DataModel::MANUAL);
      this->addPhase(ph);
    }
    idmap[newEventId] = org;
  }
  return idmap;
}

std::unordered_map<unsigned, DataModel::OriginPtr>
ScCatalog::add(const std::vector<std::string> &ids, DataSource &dataSrc)
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

  return add(origins, dataSrc);
}

std::unordered_map<unsigned, DataModel::OriginPtr>
ScCatalog::add(const std::string &idFile, DataSource &dataSrc)
{
  if (!Util::fileExists(idFile))
  {
    string msg = "File " + idFile + " does not exist";
    throw runtime_error(msg);
  }

  vector<string> ids;
  vector<unordered_map<string, string>> rows = HDD::CSV::readWithHeader(idFile);

  for (const auto &row : rows)
  {
    const string &id = row.at("seiscompId");
    ids.push_back(id);
  }

  return add(ids, dataSrc);
}

/*
 * static methods
 */

DataModel::SensorLocation *
ScCatalog::findSensorLocation(const std::string &networkCode,
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

} // namespace HDD
} // namespace Seiscomp

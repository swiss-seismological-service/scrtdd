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

#ifndef __HDD_SCCATALOG_H__
#define __HDD_SCCATALOG_H__

#include <seiscomp3/core/baseobject.h>
#include <seiscomp3/datamodel/databasequery.h>
#include <seiscomp3/datamodel/eventparameters.h>
#include <seiscomp3/datamodel/origin.h>
#include <seiscomp3/datamodel/publicobjectcache.h>
#include <seiscomp3/datamodel/station.h>

#include "catalog.h"
#include <unordered_map>
#include <vector>

namespace Seiscomp {
namespace HDD {

DEFINE_SMARTPOINTER(ScCatalog);

class ScCatalog : public HDD::Catalog
{

public:
  class DataSource
  {
  public:
    DataSource(DataModel::DatabaseQuery *query,
               DataModel::PublicObjectTimeSpanBuffer *cache)
        : _query(query), _cache(cache)
    {}

    DataSource(DataModel::EventParameters *eventParameters)
        : _eventParameters(eventParameters)
    {}

    DataSource(DataModel::DatabaseQuery *query,
               DataModel::PublicObjectTimeSpanBuffer *cache,
               DataModel::EventParameters *eventParameters)
        : _query(query), _cache(cache), _eventParameters(eventParameters)
    {}

    template <typename T>
    typename Core::SmartPointer<T>::Impl get(const std::string &publicID)
    {
      return T::Cast(getObject(T::TypeInfo(), publicID));
    }

    DataModel::PublicObject *getObject(const Seiscomp::Core::RTTI &classType,
                                       const std::string &publicID);

    void loadMagnitudes(DataModel::Origin *org,
                        bool loadStationMagnitudeContributions,
                        bool loadStationMagnitudes);

    void loadArrivals(DataModel::Origin *org);

    DataModel::Event *getParentEvent(const std::string &originID);

  private:
    DataModel::DatabaseQuery *_query;
    DataModel::PublicObjectTimeSpanBuffer *_cache;
    DataModel::EventParameters *_eventParameters;
  };

  ScCatalog()          = default;
  virtual ~ScCatalog() = default;

  ScCatalog(const ScCatalog &other) = default;
  ScCatalog &operator=(const ScCatalog &other) = default;

  ScCatalog(ScCatalog &&other) = default;
  ScCatalog &operator=(ScCatalog &&other) = default;

  ScCatalog(std::unordered_map<std::string, Catalog::Station> &&stations,
            std::map<unsigned, Catalog::Event> &&events,
            std::unordered_multimap<unsigned, Catalog::Phase> &&phases)
      : Catalog(stations, events, phases)
  {}
  ScCatalog(const std::unordered_map<std::string, Catalog::Station> &stations,
            const std::map<unsigned, Catalog::Event> &events,
            const std::unordered_multimap<unsigned, Catalog::Phase> &phases)
      : Catalog(stations, events, phases)
  {}
  ScCatalog(const std::string &stationFile,
            const std::string &eventFile,
            const std::string &phaseFile,
            bool loadRelocationInfo = false)
      : Catalog(stationFile, eventFile, phaseFile)
  {}

  // populate from seiscomp data format
  std::unordered_map<unsigned, DataModel::OriginPtr>
  add(const std::vector<DataModel::OriginPtr> &origins, DataSource &dataSrc);
  std::unordered_map<unsigned, DataModel::OriginPtr>
  add(const std::vector<std::string> &ids, DataSource &dataSrc);
  std::unordered_map<unsigned, DataModel::OriginPtr>
  add(const std::string &idFile, DataSource &dataSrc);

  //
  // static
  //
  static DataModel::SensorLocation *
  findSensorLocation(const std::string &networkCode,
                     const std::string &stationCode,
                     const std::string &locationCode,
                     const Core::Time &atTime);
};

} // namespace HDD
} // namespace Seiscomp

#endif

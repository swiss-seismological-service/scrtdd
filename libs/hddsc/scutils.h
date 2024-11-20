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

#ifndef __HDD_SCUTILS_H__
#define __HDD_SCUTILS_H__

#include "hdd/catalog.h"
#include "hdd/dd.h"
#include "hdd/timewindow.h"
#include "hdd/utctime.h"

#include <seiscomp/client/inventory.h>
#include <seiscomp/core/datetime.h>
#include <seiscomp/core/genericrecord.h>
#include <seiscomp/core/timewindow.h>
#include <seiscomp/core/typedarray.h>
#include <seiscomp/core/version.h>
#include <seiscomp/datamodel/databasequery.h>
#include <seiscomp/datamodel/eventparameters.h>
#include <seiscomp/datamodel/publicobjectcache.h>

namespace HDD {
namespace SCAdapter {

class DataSource
{
public:
  DataSource(Seiscomp::DataModel::DatabaseQuery *query,
             Seiscomp::DataModel::PublicObjectTimeSpanBuffer *cache);

  DataSource(Seiscomp::DataModel::EventParameters *eventParameters);

  DataSource(Seiscomp::DataModel::DatabaseQuery *query,
             Seiscomp::DataModel::PublicObjectTimeSpanBuffer *cache,
             Seiscomp::DataModel::EventParameters *eventParameters);

  template <typename T>
  typename Seiscomp::Core::SmartPointer<T>::Impl
  get(const std::string &publicID)
  {
    return T::Cast(getObject(T::TypeInfo(), publicID));
  }

  Seiscomp::DataModel::PublicObject *
  getObject(const Seiscomp::Core::RTTI &classType, const std::string &publicID);

  void loadArrivals(Seiscomp::DataModel::Origin *org);

  void loadMagnitudes(Seiscomp::DataModel::Origin *org,
                      bool loadStationMagnitudeContributions,
                      bool loadStationMagnitudes);

  Seiscomp::DataModel::Event *getParentEvent(const std::string &originID);

private:
  Seiscomp::DataModel::DatabaseQuery *_query;
  Seiscomp::DataModel::PublicObjectTimeSpanBuffer *_cache;
  Seiscomp::DataModel::EventParameters *_eventParameters;
};

std::unordered_map<unsigned, Seiscomp::DataModel::OriginPtr>
addToCatalog(HDD::Catalog &cat,
             const std::vector<Seiscomp::DataModel::OriginPtr> &origins,
             DataSource &dataSrc);

std::unordered_map<unsigned, Seiscomp::DataModel::OriginPtr>
addToCatalog(HDD::Catalog &cat,
             const std::vector<std::string> &ids,
             DataSource &dataSrc);

std::unordered_map<unsigned, Seiscomp::DataModel::OriginPtr>
addToCatalog(HDD::Catalog &cat, const std::string &idFile, DataSource &dataSrc);

void convertOrigin(DataSource &dataSrc,
                   const HDD::Catalog &relocatedOrg,
                   Seiscomp::DataModel::Origin *org,
                   const std::string &author,
                   const std::string &agencyID,
                   const std::string &methodID,
                   const std::string &earthModelID,
                   bool includeMagnitude,
                   bool includeExistingPicks,
                   Seiscomp::DataModel::OriginPtr &newOrg,
                   std::vector<Seiscomp::DataModel::PickPtr> &newOrgPicks);

void printEvalXcorrStats(
    const HDD::DD::XCorrEvalStats &pTotStats,
    const HDD::DD::XCorrEvalStats &sTotStats,
    const std::map<std::string, HDD::DD::XCorrEvalStats> &pStatsByStation,
    const std::map<std::string, HDD::DD::XCorrEvalStats> &sStatsByStation,
    const std::map<unsigned, HDD::DD::XCorrEvalStats> &pStatsByStaDistance,
    const std::map<unsigned, HDD::DD::XCorrEvalStats> &sStatsByStaDistance,
    const std::map<unsigned, HDD::DD::XCorrEvalStats> &pStatsByInterEvDistance,
    const std::map<unsigned, HDD::DD::XCorrEvalStats> &sStatsByInterEvDistance,
    double interEvDistStep,
    double staDistStep,
    double completionPercent);

inline Seiscomp::Core::Time toSC(const HDD::UTCTime &t)
{
  return Seiscomp::Core::Time(HDD::durToSec(t.time_since_epoch()));
}

inline HDD::UTCTime fromSC(const Seiscomp::Core::Time &t)
{
#if SC_API_VERSION < SC_API_VERSION_CHECK(17,0,0)
	return HDD::UTCTime() + HDD::secToDur(t.length());
#else
	return HDD::UTCTime() + HDD::secToDur(t.epoch());
#endif
}

inline Seiscomp::Core::TimeWindow toSC(const HDD::TimeWindow &tw)
{
  return Seiscomp::Core::TimeWindow(toSC(tw.startTime()), toSC(tw.endTime()));
}

inline HDD::TimeWindow fromSC(const Seiscomp::Core::TimeWindow &tw)
{
  return HDD::TimeWindow(fromSC(tw.startTime()), fromSC(tw.endTime()));
}

} // namespace SCAdapter
} // namespace HDD
#endif

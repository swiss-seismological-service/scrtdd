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

#ifndef __HDD_DATASRC_H__
#define __HDD_DATASRC_H__

#include <seiscomp3/core/baseobject.h>
#include <seiscomp3/datamodel/eventparameters.h>
#include <seiscomp3/datamodel/publicobjectcache.h>
#include <seiscomp3/datamodel/databasequery.h>
#include <seiscomp3/datamodel/origin.h>

namespace Seiscomp {
namespace HDD {


class DataSource {
    public:

        DataSource(DataModel::DatabaseQuery* query,
                   DataModel::PublicObjectTimeSpanBuffer* cache)
        : _query(query), _cache(cache) {}

        DataSource(DataModel::EventParameters* eventParameters)
        : _eventParameters(eventParameters) {}

        DataSource(DataModel::DatabaseQuery* query,
                   DataModel::PublicObjectTimeSpanBuffer* cache,
                   DataModel::EventParameters* eventParameters)
        : _query(query), _cache(cache), _eventParameters(eventParameters) {}

        template <typename T>
        typename Core::SmartPointer<T>::Impl
        get(const std::string& publicID) {
            return T::Cast(getObject(T::TypeInfo(), publicID));
        }

        DataModel::PublicObject* getObject(const Seiscomp::Core::RTTI& classType,
                                           const std::string& publicID);

        void loadArrivals(DataModel::Origin* org);

        DataModel::Event* getParentEvent(const std::string& originID);

    private:
        DataModel::DatabaseQuery* _query;
        DataModel::PublicObjectTimeSpanBuffer* _cache;
        DataModel::EventParameters* _eventParameters;
};



}
}

#endif


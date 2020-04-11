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
                                   

#include "datasrc.h"

#include <seiscomp3/datamodel/pick.h>
#include <seiscomp3/datamodel/amplitude.h>
#include <seiscomp3/datamodel/event.h>
#include <seiscomp3/datamodel/origin.h>

namespace Seiscomp {
namespace HDD {

DataModel::PublicObject* DataSource::getObject(const Core::RTTI& classType,
                                                const std::string& publicID)
{
    DataModel::PublicObject* ret = nullptr;

    if ( _eventParameters && ! ret)
    {
        if (classType == DataModel::Pick::TypeInfo() )
            ret = _eventParameters->findPick(publicID);
        else if (classType == DataModel::Amplitude::TypeInfo() )
            ret = _eventParameters->findAmplitude(publicID);
        else if (classType == DataModel::Origin::TypeInfo() )
            ret = _eventParameters->findOrigin(publicID);
        else if (classType == DataModel::Event::TypeInfo() )
            ret = _eventParameters->findEvent(publicID);
    }

    if ( _cache && ! ret)
    {
        ret = _cache->find(classType, publicID);
    }

   return ret;
}



void DataSource::loadArrivals(DataModel::Origin* org)
{
    bool found = false;

    if ( _eventParameters && ! found)
    {
        DataModel::Origin* epOrg = _eventParameters->findOrigin(org->publicID());
        if ( epOrg )
        {
            found = true;
            for (size_t i = 0; i < epOrg->arrivalCount(); i++)
                org->add(DataModel::Arrival::Cast(epOrg->arrival(i)->clone()));
        }
    }

    if ( _query && ! found )
    {
        found = true;
        _query->loadArrivals(org);
    }
}



DataModel::Event* DataSource::getParentEvent(const std::string& originID)
{
    DataModel::Event* ret = nullptr;

    if ( _eventParameters && ! ret )
    {
        for (size_t i = 0; i <  _eventParameters->eventCount() && !ret; i++)
        {
            DataModel::Event* ev = _eventParameters->event(i);
            for (size_t j = 0; j < ev->originReferenceCount() && !ret; j++)
            {
                DataModel::OriginReference* orgRef = ev->originReference(j);
                if (orgRef->originID() == originID)
                    ret = ev;
            }
        }
    }

    if ( _query && ! ret)
    {
        ret = _query->getEvent(originID);
    }

    return ret;

}



}
}


/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 *   You can redistribute and/or modify this program under the             *
 *   terms of the "SED Public License for Seiscomp Contributions"          *
 *                                                                         *
 *   You should have received a copy of the "SED Public License for        *
 *   Seiscomp Contributions" with this. If not, you can find it at         *
 *   http://www.seismo.ethz.ch/static/seiscomp_contrib/license.txt         *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   "SED Public License for Seiscomp Contributions" for more details      *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/


#include "catalog.h"
#include "csvreader.h"

#include <seiscomp3/core/datetime.h>
#include <seiscomp3/core/strings.h>
#include <seiscomp3/client/inventory.h>
#include <seiscomp3/datamodel/station.h>
#include <seiscomp3/datamodel/event.h>
#include <seiscomp3/datamodel/pick.h>
#include <seiscomp3/datamodel/origin.h>
#include <seiscomp3/datamodel/magnitude.h>
#include <seiscomp3/datamodel/amplitude.h>
#include <seiscomp3/utils/files.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <boost/filesystem.hpp>
#include <boost/bind.hpp>
#include <boost/range/iterator_range_core.hpp>

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp3/logging/log.h>
 
using namespace std;
using namespace Seiscomp;
using Seiscomp::Core::stringify;

namespace {


template <class Map, class Val> typename Map::iterator
searchByValue(Map & SearchMap, const Val & SearchVal)
{
    typename Map::iterator iRet = SearchMap.end();
    for (typename Map::iterator iTer = SearchMap.begin(); iTer != SearchMap.end(); iTer ++)
    {
        if (iTer->second == SearchVal)
        {
            iRet = iTer;
            break;
        }
    }
    return iRet;
}



template <class Map, class Val> typename Map::const_iterator
searchByValue(const Map & SearchMap, const Val & SearchVal)
{
    typename Map::const_iterator iRet = SearchMap.end();
    for (typename Map::const_iterator iTer = SearchMap.begin(); iTer != SearchMap.end(); iTer ++)
    {
        if (iTer->second == SearchVal)
        {
            iRet = iTer;
            break;
        }
    }
    return iRet;
}



std::pair<double,double> getPickUncertainty(DataModel::Pick *pick)
{
    pair<double,double> uncertainty(-1, -1); // secs
    try
    {
        // Symmetric uncertainty
        uncertainty.first = uncertainty.second = pick->time().uncertainty();
    }
    catch ( Core::ValueException& )
    {
        // Unsymmetric uncertainty
        try
        {
          uncertainty.first = pick->time().lowerUncertainty();
          uncertainty.second = pick->time().upperUncertainty();
        } catch ( Core::ValueException& ) {}
    }

    if ( uncertainty.first < 0 && uncertainty.second < 0 )
    {
        try {
           uncertainty.first = uncertainty.second = (pick->evaluationMode() == Seiscomp::DataModel::MANUAL)
                                                  ? HDD::Catalog::DEFAULT_MANUAL_PICK_UNCERTAINTY
                                                  : HDD::Catalog::DEFAULT_AUTOMATIC_PICK_UNCERTAINTY;
        } catch ( Core::ValueException& )  {
           uncertainty.first = uncertainty.second = HDD::Catalog::DEFAULT_AUTOMATIC_PICK_UNCERTAINTY;
        }
    }

    return uncertainty;
}



bool strToBool(const std::string& s) {
  return s == "1" || s == "true" || s == "True" || s == "TRUE";
}


}


namespace Seiscomp {
namespace HDD {


/*
 * static methods
 */
 
 
DataModel::Station* Catalog::findStation(const string& netCode,
                                         const string& stationCode,
                                         const Core::Time& atTime)
{
    DataModel::Inventory *inv = Client::Inventory::Instance()->inventory();
    if ( ! inv )
    {
        SEISCOMP_ERROR("Inventory not available");
        return nullptr;
    }

    DataModel::InventoryError error;
    DataModel::Station* station = DataModel::getStation(inv, netCode, stationCode, atTime, &error);
    if ( ! station ) 
    {
        SEISCOMP_ERROR("Cannot find station %s.%s: %s",
             netCode.c_str(), stationCode.c_str(), error.toString() );
        return nullptr;
    }
    return station;
}


/*
 * DataSource class
 */

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


/*
 * Catalog class
 */



Catalog::Catalog() : Catalog(map<string,Station>(),
                             map<unsigned,Event>(),
                             multimap<unsigned,Phase>())
{
}



Catalog::Catalog(const map<string,Station>& stations,
                 const map<unsigned,Event>& events,
                 const multimap<unsigned,Phase>& phases)
        : _stations(stations),
          _events(events),
          _phases(phases) 
{
}


Catalog::Catalog(map<string,Station>&& stations,
                 map<unsigned,Event>&& events,
                 multimap<unsigned,Phase>&& phases)
        : _stations(stations),
          _events(events),
          _phases(phases) 
{
}



Catalog::Catalog(const Catalog& other) : Catalog(other._stations,
                                                 other._events,
                                                 other._phases)
{ }



Catalog::Catalog(Catalog&& other) : Catalog(std::move(other._stations),
                                            std::move(other._events),
                                            std::move(other._phases) )
{ }


Catalog& Catalog::operator=(const Catalog& other)
{
    if (this == &other) return *this;
    _stations = other._stations;
    _events   = other._events;
    _phases   = other._phases;
    return *this;
}



Catalog& Catalog::operator=(Catalog&& other)
{
    if (this == &other) return *this;
    _stations = std::move(other._stations);
    _events   = std::move(other._events);
    _phases   = std::move(other._phases);
    return *this;
}



Catalog::Catalog(const string& stationFile,
                 const string& eventFile,
                 const string& phaFile,
                 bool loadRelocationInfo)
{
    if ( !Util::fileExists(stationFile) )
    {
        string msg = "File " + stationFile + " does not exist";
        throw runtime_error(msg);
    }

    if ( !Util::fileExists(eventFile) )
    {
        string msg = "File " + eventFile + " does not exist";
        throw runtime_error(msg);
    }

    if ( !Util::fileExists(phaFile) )
    {
        string msg = "File " + phaFile + " does not exist";
        throw runtime_error(msg);
    }

    vector<map<string,string> > stations = CSV::readWithHeader(stationFile);

    for (const auto& row : stations )
    {
        Station sta;
        sta.id = row.at("id");
        sta.latitude = std::stod(row.at("latitude"));
        sta.longitude = std::stod(row.at("longitude"));
        sta.elevation = std::stod(row.at("elevation"));
        sta.networkCode = row.at("networkCode");
        sta.stationCode = row.at("stationCode");
        _stations[sta.id] = sta;
    }

    vector<map<string,string> >events = CSV::readWithHeader(eventFile);

    for (const auto& row : events )
    {
        Event ev;
        ev.id          = std::stoul(row.at("id"));
        ev.time        = Core::Time::FromString(row.at("isotime").c_str(), "%FT%T.%fZ"); //iso format
        ev.latitude    = std::stod(row.at("latitude"));
        ev.longitude   = std::stod(row.at("longitude"));
        ev.depth       = std::stod(row.at("depth"));
        ev.magnitude   = std::stod(row.at("magnitude"));
        ev.horiz_err   = std::stod(row.at("horizontal_err"));
        ev.vert_err    = std::stod(row.at("vertical_err"));
        ev.rms         = std::stod(row.at("rms"));
        ev.relocInfo.isRelocated = false;
        if ( loadRelocationInfo && (row.count("relocated") != 0) && strToBool(row.at("relocated")) )
        {
            ev.relocInfo.isRelocated = true;
            ev.relocInfo.lonUncertainty   = std::stod(row.at("lonUncertainty"));
            ev.relocInfo.latUncertainty   = std::stod(row.at("latUncertainty"));
            ev.relocInfo.depthUncertainty = std::stod(row.at("depthUncertainty"));
            ev.relocInfo.numCCp           = std::stoi(row.at("numCCp"));
            ev.relocInfo.numCCs           = std::stoi(row.at("numCCs"));
            ev.relocInfo.numCTp           = std::stoi(row.at("numCTp"));
            ev.relocInfo.numCTs           = std::stoi(row.at("numCTs"));
            ev.relocInfo.rmsResidualCC    = std::stod(row.at("residualCC"));
            ev.relocInfo.rmsResidualCT    = std::stod(row.at("residualCT"));
            ev.relocInfo.numNeighbours    = std::stod(row.at("numNeighbours"));
        }
        _events[ev.id] = ev;
    }

    vector<map<string,string> >phases = CSV::readWithHeader(phaFile);

    for (const auto& row : phases )
    {
        Phase ph;
        ph.eventId          = std::stoul(row.at("eventId"));
        ph.stationId        = row.at("stationId");
        ph.time             = Core::Time::FromString(row.at("isotime").c_str(), "%FT%T.%fZ"); //iso format
        ph.lowerUncertainty = std::stod(row.at("lowerUncertainty"));
        ph.upperUncertainty = std::stod(row.at("upperUncertainty"));
        ph.type             = row.at("type");
        ph.networkCode      = row.at("networkCode");
        ph.stationCode      = row.at("stationCode");
        ph.locationCode     = row.at("locationCode");
        ph.channelCode      = row.at("channelCode");
        ph.isManual         = row.at("evalMode") == "manual";
        ph.relocInfo.isRelocated = false;
        if ( loadRelocationInfo && (row.count("usedInReloc") != 0) && strToBool(row.at("usedInReloc")) )
        {
            ph.relocInfo.isRelocated  = true;
            ph.procInfo.weight        = std::stod(row.at("initialWeight"));
            ph.relocInfo.finalWeight  = std::stod(row.at("finalWeight"));
            ph.relocInfo.residual     = std::stod(row.at("residual"));
        }
        _phases.emplace(ph.eventId, ph);
    }
}



void Catalog::add(const std::vector<DataModel::OriginPtr>& origins,
                  DataSource& dataSrc)
{
    for(DataModel::OriginPtr org : origins)
    {
        if ( org->arrivalCount() == 0)
            dataSrc.loadArrivals(org.get()); // try to load arrivals

        if ( org->arrivalCount() == 0)
        {
            SEISCOMP_WARNING("Origin %s doesn't have any arrival. Skip it.", org->publicID().c_str());
            continue;
        }


        // Add event
        Event ev;
        ev.id          = 0;
        ev.time        = org->time().value();
        ev.latitude    = org->latitude();
        ev.longitude   = org->longitude();
        ev.depth       = org->depth(); // km
        ev.horiz_err   = 0;
        ev.vert_err    = 0;
        try {
            ev.rms = org->quality().standardError();
        } catch ( ... ) {  ev.rms = 0; }

        DataModel::MagnitudePtr mag;
        // try to fetch preferred magnitude stored in the event
        DataModel::EventPtr parentEvent = dataSrc.getParentEvent(org->publicID());
        if ( parentEvent )
        {
            mag = dataSrc.get<DataModel::Magnitude>(parentEvent->preferredMagnitudeID());
        }
        if ( mag )
        {
            ev.magnitude = mag->magnitude();
        }
        else
        {
            SEISCOMP_WARNING("Origin %s: cannot load preferred magnitude from parent event, set it to 0",
                             org->publicID().c_str());
            ev.magnitude = 0.;
        }

        SEISCOMP_DEBUG("Adding origin '%s' to the catalog", org->publicID().c_str());

        // make sure we didn't ger duplicated events
        if ( this->searchEvent(ev) != _events.end())
        {
            SEISCOMP_WARNING("Skipping duplicated origin %s (%s)",
                             org->publicID().c_str(), string(ev).c_str());
            continue;
        }

        this->addEvent(ev, false, false);
        Event& newEvent = searchByValue(this->_events, ev)->second;  // fetch the id

        // Add Phases
        for ( size_t i = 0; i < org->arrivalCount(); ++i )
        {
            DataModel::Arrival *orgArr = org->arrival(i);
            const DataModel::Phase& orgPh = orgArr->phase();

            DataModel::PickPtr pick = dataSrc.get<DataModel::Pick>(orgArr->pickID());
            if ( !pick )
            {
                SEISCOMP_ERROR("Cannot load pick '%s' (origin %s)",
                               orgArr->pickID().c_str(), org->publicID().c_str());
                continue;
            }

            // find the station
            Station sta;
            sta.networkCode = pick->waveformID().networkCode();
            sta.stationCode = pick->waveformID().stationCode();
 
            // skip not selected picks/phases for manaul origins only
            if ( org->evaluationMode() == Seiscomp::DataModel::MANUAL &&
                 ! orgArr->timeUsed() )
            {
                SEISCOMP_DEBUG("Skip not used %s phase %s.%s in manual origin",
                               orgPh.code().c_str(), sta.networkCode.c_str(), sta.stationCode.c_str()  );
                continue;
            }

            // add station if not already there
            if (searchStation(sta) == _stations.end())
            {
                DataModel::Station* orgArrStation = findStation(sta.networkCode, sta.stationCode,
                                                                pick->time());
                if ( !orgArrStation )
                {
                    SEISCOMP_ERROR("Cannot load station %s.%s information for arrival '%s' (origin '%s')."
                                   "All picks associated with this station will not be used.",
                                    sta.networkCode.c_str(), sta.stationCode.c_str(),
                                    orgArr->pickID().c_str(), org->publicID().c_str());
                    continue;
                }
                sta.latitude = orgArrStation->latitude();
                sta.longitude = orgArrStation->longitude();
                sta.elevation = orgArrStation->elevation(); // meter
                this->addStation(sta, false);
            }
            // the station has to be there at this point
            sta = this->searchStation(sta)->second;

            // get uncertainty
            pair<double,double> uncertainty = getPickUncertainty(pick.get());

            Phase ph;
            ph.eventId          = newEvent.id;
            ph.stationId        = sta.id;
            ph.time             = pick->time().value();
            ph.lowerUncertainty = uncertainty.first;
            ph.upperUncertainty = uncertainty.second;
            ph.type             = orgPh.code();
            ph.networkCode      = pick->waveformID().networkCode();
            ph.stationCode      = pick->waveformID().stationCode();
            ph.locationCode     = pick->waveformID().locationCode();
            ph.channelCode      = pick->waveformID().channelCode();
            ph.isManual         = (pick->evaluationMode() == Seiscomp::DataModel::MANUAL);
            this->addPhase(ph, true, false);
        }
    }
}



void Catalog::add(const std::vector<std::string>& ids, DataSource& dataSrc)
{
    vector<DataModel::OriginPtr> origins;

    for(const string& id : ids)
    {
        DataModel::OriginPtr org = dataSrc.get<DataModel::Origin>(id);
        if ( !org )
        {
            SEISCOMP_ERROR("Cannot find origin with id %s", id.c_str());
            continue;
        }
        origins.push_back(org);
    }

    add(origins, dataSrc);
}



void Catalog::add(const std::string& idFile, DataSource& dataSrc)
{
    if ( !Util::fileExists(idFile) )
    {
        string msg = "File " + idFile + " does not exist";
        throw runtime_error(msg);
    }

    vector<string> ids;
    vector< map<string,string> > rows = CSV::readWithHeader(idFile);

    for(const auto& row : rows)
    {
        const string& id = row.at("seiscompId");
        ids.push_back(id);
    }

    add(ids, dataSrc);
}



void Catalog::add(const Catalog& other, bool keepEvId)
{
    for (const auto& kv :  other.getEvents() )
    {
        const Catalog::Event& event = kv.second;
        if ( keepEvId && _events.find(event.id) != _events.end() )
        {
            SEISCOMP_DEBUG("Skipping duplicated event id %u", event.id);
            continue;
        }
        this->add(event.id, other, keepEvId);
    }
}



CatalogPtr Catalog::extractEvent(unsigned eventId, bool keepEvId) const
{
    CatalogPtr eventToExtract = new Catalog();

    auto search = this->getEvents().find(eventId);
    if (search == this->getEvents().end())
    {
        string msg = stringify("Cannot find event id %u in the catalog.", eventId);
        throw runtime_error(msg);
    }

    const Catalog::Event& event = search->second;
    unsigned newEventId;

    if ( keepEvId )
    {
        eventToExtract->_events[event.id] = event;
        newEventId = event.id;
    }
    else
    {
        eventToExtract->addEvent(event, false, false);
        newEventId = eventToExtract->searchEvent(event)->first;
    }

    auto eqlrng = this->getPhases().equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
        Catalog::Phase phase = it->second;

        auto search = this->getStations().find(phase.stationId);
        if (search == this->getStations().end())
        {
            string msg = stringify("Malformed catalog: cannot find station '%s' "
                                   " referenced by phase '%s' for event '%s'",
                                   phase.stationId.c_str(), string(phase).c_str(), string(event).c_str());
            throw runtime_error(msg);
        }
        const Catalog::Station& station = search->second;
        eventToExtract->addStation(station, true);

        phase.eventId = newEventId;
        eventToExtract->addPhase(phase, false, false);
    } 

    return eventToExtract;
}



unsigned Catalog::add(unsigned evId, const Catalog& evCat, bool keepEvId)
{
    unsigned newEventId;

    const Catalog::Event& event = evCat._events.find(evId)->second;

    if ( keepEvId )
    {
        if  (_events.find(event.id) != _events.end() )
            throw runtime_error("Cannot add event, internal logic error");
        _events[event.id] = event;
        newEventId = event.id;
    }
    else
    {
        // don't add an event with same values, but keep merging phases
        addEvent(event, false, false);
        newEventId = searchEvent(event)->first;
    }

    auto eqlrng = evCat._phases.equal_range(event.id);
    for (auto it = eqlrng.first; it != eqlrng.second; ++it)
    {
        Catalog::Phase phase = it->second;

        const Catalog::Station& station = evCat._stations.at(phase.stationId); 
        addStation(station, true);

        phase.eventId = newEventId;
        addPhase(phase, false, false);
    }

    return newEventId;
}



void Catalog::removeEvent(const Event& event)
{
    map<unsigned,Catalog::Event>::const_iterator it = searchEvent(event);
    if ( it != _events.end() )
    {
        _events.erase(it);
        auto eqlrng = _phases.equal_range(it->second.id);
        _phases.erase(eqlrng.first, eqlrng.second);
    }
}


void Catalog::removeEvent(unsigned eventId)
{
    map<unsigned,Catalog::Event>::const_iterator it = _events.find(eventId);
    if ( it != _events.end() )
    {
        _events.erase(it);
    }
    auto eqlrng = _phases.equal_range(eventId);
    _phases.erase(eqlrng.first, eqlrng.second);
}


void Catalog::removePhase(const Phase& phase)
{
    std::map<unsigned,Phase>::const_iterator it = searchPhase(phase);
    if ( it != _phases.end() )
    {
        _phases.erase(it);
    }
}


void Catalog::removePhase(unsigned eventId, const std::string& stationId, const string& type)
{
    auto eqlrng = _phases.equal_range(eventId);
    auto it = eqlrng.first;
    while ( it != eqlrng.second )
    {
        const Catalog::Phase& ph = it->second;
        if ( ph.stationId == stationId && ph.type == type )
            _phases.erase(it++);
        else
            ++it;
    }
}


bool Catalog::updateStation(const Station& newStation)
{
    map<string,Catalog::Station>::iterator it = _stations.find(newStation.id);
    if ( it != _stations.end() )
    {
        it->second = newStation;
        return true;
    }
    return false;
}


bool Catalog::updateEvent(const Event& newEv)
{
    map<unsigned,Catalog::Event>::iterator it = _events.find(newEv.id);
    if ( it != _events.end() )
    {
        it->second = newEv;
        return true;
    }
    return false;
}


bool Catalog::updatePhase(const Phase& newPh)
{
    auto eqlrng = _phases.equal_range(newPh.eventId);
    auto it = eqlrng.first;
    bool updated = false;
    while ( it != eqlrng.second )
    {
        Catalog::Phase& oldPh = it->second;
        if ( oldPh.stationId == newPh.stationId )
        {
            oldPh = newPh;
            updated = true;
        }
    }
    return updated;
}


map<string,Catalog::Station>::const_iterator Catalog::searchStation(const Station& station) const
{
    return searchByValue(_stations, station);
}


map<unsigned,Catalog::Event>::const_iterator Catalog::searchEvent(const Event& event) const
{
    return searchByValue(_events, event);
}


map<unsigned,Catalog::Phase>::const_iterator Catalog::searchPhase(const Phase& phase) const
{
    return searchByValue(_phases, phase);
}


bool Catalog::addStation(const Station& station, bool checkDuplicate)
{
    if ( checkDuplicate && _stations.find(station.id) != _stations.end() )
    {
        return false;
    }

    Station newStation = station;
    newStation.id = newStation.networkCode + newStation.stationCode;
    _stations[newStation.id] = newStation;
    return true;
}


bool Catalog::addEvent(const Event& event, bool checkDuplicateValue, bool checkDuplicateId)
{
    if ( checkDuplicateValue && searchEvent(event) != _events.end() )
    {
        return false;
    }

    if ( checkDuplicateId && _events.find(event.id) != _events.end() )
    {
        return false;
    }

    decltype(_events)::key_type maxKey = _events.empty() ? 0 : _events.rbegin()->first;
    Event newEvent = event;
    newEvent.id = maxKey + 1;
    _events[newEvent.id] = newEvent;
    return true;
}


bool Catalog::addPhase(const Phase& phase, bool checkDuplicateValue, bool checkDuplicateId)
{
    if ( checkDuplicateValue && searchPhase(phase) != _phases.end() )
    {
        return false;
    }

    if ( checkDuplicateId )
    {
        auto eqlrng = _phases.equal_range(phase.eventId);
        for (auto it = eqlrng.first; it != eqlrng.second; ++it)
        {
            const Catalog::Phase& currPh = it->second;
            if ( phase.stationId == currPh.stationId )
                return false;
        }
    }

    _phases.emplace(phase.eventId, phase);
    return true;
}

void Catalog::writeToFile(string eventFile, string phaseFile, string stationFile) const
{
    /*
     * Write Events
     * */
    stringstream evStreamNoReloc;
    stringstream evStreamReloc;

    evStreamNoReloc << "id,isotime,latitude,longitude,depth,magnitude,horizontal_err,vertical_err,rms"; 
    evStreamReloc << evStreamNoReloc.str() << ",relocated,numNeighbours,numCCp,numCCs,numCTp,numCTs,residualCC,residualCT,lonUncertainty,latUncertainty,depthUncertainty" << endl;
    evStreamNoReloc << endl;

    bool relocInfo = false;
    for (const auto& kv : _events )
    {
        const Catalog::Event& ev = kv.second;

        stringstream evStream;
        evStream << stringify("%u,%s,%.6f,%.6f,%.4f,%.2f,%.4f,%.4f,%.4f",
                              ev.id,ev.time.iso().c_str(),
                              ev.latitude,ev.longitude,ev.depth,ev.magnitude,
                              ev.horiz_err,ev.vert_err,ev.rms);

        evStreamNoReloc << evStream.str() << endl;
        evStreamReloc   << evStream.str();

        if ( ! ev.relocInfo.isRelocated )
        {
            evStreamReloc << ",false,,,,,,,,,,";
        }
        else
        {
            relocInfo = true;
            evStreamReloc <<  stringify(",true,%d,%d,%d,%d,%d,%.4f,%.4f,%.6f,%.6f,%.6f",
                                  ev.relocInfo.numNeighbours, ev.relocInfo.numCCp, ev.relocInfo.numCCs,
                                  ev.relocInfo.numCTp, ev.relocInfo.numCTs,
                                  ev.relocInfo.rmsResidualCC, ev.relocInfo.rmsResidualCT, 
                                  ev.relocInfo.lonUncertainty, ev.relocInfo.latUncertainty,
                                  ev.relocInfo.depthUncertainty );
        }
        evStreamReloc << endl;
    }

    ofstream evStream(eventFile);
    evStream << ( relocInfo ? evStreamReloc.str() : evStreamNoReloc.str() );

    /*
     * Write Phases
     * */
    ofstream phStream(phaseFile);

    phStream << "eventId,stationId,isotime,lowerUncertainty,upperUncertainty,type,networkCode,stationCode,locationCode,channelCode,evalMode";
    if (relocInfo)
    {
        phStream << ",usedInReloc,initialWeight,finalWeight,residual";
    }
    phStream << endl;

    for (const auto& kv : _phases )
    {
        const Catalog::Phase& ph = kv.second;
        phStream << stringify("%u,%s,%s,%.3f,%.3f,%s,%s,%s,%s,%s,%s",
                              ph.eventId, ph.stationId.c_str(), ph.time.iso().c_str(),
                              ph.lowerUncertainty, ph.upperUncertainty, ph.type.c_str(),
                              ph.networkCode.c_str(), ph.stationCode.c_str(),
                              ph.locationCode.c_str(), ph.channelCode.c_str(),
                              (ph.isManual ? "manual" : "automatic"));

        if (relocInfo)
        {
            if ( ! ph.relocInfo.isRelocated )
            {
                phStream << ",false,,,";
            }
            else
            {
                phStream << stringify(",true,%.4f,%.3f,%.3f", ph.procInfo.weight,
                                      ph.relocInfo.finalWeight, ph.relocInfo.residual);
            }
        }
        phStream << endl;
    }

    /*
     * Write Stations
     * */
    ofstream staStream(stationFile);
    staStream << "id,latitude,longitude,elevation,networkCode,stationCode" << endl;
    for (const auto& kv : _stations )
    {
        const Catalog::Station& sta = kv.second;
        staStream << stringify("%s,%.6f,%.6f,%.1f,%s,%s",
                               sta.id.c_str(), sta.latitude, sta.longitude, sta.elevation,
                               sta.networkCode.c_str(), sta.stationCode.c_str())
                  << endl;
    }
}


} // HDD
} // Seiscomp

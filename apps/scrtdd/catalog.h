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

#ifndef __RTDD_APPLICATIONS_CATALOG_H__
#define __RTDD_APPLICATIONS_CATALOG_H__

#include <seiscomp3/core/baseobject.h>
#include <seiscomp3/datamodel/eventparameters.h>
#include <seiscomp3/datamodel/publicobjectcache.h>
#include <seiscomp3/datamodel/databasequery.h>
#include <seiscomp3/datamodel/origin.h>

#include <map>
#include <vector>

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



DEFINE_SMARTPOINTER(Catalog);

// DD background catalog
class Catalog : public Core::BaseObject {

    public:

        struct Station {
            std::string id;
            double latitude;
            double longitude;
            double elevation; // meter
            std::string networkCode;
            std::string stationCode;

            // this equality works between multiple catalogs (same id is not required)
            bool operator==(const Station& other) const
            {
             return (networkCode == other.networkCode) &&
                    (stationCode == other.stationCode);
            }
            bool operator!=(const Station& other) const
            {
                return !operator==(other);
            }
            operator std::string() const
            {
                return id;
            }
        };

        struct Event {
            unsigned id; // makes it unique in the catalog 
            Core::Time time;
            double latitude;
            double longitude;
            double depth;   // km
            double magnitude;
            double horiz_err;
            double vert_err;
            double rms;
            struct {
                bool isRelocated = false;
                double lonUncertainty;
                double latUncertainty;
                double depthUncertainty;
                int numCCp;
                int numCCs;
                int numCTp;
                int numCTs;
                double rmsResidualCC;
                double rmsResidualCT;
            } relocInfo;

            // this equality works between multiple catalogs (same id is not required)
            bool operator==(const Event& other) const
            {
             return (time == other.time) &&
                    (latitude == other.latitude) &&
                    (longitude == other.longitude) &&
                    (depth == other.depth) &&
                    (magnitude == other.magnitude);
            }
            bool operator!=(const Event& other) const
            {
                return !operator==(other);
            }
            operator std::string() const
            {
                return std::to_string(id);
            }
        };

        struct Phase {
            unsigned eventId;
            std::string stationId;
            Core::Time time;
            double lowerUncertainty;
            double upperUncertainty;
            std::string type;
            std::string networkCode;
            std::string stationCode;
            std::string locationCode;
            std::string channelCode;
            bool isManual;
            struct {
                bool isRelocated = false;
                double weight;       // 0-1 interval
                double finalWeight;
                double residual;
            } relocInfo;

            // this equality works between multiple catalogs (same id is not required) 
            bool operator==(const Phase& other) const
            {
                return (eventId == other.eventId) &&
                       (stationId == other.stationId) &&
                       (time == other.time) &&
                       (type == other.type) &&
                       (networkCode == other.networkCode) &&
                       (stationCode == other.stationCode) &&
                       (locationCode == other.locationCode) &&
                       (channelCode == other.channelCode) &&
                       (isManual == other.isManual);
            }
            bool operator!=(const Phase& other) const
            {
                return !operator==(other);
            }
            operator std::string() const
            {
                return networkCode + "." +  stationCode + "." +
                       locationCode + "." + channelCode + " type " +
                       type + (isManual ? " (manual) " : " (auto) " ) + time.iso() +
                       " evId " + std::to_string(eventId)  + " staId " + stationId;
            }
        };

        Catalog();
        virtual ~Catalog() { }

        Catalog(const Catalog& other);

        // custom data format constructors
        Catalog(const std::map<std::string,Station>& stations,
                const std::map<unsigned,Event>& events,
                const std::multimap<unsigned,Phase>& phases);
        Catalog(const std::string& stationFile,
                const std::string& catalogFile,
                const std::string& phaFile);

        // populate from seiscomp data format
        void add(const std::vector<DataModel::OriginPtr>& origins, DataSource& dataSrc);
        void add(const std::vector<std::string>& ids, DataSource& dataSrc);
        void add(const std::string& idFile, DataSource& dataSrc);

        CatalogPtr merge(const CatalogCPtr& other, bool skipExistingId) const;
        CatalogPtr extractEvent(unsigned eventId) const;
        bool copyEvent(const Catalog::Event& event, const CatalogCPtr& eventCatalog, bool keepEvId);

        void removeEvent(const Event& event);
        void removeEvent(unsigned eventId);
        void removePhase(const Phase& phase);
        void removePhase(unsigned eventId, const std::string& stationId, const std::string& type);

        bool addStation(const Station&, bool checkDuplicateValue, bool checkDuplicateId);
        bool addEvent(const Event&, bool checkDuplicateValue, bool checkDuplicateId);
        bool addPhase(const Phase&, bool checkDuplicateValue, bool checkDuplicateId);

        bool updateStation(const Station& newStation);
        bool updateEvent(const Event& newEv);
        bool updatePhase(const Phase& newPh);

        const std::map<std::string,Station>& getStations() const { return _stations;}
        const std::map<unsigned,Event>& getEvents() const { return _events;}
        const std::multimap<unsigned,Phase>& getPhases() const { return _phases;}

        std::map<std::string,Station>::const_iterator searchStation(const Station&) const;
        std::map<unsigned,Event>::const_iterator searchEvent(const Event&) const;
        std::map<unsigned,Phase>::const_iterator searchPhase(const Phase&) const;

        void writeToFile(std::string eventFile,
                         std::string phaseFile,
                         std::string stationFile) const;

    /*
     *  static
     */

        static DataModel::Station* findStation(const std::string& netCode,
                                               const std::string& stationCode,
                                               const Core::Time& atTime);

    private:
        Catalog& operator=(const Catalog& other);

        std::map<std::string,Station> _stations; // indexed by station id
        std::map<unsigned,Event> _events; //indexed by event id
        std::multimap<unsigned,Phase> _phases; //indexed by event id
};


}
}

#endif

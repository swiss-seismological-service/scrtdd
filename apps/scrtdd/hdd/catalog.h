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

#ifndef __RTDD_APPLICATIONS_CATALOG_H__
#define __RTDD_APPLICATIONS_CATALOG_H__

#include "datasrc.h"

#include <seiscomp3/core/baseobject.h>
#include <seiscomp3/datamodel/eventparameters.h>
#include <seiscomp3/datamodel/publicobjectcache.h>
#include <seiscomp3/datamodel/databasequery.h>
#include <seiscomp3/datamodel/origin.h>

#include <unordered_map>
#include <vector>

namespace Seiscomp {
namespace HDD {

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
            std::string locationCode;

            // search by value when the Id is not known (works between multiple catalogs )
            bool operator==(const Station& other) const
            {
                return (networkCode == other.networkCode) &&
                        (stationCode == other.stationCode) &&
                        (locationCode == other.locationCode) &&
                        (latitude == other.latitude) &&
                        (longitude == other.longitude) &&
                        (elevation == other.elevation);
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
            double rms;

            struct {
                bool isRelocated = false;
                unsigned numNeighbours;
                unsigned usedP;
                unsigned usedS;
                unsigned numCCp;
                unsigned numCCs;
                unsigned numCTp;
                unsigned numCTs;
                double meanObsWeight;
                double meanFinalObsWeight;
            } relocInfo;

            // search by value when the Id is not known (works between multiple catalogs )
            bool operator==(const Event& other) const
            {
             return (time == other.time) &&
                    (latitude == other.latitude) &&
                    (longitude == other.longitude) &&
                    (depth == other.depth) &&
                    (magnitude == other.magnitude) &&
                    (rms == other.rms);
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

            enum class Type : char { P='P', S='S' };

            enum class Source { CATALOG, RT_EVENT, THEORETICAL, XCORR };

            struct {
                Type type;
                double weight;       // 0-1 interval
                Source source;
            } procInfo;

            struct {
                bool isRelocated = false;
                double residual;
                double finalWeight;
                unsigned numObservs;
                unsigned numXcorrObservs;
                double meanObsWeight;
                double meanFinalObsWeight;
            } relocInfo;


            // search by value when the Id is not known (works between multiple catalogs )
            bool operator==(const Phase& other) const
            {
                return (time == other.time) &&
                       (lowerUncertainty == other.lowerUncertainty) &&
                       (upperUncertainty == other.upperUncertainty) &&
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
                return "\"Phase " + type + (isManual ? " (manual) " : " (auto) " ) +
                        networkCode + "." +  stationCode + "." + locationCode + "." + channelCode +
                        " " + time.iso() + " evId " + std::to_string(eventId)  + " staId " + stationId + "\"";
            }
        };

        Catalog();
        virtual ~Catalog() { }

        // copy constructor / assignment operator
        Catalog(const Catalog& other);
        Catalog& operator=(const Catalog& other);

        // move constructor / assignment operator
        Catalog(Catalog&& other);
        Catalog& operator=(Catalog&& other);

        // custom data format constructors
        Catalog(std::unordered_map<std::string,Station>&& stations,
                std::map<unsigned,Event>&& events,
                std::unordered_multimap<unsigned,Phase>&& phases);
        Catalog(const std::unordered_map<std::string,Station>& stations,
                const std::map<unsigned,Event>& events,
                const std::unordered_multimap<unsigned,Phase>& phases);
        Catalog(const std::string& stationFile,
                const std::string& eventFile,
                const std::string& phaseFile,
                bool loadRelocationInfo=false);

        // populate from seiscomp data format
        void add(const std::vector<DataModel::OriginPtr>& origins, DataSource& dataSrc);
        void add(const std::vector<std::string>& ids, DataSource& dataSrc);
        void add(const std::string& idFile, DataSource& dataSrc);

        void add(const Catalog& other, bool keepEvId);
        unsigned add(unsigned evId, const Catalog& eventCatalog, bool keepEvId);
        CatalogPtr extractEvent(unsigned eventId, bool keepEvId) const;

        void removeEvent(unsigned eventId);
        void removePhase(unsigned eventId, const std::string& stationId, const Phase::Type& type);

        std::string addStation(const Station&);
        unsigned addEvent(const Event&);
        void addPhase(const Phase&);

        bool updateStation(const Station& newStation, bool addIfMissing=false);
        bool updateEvent(const Event& newEv, bool addIfMissing=false);
        bool updatePhase(const Phase& newPh, bool addIfMissing=false);

        const std::unordered_map<std::string,Station>& getStations() const { return _stations;}
        const std::map<unsigned,Event>& getEvents() const { return _events;}
        const std::unordered_multimap<unsigned,Phase>& getPhases() const { return _phases;}

        std::unordered_map<std::string,Station>::const_iterator
        searchStation(const std::string& networkCode,
                      const std::string& stationCode,
                      const std::string& locationCode) const;
        std::map<unsigned,Event>::const_iterator searchEvent(const Event&) const;
        std::unordered_map<unsigned,Phase>::const_iterator
        searchPhase(unsigned eventId, const std::string& stationId, const Phase::Type& type) const;

        void writeToFile(std::string eventFile,
                         std::string phaseFile,
                         std::string stationFile) const;

        //
        //  static
        //
        static double computePickWeight(double uncertainty);
        static double computePickWeight(const Catalog::Phase& phase);
        static CatalogPtr filterPhasesAndSetWeights(const CatalogCPtr& catalog,
                                             const Catalog::Phase::Source& source,
                                             const std::vector<std::string>& PphaseToKeep,
                                             const std::vector<std::string>& SphaseToKeep);

        static DataModel::Station* findStation(const std::string& netCode,
                                               const std::string& stationCode,
                                               const Core::Time& atTime);

        static DataModel::SensorLocation* findSensorLocation(const std::string &networkCode,
                                                      const std::string &stationCode,
                                                      const std::string &locationCode,
                                                      const Core::Time &atTime);

        static constexpr double DEFAULT_MANUAL_PICK_UNCERTAINTY    = 0.030;
        static constexpr double DEFAULT_AUTOMATIC_PICK_UNCERTAINTY = 0.100;

    private:

        std::unordered_map<std::string,Station> _stations; // indexed by station id
        std::map<unsigned,Event> _events; //indexed by event id
        std::unordered_multimap<unsigned,Phase> _phases; //indexed by event id
};


}
}

#endif

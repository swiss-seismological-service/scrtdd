/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 *   This program is free software: you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation, either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#ifndef __SEISCOMP_APPLICATIONS_HYPODD_H__
#define __SEISCOMP_APPLICATIONS_HYPODD_H__

#include <seiscomp3/core/baseobject.h>
#include <seiscomp3/core/genericrecord.h>
#include <seiscomp3/core/recordsequence.h>
#include <seiscomp3/datamodel/eventparameters.h>
#include <seiscomp3/datamodel/publicobjectcache.h>
#include <seiscomp3/datamodel/databasequery.h>
#include <seiscomp3/datamodel/origin.h>

#include <map>
#include <vector>

#include "hypodd.h"

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
            std::string type;
            double weight;       // 0-1 interval
            std::string networkCode;
            std::string stationCode;
            std::string locationCode;
            std::string channelCode;
            struct {
                bool isRelocated = false;
                double residual;
                double finalWeight;
            } relocInfo;

            // this equality works between multiple catalogs (same id is not required) 
            bool operator==(const Phase& other) const
            {
                return (time == other.time) &&
                       (type == other.type) &&
                       (networkCode == other.networkCode) &&
                       (stationCode == other.stationCode) &&
                       (locationCode == other.locationCode) &&
                       (channelCode == other.channelCode);
            }
            bool operator!=(const Phase& other) const
            {
                return !operator==(other);
            }
            operator std::string() const
            {
                return type + " " + time.iso() + " " + networkCode + "." +
                       stationCode + "." + locationCode + "." + channelCode + 
                       " evId " + std::to_string(eventId)  + " staId " + stationId;
            }
        };

        Catalog();
        virtual ~Catalog() { }

        // custom data format constructors
        Catalog(const std::map<std::string,Station>& stations,
                const std::map<unsigned,Event>& events,
                const std::multimap<unsigned,Phase>& phases);
        Catalog(const std::string& stationFile,
                const std::string& catalogFile,
                const std::string& phaFile);

        // add from seiscomp data format
        void add(const std::vector<DataModel::Origin*>& origins, DataSource& dataSrc);
        void add(const std::vector<std::string>& ids, DataSource& dataSrc);
        void add(const std::string& idFile, DataSource& dataSrc);

        CatalogPtr merge(const CatalogCPtr& other) const;
        bool copyEvent(const Catalog::Event& event, const CatalogCPtr& other, bool keepEvId);
        void removeEvent(const Event& event);
        void removeEvent(unsigned eventId);

        bool addStation(const Station&, bool checkDuplicate);
        bool addEvent(const Event&, bool checkDuplicate);
        bool addPhase(const Phase&, bool checkDuplicate);

        const std::map<std::string,Station>& getStations() const { return _stations;}
        const std::map<unsigned,Event>& getEvents() const { return _events;}
        const std::multimap<unsigned,Phase>& getPhases() const { return _phases;}

        std::map<std::string,Station>::const_iterator searchStation(const Station&) const;
        std::map<unsigned,Event>::const_iterator searchEvent(const Event&) const;
        std::map<unsigned,Phase>::const_iterator searchPhase(const Phase&) const;

        void writeToFile(std::string eventFile,
                         std::string phaseFile,
                         std::string stationFile) const;

    private:
        std::map<std::string,Station> _stations; // indexed by station id
        std::map<unsigned,Event> _events; //indexed by event id
        std::multimap<unsigned,Phase> _phases; //indexed by event id
};



struct Config {

    std::vector<std::string> validPphases = {"Pg,P"};
    std::vector<std::string> validSphases = {"Sg,S"};

    // hypodd executable specific
    struct {
        std::string exec = "hypodd";
        std::string step1CtrlFile;
        std::string step2CtrlFile;
    } hypodd;

    // ph2dt executable specific
    struct {
        std::string exec = "ph2dt";
        std::string ctrlFile;
    } ph2dt;

    // differential travel time specific
    struct {
        double minWeight      = 0;  // Min weight of phases required (0-1)
        double minEStoIEratio = 0;  // Min epi-sta to  interevent distance ration required
        double minESdist      = 0;  // Min epi-sta distance required
        double maxESdist      =-1;  // Max epi-sta distance allowed
        double maxIEdist      =-1;  // Max interevent-distance allowed (km)
        int minNumNeigh       = 1;  // Min neighbors required
        int maxNumNeigh       =-1;  // Max neighbors allowed (furthest events are discarded)
        int minDTperEvt       = 1;  // Min differential times per event pair required (Including P+S)
    } dtct;

    // cross correlation specific
    struct {
        std::string recordStreamURL;

        double minWeight      = 0;  // Min weight of phases required (0-1)
        double minEStoIEratio = 0;  // Min epi-sta to  interevent distance ration required
        double minESdist      = 0;  // Min epi-sta distance required
        double maxESdist      =-1;  // Max epi-sta distance allowed
        double maxIEdist      =-1;  // Max interevent-distance allowed (km)
        int minNumNeigh       = 1;  // Min neighbors required
        int maxNumNeigh       =-1;  // Max neighbors allowed (furthest events are discarded)
        int minDTperEvt       = 1;  // Min differential times per event pair required (Including P+S)
    } dtcc;

    struct {
        double minCoef        = 0;    // Min xcorr coefficient required (0-1)

        double timeBeforePick; // secs
        double timeAfterPick;  // secs
        double maxDelay;       //secs
    } xcorr;

    struct {
        int filterOrder;
        double filterFmin     = 0;
        double filterFmax     = 0;

        double resampleFreq = 0;
    } wfFilter;

    struct {
        double minSnr = 0;
        double noiseStart = 0;
        double noiseEnd = 0;
        double signalStart = 0;
        double signalEnd = 0;
    } snr;
};



DEFINE_SMARTPOINTER(HypoDD);

class HypoDD : public Core::BaseObject {

    public:
        HypoDD(const CatalogCPtr& catalog, const Config& cfg, const std::string& workingDir);
        virtual ~HypoDD();

        void preloadData();

        CatalogCPtr getCatalog() { return _ddbgc; }
        void setCatalog(const CatalogCPtr& catalog);

        CatalogPtr relocateCatalog(bool force = true, bool usePh2dt = false);
        CatalogPtr relocateSingleEvent(const CatalogCPtr& orgToRelocate);

        void setWorkingDirCleanup(bool cleanup) { _workingDirCleanup = cleanup; }
        bool workingDirCleanup() { return _workingDirCleanup; }

        void setUseCatalogDiskCache(bool cache) { _useCatalogDiskCache = cache; }
        bool useCatalogDiskCache() { return _useCatalogDiskCache; }

        void setUseSingleEvDiskCache(bool cache) { _useSingleEvDiskCache = cache; }
        bool useSingleEvDiskCache() { return _useSingleEvDiskCache; }


    private:
        CatalogPtr filterOutPhases(const CatalogCPtr& catalog,
                                   const std::vector<std::string>& PphaseToKeep,
                                   const std::vector<std::string>& SphaseToKeep) const;
        void createStationDatFile(const std::string& staFileName, const CatalogCPtr& catalog) const;
        void createPhaseDatFile(const std::string& phaseFileName, const CatalogCPtr& catalog) const;
        void createEventDatFile(const std::string& eventFileName, const CatalogCPtr& catalog) const;
        void createDtCtCatalog(const CatalogCPtr& catalog,
                               const std::string& dtctFile) const;
        void createDtCtSingleEvent(const CatalogCPtr& catalog,
                                   unsigned evToRelocateId,
                                   const std::string& dtctFile) const;
        void buildAbsTTimePairs(const CatalogCPtr& catalog,
                                 unsigned evToRelocateId,
                                 std::ofstream& outStream) const;
        void createDtCcPh2dt(const std::string& dtctFile, const std::string& dtccFile);
        void createDtCcCatalog(const CatalogCPtr& catalog, const std::string& dtccFile);
        void createDtCcSingleEvent(const CatalogCPtr& catalog,
                                   unsigned evToRelocateId,
                                   const std::string& dtccFile);
        void buildXcorrDiffTTimePairs(const CatalogCPtr& catalog,
                                      unsigned evToRelocateId,
                                      std::ofstream& outStream,
                                      std::map<std::string,GenericRecordPtr>& catalogCache,
                                      bool useDiskCacheCatalog,
                                      std::map<std::string,GenericRecordPtr>& refEvCache,
                                      bool useDiskCacheRefEv);
        bool xcorr(const Catalog::Event& event1, const Catalog::Phase& phase1,
                   const Catalog::Event& event2, const Catalog::Phase& phase2,
                   double& dtccOut, double& weightOut,
                   std::map<std::string,GenericRecordPtr>& cache1,  bool useDiskCache1,
                   std::map<std::string,GenericRecordPtr>& cache2,  bool useDiskCache2);

        void runHypodd(const std::string& workingDir, const std::string& dtccFile,
                       const std::string& dtctFile, const std::string& eventFile,
                       const std::string& stationFile, const std::string& ctrlFile) const;
        void runPh2dt(const std::string& workingDir,
                      const std::string& stationFile,
                      const std::string& phaseFile) const;
        CatalogPtr loadRelocatedCatalog(const CatalogCPtr& originalCatalog,
                                        const std::string& ddrelocFile,
                                        const std::string& ddresidualFile="") const;
        CatalogPtr selectNeighbouringEvents(const CatalogCPtr& catalog, const Catalog::Event& refEv,
                                            double minPhaseWeight = 0, double minESdis=0,
                                            double maxESdis=-1, double minEStoIEratio=0,
                                            double maxIEdis=-1, int minDTperEvt=1,
                                            int minNumNeigh=1, int maxNumNeigh=-1) const;
        std::map<unsigned,CatalogPtr> 
        selectNeighbouringEventsCatalog(const CatalogCPtr& catalog, double minPhaseWeight,
                                        double minESdis, double maxESdis,
                                        double minEStoIEratio, double maxIEdis,
                                        int minDTperEvt, int minNumNeigh,
                                        int maxNumNeigh) const;
        bool xcorr(const GenericRecordCPtr& tr1, const GenericRecordCPtr& tr2, double maxDelay,
                   double& delayOut, double& coeffOut) const;
        double S2Nratio(const GenericRecordCPtr& tr, const Core::Time& guidingPickTime,
                        double noiseOffsetStart, double noiseOffsetEnd,
                        double signalOffsetStart, double signalOffsetEnd) const;
        CatalogPtr extractEvent(const CatalogCPtr& catalog, unsigned eventId) const;
        GenericRecordPtr getWaveform(const Core::TimeWindow& tw,
                                     const Catalog::Event& ev,
                                     const Catalog::Phase& ph,
                                     std::map<std::string,GenericRecordPtr>& memCache,
                                     bool useDiskCache) const;
        GenericRecordPtr loadWaveform(const Core::TimeWindow& tw,
                                      const std::string& networkCode,
                                      const std::string& stationCode,
                                      const std::string& locationCode,
                                      const std::string& channelCode,
                                      bool useDiskCache) const;
        GenericRecordPtr readWaveformFromRecordStream(const Core::TimeWindow& tw,
                                                      const std::string& networkCode,
                                                      const std::string& stationCode,
                                                      const std::string& locationCode,
                                                      const std::string& channelCode) const;
        bool merge(GenericRecord &trace, const RecordSequence& seq) const;
        bool trim(GenericRecord &trace, const Core::TimeWindow& tw) const;
        void filter(GenericRecord &trace, bool demeaning=true,
                    int order=3, double fmin=-1, double fmax=-1, double resampleFreq=0) const;
        void resample(GenericRecord& trace, double sf, bool average) const;
        std::string generateWorkingSubDir(const Catalog::Event& ev) const;
        std::string waveformId(const Catalog::Phase& ph, const Core::TimeWindow& tw) const;
        std::string waveformId(const std::string& networkCode, const std::string& stationCode,
                               const std::string& locationCode, const std::string& channelCode,
                               const Core::TimeWindow& tw) const;

    private:
        std::string _workingDir;
        std::string _cacheDir;
        CatalogPtr _ddbgc;
        Config _cfg;
        bool _workingDirCleanup = true;
        bool _useCatalogDiskCache = false;
        bool _useSingleEvDiskCache = false;
        std::map<std::string, GenericRecordPtr> _wfCache;
        std::map<std::string, bool> _validWfs;
};

}
}

#endif

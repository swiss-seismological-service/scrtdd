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

#ifndef __RTDD_APPLICATIONS_HYPODD_H__
#define __RTDD_APPLICATIONS_HYPODD_H__

#include <seiscomp3/core/baseobject.h>
#include <seiscomp3/core/genericrecord.h>
#include <seiscomp3/core/recordsequence.h>
#include <seiscomp3/datamodel/origin.h>
#include <seiscomp3/datamodel/utils.h>

#include <set>
#include <map>
#include <vector>

#include "catalog.h"
#include "wfmngr.h"

namespace Seiscomp {
namespace HDD {

struct Config {

    std::vector<std::string> validPphases = {"Pg","P","Px"};
    std::vector<std::string> validSphases = {"Sg","S","Sx"};

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
        int minNumNeigh       = 1;  // Min neighbors required
        int maxNumNeigh       =-1;  // Max neighbors allowed (furthest events are discarded)
        int minDTperEvt       = 1;  // Min differential times per event pair required (Including P+S)
        int maxDTperEvt       =-1;  // Max differential times per event pair required (Including P+S)
        // From Waldhauser 2009: to assure a spatially homogeneous subsampling, reference
        // events are selected within each of five concentric, vertically longated
        // ellipsoidal layers of increasing thickness. Each layer has 8 quadrants.
        int numEllipsoids       = 5;
        double maxEllipsoidSize = 10; // km
    } step1Clustering;

    // cross correlation specific
    struct {
        std::string recordStreamURL;

        double minWeight      = 0;  // Min weight of phases required (0-1)
        double minEStoIEratio = 0;  // Min epi-sta to  interevent distance ration required
        double minESdist      = 0;  // Min epi-sta distance required
        double maxESdist      =-1;  // Max epi-sta distance allowed
        int minNumNeigh       = 1;  // Min neighbors required
        int maxNumNeigh       =-1;  // Max neighbors allowed (furthest events are discarded)
        int minDTperEvt       = 1;  // Min differential times per event pair required (Including P+S)
        int maxDTperEvt       =-1;  // Max differential times per event pair required (Including P+S)
        // From Waldhauser 2009: to assure a spatially homogeneous subsampling, reference
        // events are selected within each of five concentric, vertically longated
        // ellipsoidal layers of increasing thickness. Each layer has 8 quadrants.
        int numEllipsoids       = 5;
        double maxEllipsoidSize = 10; // km
    } step2Clustering;

    struct XCorr {
        double minCoef;     // Min xcorr coefficient required (0-1)
        double startOffset; // secs
        double endOffset;   // secs
        double maxDelay;    //secs
        std::vector<std::string> components; // priority list of components to use
    };
    std::map<std::string,struct XCorr> xcorr = {
        {"P", {}},
        {"S", {}}
    };

    // artificial phases
    struct {
        bool enable           = true;
    } artificialPhases;

    struct {
        std::string filterStr = "";
        double resampleFreq = 0;
        bool dump = false;
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

        CatalogCPtr getCatalog() { return _srcCat; }
        void setCatalog(const CatalogCPtr& catalog);

        CatalogPtr relocateCatalog(bool force = true, bool usePh2dt = false);
        CatalogPtr relocateSingleEvent(const CatalogCPtr& orgToRelocate);

        void setWorkingDirCleanup(bool cleanup) { _workingDirCleanup = cleanup; }
        bool workingDirCleanup() { return _workingDirCleanup; }

        void setUseCatalogDiskCache(bool cache) { _useCatalogDiskCache = cache; }
        bool useCatalogDiskCache() { return _useCatalogDiskCache; }

        static std::string relocationReport(const CatalogCPtr& relocatedEv);

    private:
        class XCorrCache;

        std::string generateWorkingSubDir(const Catalog::Event& ev) const;

        CatalogPtr relocateEventSingleStep(const CatalogCPtr& evToRelocateCat, const std::string& workingDir,
                                bool doXcorr, bool computeTheoreticalPhases,
                                std::string hypoddCtrlFile, double minPhaseWeight,
                                double minESdist, double maxESdist, double minEStoIEratio,
                                int minDTperEvt, int maxDTperEvt, int minNumNeigh, int maxNumNeigh,
                                int numEllipsoids, double maxEllipsoidSize);

        double computePickWeight(double uncertainty) const;
        double computePickWeight(const Catalog::Phase& phase) const;
        CatalogPtr filterPhasesAndSetWeights(const CatalogCPtr& catalog, const Catalog::Phase::Source& source,
                                   const std::vector<std::string>& PphaseToKeep,
                                   const std::vector<std::string>& SphaseToKeep) const;

        CatalogPtr selectNeighbouringEvents(const CatalogCPtr& catalog,
                                            const Catalog::Event& refEv,
                                            const CatalogCPtr& refEvCatalog,
                                            double minPhaseWeight = 0, double minESdis=0,
                                            double maxESdis=-1, double minEStoIEratio=0,
                                            int minDTperEvt=1, int maxDTperEvt=-1,
                                            int minNumNeigh=1, int maxNumNeigh=-1,
                                            int numEllipsoids=5, double maxEllipsoidSize=10,
                                            bool keepUnmatched=false, int *numNeigh=nullptr) const;
        std::map<unsigned,CatalogPtr> 
        selectNeighbouringEventsCatalog(const CatalogCPtr& catalog, double minPhaseWeight,
                                        double minESdis, double maxESdis, double minEStoIEratio,
                                        int minDTperEvt, int maxDTperEvt,
                                        int minNumNeigh, int maxNumNeigh,
                                        int numEllipsoids, double maxEllipsoidSize,
                                        bool keepUnmatched) const;

        void addMissingEventPhases(const CatalogCPtr& searchCatalog,
                                   const Catalog::Event& refEv,
                                   CatalogPtr& refEvCatalog);
        std::vector<Catalog::Phase> findMissingEventPhases(const CatalogCPtr& searchCatalog,
                                                           const Catalog::Event& refEv,
                                                           const CatalogPtr& refEvCatalog);
        typedef std::pair<std::string,std::string> MissingStationPhase;
        std::vector<MissingStationPhase> getMissingPhases(const CatalogCPtr& searchCatalog,
                                                          const Catalog::Event& refEv,
                                                          const CatalogPtr& refEvCatalog) const;
        typedef std::pair<Catalog::Event, Catalog::Phase> PhasePeer;
        std::vector<PhasePeer> findPhasePeers(const Catalog::Station& station, const std::string& phaseType,
                                              const CatalogCPtr& searchCatalog) const;
        Catalog::Phase createThoreticalPhase(const Catalog::Station& station,
                                             const std::string& phaseType,
                                             const Catalog::Event& refEv,
                                             const std::vector<HypoDD::PhasePeer>& peers,
                                             double phaseVelocity); 

        // Used to save xcorr results
        class XCorrCache {

        public:

            struct XCorrCacheEntry {

                double mean_coeff;
                double mean_lag;
                double min_lag;
                double max_lag;
                unsigned ccCount;

                typedef struct {
                    double coeff, lag, dtcc, weight, lowerUncertainty, upperUncertainty;
                } PeerInfo;
                std::map<unsigned,const PeerInfo> peers;

                std::string peersStr; // debug

                void update(const Catalog::Event& event, const Catalog::Phase& phase,
                            double coeff, double lag, double dtcc, double weight)
                {
                    PeerInfo pi= {coeff, lag, dtcc, weight, phase.lowerUncertainty, phase.upperUncertainty};
                    peers.insert( std::pair<unsigned,const PeerInfo>(event.id, pi) );
                    peersStr   += std::string(event) + " ";
                }

                void computeStats()
                {
                    ccCount = peers.size();
                    mean_coeff = 0;
                    mean_lag   = 0;
                    min_lag    = 0;
                    max_lag    = 0;
                    for ( auto& pair : peers )
                    {
                        const PeerInfo& data = pair.second;
                        mean_coeff += std::abs(data.coeff);
                        mean_lag   += data.lag;
                        min_lag    += data.lag - data.lowerUncertainty;
                        max_lag    += data.lag + data.upperUncertainty; 
                    }
                    mean_coeff /= ccCount;
                    mean_lag   /= ccCount;
                    min_lag    /= ccCount;
                    max_lag    /= ccCount; 
                }
            };

            XCorrCacheEntry& getForUpdate(unsigned evId, const std::string& stationId, const std::string& type)
            {
                std::string key = make_key(evId, stationId, type);
                return resultsByPhase[key];
            }

            void computeStats()
            {
                for ( auto& pair : resultsByPhase )  pair.second.computeStats();
            }

            bool has(unsigned evId, const std::string& stationId, const std::string& type ) const
            {
                std::string key = make_key(evId, stationId, type);
                return resultsByPhase.count(key) != 0;
            }

            const XCorrCacheEntry& get(unsigned evId, const std::string& stationId, const std::string& type ) const
            {
                std::string key = make_key(evId, stationId, type);
                return resultsByPhase.at(key);
            }

            bool has(unsigned evId1, unsigned evId2, const std::string& stationId, const std::string& type ) const
            {
                return has(evId1, stationId, type) && (get(evId1, stationId, type).peers.count(evId2) != 0);
            }

            const XCorrCacheEntry::PeerInfo& get(unsigned evId1, unsigned evId2, const std::string& stationId, const std::string& type ) const
            {
                return get(evId1, stationId, type).peers.at(evId2);
            } 

        private:
            static std::string make_key(unsigned evId, const std::string& stationId, const std::string& type )
            {
                return std::to_string(evId) + "." + stationId + "." + type; 
            }

            // cache of computed xcorr
            std::map<std::string, XCorrCacheEntry> resultsByPhase;

        };

        struct PhaseXCorrCfg {
            bool useDiskCache;
            WfMngr::WfCache* cache;
            bool allowSnrCheck;
        };

        XCorrCache buildXCorrCache(std::map<unsigned,CatalogPtr>& neighbourCats,
                                   bool computeTheoreticalPhases);
        XCorrCache buildXCorrCache(CatalogPtr& catalog, unsigned evToRelocateId,
                                   bool computeTheoreticalPhases);
        void buildXcorrDiffTTimePairs(CatalogPtr& catalog, const Catalog::Event& refEv,
                                      XCorrCache& xcorr);

        void fixPhases(CatalogPtr& catalog, const Catalog::Event& refEv, XCorrCache& xcorr);
 
        bool xcorrPhases(const Catalog::Event& event1, const Catalog::Phase& phase1, PhaseXCorrCfg& phCfg1,
                         const Catalog::Event& event2, const Catalog::Phase& phase2, PhaseXCorrCfg& phCfg2,
                         double& coeffOut, double& lagOut, double& diffTimeOut, double& weightOut);
        bool _xcorrPhases(const Catalog::Event& event1, const Catalog::Phase& phase1, PhaseXCorrCfg& phCfg1,
                          const Catalog::Event& event2, const Catalog::Phase& phase2, PhaseXCorrCfg& phCfg2,
                          double& coeffOut, double& lagOut, double& diffTimeOut, double& weightOut);
        bool xcorr(const GenericRecordCPtr& tr1, const GenericRecordCPtr& tr2, double maxDelay,
                   bool qualityCheck, double& delayOut, double& coeffOut) const;
        Core::TimeWindow xcorrTimeWindowLong(const Catalog::Phase& phase) const;
        Core::TimeWindow xcorrTimeWindowShort(const Catalog::Phase& phase) const;

        void printCounters();

        //
        // External HypoDD program specific
        //
        void runHypodd(const std::string& workingDir, const std::string& dtccFile,
                       const std::string& dtctFile, const std::string& eventFile,
                       const std::string& stationFile, const std::string& ctrlFile) const;
        void runPh2dt(const std::string& workingDir,
                      const std::string& stationFile,
                      const std::string& phaseFile) const;

        void createStationDatFile(const CatalogCPtr& catalog, const std::string& staFileName) const;
        void createPhaseDatFile(const CatalogCPtr& catalog, const std::string& phaseFileName) const;
        void createEventDatFile(const CatalogCPtr& catalog, const std::string& eventFileName) const;

        void createDtCt(std::map<unsigned,CatalogPtr>& neighbourCats, const std::string& dtctFile) const;
        void createDtCt(CatalogPtr& catalog, unsigned evToRelocateId, const std::string& dtctFile) const;
        void writeAbsTTimePairs(const CatalogCPtr& catalog, unsigned evToRelocateId, std::ofstream& outStream) const;

        std::set<unsigned> createDtCcPh2dt(const CatalogCPtr& catalog,
                                           const std::string& dtctFile,
                                           const std::string& dtccFile);
        void createDtCc(std::map<unsigned,CatalogPtr>& neighbourCats,
                        const std::string& dtccFile, const XCorrCache& xcorr);
        void createDtCc(CatalogPtr& catalog, unsigned evToRelocateId,
                         const std::string& dtccFile, const XCorrCache& xcorr);
        void writeXcorrDiffTTimePairs(CatalogPtr& catalog,
                                      unsigned evToRelocateId,
                                      const XCorrCache& xcorr,
                                      std::ofstream& outStream);

        CatalogPtr loadRelocatedCatalog(const CatalogCPtr& originalCatalog,
                                        const std::string& ddrelocFile,
                                        const std::string& ddresidualFile="") const;
    private:
        bool _workingDirCleanup = true;
        std::string _workingDir;
        std::string _cacheDir;
        std::string _wfDebugDir;

        CatalogCPtr _srcCat;
        CatalogCPtr _ddbgc;

        const Config _cfg;

        WfMngrPtr  _wf;
        bool _useCatalogDiskCache = false;
        WfMngr::WfCache _wfCache;

        struct {
            unsigned xcorr_performed;
            unsigned xcorr_performed_theo;
            unsigned xcorr_performed_s;
            unsigned xcorr_performed_s_theo;
            unsigned xcorr_good_cc;
            unsigned xcorr_good_cc_theo;
            unsigned xcorr_good_cc_s;
            unsigned xcorr_good_cc_s_theo;
        } mutable _counters;
};

}
}

#endif

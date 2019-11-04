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


namespace Seiscomp {
namespace HDD {

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
        double maxIEdist      =-1;  // Max interevent-distance allowed (km)
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
    };
    std::map<std::string,struct XCorr> xcorr = {
        {"P", {}},
        {"S", {}}
    };

    // artificial phases
    struct {
        bool enable           = false;
        bool fixAutoPhase     = false;
        double maxIEdist      = 5;
        unsigned numCC        = 2;
        double maxCCtw        = 10;
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

    private:
        double computePickWeight(double uncertainty) const;
        double computePickWeight(const Catalog::Phase& phase) const;
        CatalogPtr filterPhasesAndSetWeights(const CatalogCPtr& catalog,
                                   const std::vector<std::string>& PphaseToKeep,
                                   const std::vector<std::string>& SphaseToKeep) const;
        CatalogPtr createMissingPhases(const CatalogCPtr& catalog);
        void addMissingEventPhases(CatalogPtr& catalog, const Catalog::Event& refEv);
        std::vector<Catalog::Phase> findMissingEventPhases(const CatalogCPtr& catalog,
                                                           const Catalog::Event& refEv);
        void createStationDatFile(const CatalogCPtr& catalog, const std::string& staFileName) const;
        void createPhaseDatFile(const CatalogCPtr& catalog, const std::string& phaseFileName) const;
        void createEventDatFile(const CatalogCPtr& catalog, const std::string& eventFileName) const;
        void createDtCtCatalog(const CatalogCPtr& catalog,
                               const std::string& dtctFile) const;
        void createDtCtSingleEvent(const CatalogCPtr& catalog,
                                   unsigned evToRelocateId,
                                   const std::string& dtctFile) const;
        void buildAbsTTimePairs(const CatalogCPtr& catalog,
                                 unsigned evToRelocateId,
                                 std::ofstream& outStream) const;
        void createDtCcPh2dt(const CatalogCPtr& catalog,
                             const std::string& dtctFile,
                             const std::string& dtccFile);
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
        Core::TimeWindow xcorrTimeWindowLong(const Catalog::Phase& phase) const;
        Core::TimeWindow xcorrTimeWindowShort(const Catalog::Phase& phase) const;
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
                                            double minPhaseWeight = 0, double minESdis=0, double maxESdis=-1,
                                            double minEStoIEratio=0, double maxIEdis=-1,
                                            int minDTperEvt=1, int maxDTperEvt=-1,
                                            int minNumNeigh=1, int maxNumNeigh=-1,
                                            int numEllipsoids=5, int maxEllipsoidSize=0) const;
        std::map<unsigned,CatalogPtr> 
        selectNeighbouringEventsCatalog(const CatalogCPtr& catalog, double minPhaseWeight,
                                        double minESdis, double maxESdis,
                                        double minEStoIEratio, double maxIEdis,
                                        int minDTperEvt, int maxDTperEvt,
                                        int minNumNeigh, int maxNumNeigh,
                                        int numEllipsoids, int maxEllipsoidSize) const;
        bool xcorr(const GenericRecordCPtr& tr1, const GenericRecordCPtr& tr2, double maxDelay,
                   bool qualityCheck, double& delayOut, double& coeffOut) const;
        double S2Nratio(const GenericRecordCPtr& tr, const Core::Time& guidingPickTime,
                        double noiseOffsetStart, double noiseOffsetEnd,
                        double signalOffsetStart, double signalOffsetEnd) const;
        GenericRecordPtr getWaveform(const Core::TimeWindow& tw,
                                     const Catalog::Event& ev,
                                     const Catalog::Phase& ph,
                                     std::map<std::string,GenericRecordPtr>& memCache,
                                     bool useDiskCache,
                                     bool checkSnr);
        GenericRecordPtr loadProjectWaveform(const Core::TimeWindow& tw,
                                             const Catalog::Event& ev,
                                             const Catalog::Phase& ph,
                                             const DataModel::ThreeComponents& tc,
                                             const DataModel::SensorLocation *loc,
                                             bool useDiskCache) const;
        Core::TimeWindow traceTimeWindowToLoad(const Catalog::Phase& ph,
                                               const Core::TimeWindow& neededTW) const;
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
        void filter(GenericRecord &trace, bool demeaning=true, const std::string& filterStr="", double resampleFreq=0) const;
        void resample(GenericRecord& trace, double sf, bool average) const;
        std::string generateWorkingSubDir(const Catalog::Event& ev) const;
        std::string waveformFilename(const Catalog::Phase& ph, const Core::TimeWindow& tw) const;
        std::string waveformFilename(const std::string& networkCode, const std::string& stationCode,
                                     const std::string& locationCode, const std::string& channelCode,
                                     const Core::TimeWindow& tw) const;
        std::string waveformId(const Catalog::Phase& ph, const Core::TimeWindow& tw) const;
        std::string waveformId(const std::string& networkCode, const std::string& stationCode,
                               const std::string& locationCode, const std::string& channelCode,
                               const Core::TimeWindow& tw) const;

    private:
        std::string _workingDir;
        std::string _cacheDir;
        CatalogCPtr _srcCat;
        CatalogPtr _ddbgc;
        Config _cfg;
        bool _workingDirCleanup = true;
        bool _useCatalogDiskCache = false;
        std::map<std::string, GenericRecordPtr> _wfCache;
        std::set<std::string> _excludedWfs;

        struct {
            unsigned xcorr_tot;
            unsigned xcorr_performed;
            unsigned xcorr_cc_good;
            unsigned xcorr_cc_low;
            unsigned snr_low;
            unsigned wf_no_avail;
        } _counters;
};

}
}

#endif

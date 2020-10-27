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

#ifndef __HDD_HYPODD_H__
#define __HDD_HYPODD_H__

#include "catalog.h"
#include "clustering.h"
#include "solver.h"
#include "waveform.h"
#include "xcorrcache.ipp"

#include <seiscomp3/core/baseobject.h>
#include <seiscomp3/seismology/ttt.h>

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace Seiscomp {
namespace HDD {

struct Config
{

  std::vector<std::string> validPphases = {"Pg", "P", "Px"};
  std::vector<std::string> validSphases = {"Sg", "S", "Sx"};

  // Absolute travel time difference observations only
  struct
  {
    double minWeight = 0; // Min weight of phases required (0-1)
    double minEStoIEratio =
        0; // Min epi-sta to  interevent distance ration required
    double minESdist     = 0;  // Min epi-sta distance required
    double maxESdist     = -1; // Max epi-sta distance allowed
    unsigned minNumNeigh = 1;  // Min neighbors required
    unsigned maxNumNeigh =
        0; // Max neighbors allowed (furthest events are discarded)
    unsigned minDTperEvt =
        1; // Min differential times per event pair required (Including P+S)
    unsigned maxDTperEvt =
        0; // Max differential times per event pair required (Including P+S)
    // From Waldhauser 2009: to assure a spatially homogeneous subsampling,
    // reference events are selected within each of five concentric, vertically
    // longated ellipsoidal layers of increasing thickness. Each layer has 8
    // quadrants.
    unsigned numEllipsoids  = 5;
    double maxEllipsoidSize = 10; // km
  } ddObservations1;

  // Absolute travel time difference AND cross-correlation observations
  struct
  {
    double minWeight = 0; // Min weight of phases required (0-1)
    double minEStoIEratio =
        0; // Min epi-sta to  interevent distance ration required
    double minESdist     = 0;  // Min epi-sta distance required
    double maxESdist     = -1; // Max epi-sta distance allowed
    unsigned minNumNeigh = 1;  // Min neighbors required
    unsigned maxNumNeigh =
        0; // Max neighbors allowed (furthest events are discarded)
    unsigned minDTperEvt =
        1; // Min differential times per event pair required (Including P+S)
    unsigned maxDTperEvt =
        0; // Max differential times per event pair required (Including P+S)
    // From Waldhauser 2009: to assure a spatially homogeneous subsampling,
    // reference events are selected within each of five concentric, vertically
    // longated ellipsoidal layers of increasing thickness. Each layer has 8
    // quadrants.
    unsigned numEllipsoids  = 5;
    double maxEllipsoidSize = 10; // km

    //  cross-correlation specific
    double xcorrMaxEvStaDist   = -1; // max event to station distance
    double xcorrMaxInterEvDist = -1; // max inter-event distance
    std::string recordStreamURL;
  } ddObservations2;

  struct XCorr
  {
    double minCoef;                      // Min xcorr coefficient required (0-1)
    double startOffset;                  // secs
    double endOffset;                    // secs
    double maxDelay;                     // secs
    std::vector<std::string> components; // priority list of components to use
  };
  std::map<Catalog::Phase::Type, struct XCorr> xcorr = {
      {Catalog::Phase::Type::P, {}}, {Catalog::Phase::Type::S, {}}};

  struct
  {
    std::string filterStr = "";
    double resampleFreq   = 0;
  } wfFilter;

  struct
  {
    double minSnr      = 0;
    double noiseStart  = 0;
    double noiseEnd    = 0;
    double signalStart = 0;
    double signalEnd   = 0;
  } snr;

  struct
  {
    std::string type  = "LOCSAT";
    std::string model = "iasp91";
  } ttt;

  struct
  {
    std::string type                               = "LSMR"; // LSMR or LSQR
    bool L2normalization                           = true;
    unsigned solverIterations                      = 0;
    unsigned algoIterations                        = 20;
    double dampingFactorStart                      = 0.;
    double dampingFactorEnd                        = 0.;
    std::array<double, 4> meanShiftConstraintStart = {
        {0., 0., 0., 0.}}; // lon, lat, depth, time
    std::array<double, 4> meanShiftConstraintEnd = {
        {0., 0., 0., 0.}}; // lon, lat, depth, time
    double downWeightingByResidualStart = 0.;
    double downWeightingByResidualEnd   = 0.;
    bool usePickUncertainty             = false;
    double absTTDiffObsWeight           = 1.0;
    double xcorrObsWeight               = 1.0;
  } solver;
};

DEFINE_SMARTPOINTER(HypoDD);

class HypoDD : public Core::BaseObject
{

public:
  HypoDD(const CatalogCPtr &catalog,
         const Config &cfg,
         const std::string &workingDir);
  virtual ~HypoDD();

  void preloadData();

  CatalogCPtr getCatalog() { return _srcCat; }
  void setCatalog(const CatalogCPtr &catalog);

  CatalogPtr relocateCatalog();
  CatalogPtr relocateSingleEvent(const CatalogCPtr &orgToRelocate);
  void evalXCorr();

  void setWorkingDirCleanup(bool cleanup) { _workingDirCleanup = cleanup; }
  bool workingDirCleanup() const { return _workingDirCleanup; }

  void setUseCatalogWaveformDiskCache(bool cache);
  bool useCatalogWaveformDiskCache() const
  {
    return _useCatalogWaveformDiskCache;
  }

  void setWaveformCacheAll(bool all) { _waveformCacheAll = all; }
  bool waveformCacheAll() const { return _waveformCacheAll; }

  void setWaveformDebug(bool debug);
  bool waveformDebug() const { return _waveformDebug; }

  void setUseArtificialPhases(bool use) { _useArtificialPhases = use; }
  bool useArtificialPhases() const { return _useArtificialPhases; }

  static std::string relocationReport(const CatalogCPtr &relocatedEv);

private:
  void createWaveformCache();

  std::string generateWorkingSubDir(const Catalog::Event &ev) const;

  CatalogPtr relocateEventSingleStep(const CatalogCPtr bgCat,
                                     const CatalogCPtr &evToRelocateCat,
                                     const std::string &workingDir,
                                     bool doXcorr,
                                     bool computeTheoreticalPhases,
                                     double minPhaseWeight,
                                     double minESdist,
                                     double maxESdist,
                                     double minEStoIEratio,
                                     int minDTperEvt,
                                     int maxDTperEvt,
                                     int minNumNeigh,
                                     int maxNumNeigh,
                                     int numEllipsoids,
                                     double maxEllipsoidSize);

  CatalogPtr relocate(const CatalogCPtr &catalog,
                      const std::list<NeighboursPtr> &neighbourCats,
                      bool keepNeighboursFixed,
                      const XCorrCache &xcorr) const;

  struct ObservationParams
  {
    struct Entry
    {
      Catalog::Event event;
      Catalog::Station station;
      char phaseType;
      double travelTime;
      double takeOffAngle;
      double velocityAtSrc;
    };
    void add(TravelTimeTableInterfacePtr ttt,
             const Catalog::Event &event,
             const Catalog::Station &station,
             char phaseType);
    const Entry &
    get(unsigned eventId, const std::string stationId, char phaseType) const;
    void addToSolver(Solver &solver) const;

  private:
    std::unordered_map<std::string, Entry> _entries;
  };

  void addObservations(Solver &solver,
                       double absTTDiffObsWeight,
                       double xcorrObsWeight,
                       const CatalogCPtr &catalog,
                       const NeighboursPtr &neighbours,
                       bool keepNeighboursFixed,
                       const XCorrCache &xcorr,
                       ObservationParams &obsparams) const;

  CatalogPtr
  updateRelocatedEvents(const Solver &solver,
                        const CatalogCPtr &catalog,
                        const std::list<NeighboursPtr> &neighbourCats,
                        ObservationParams &obsparams) const;

  void updateRelocatedEventsFinalStats(
      const CatalogCPtr &startingCatalog,
      CatalogPtr &finalCatalog,
      const std::list<NeighboursPtr> &neighCluster) const;

  void addMissingEventPhases(const Catalog::Event &refEv,
                             CatalogPtr &refEvCatalog,
                             const CatalogCPtr &searchCatalog,
                             const NeighboursPtr &neighbours);

  std::vector<Catalog::Phase>
  findMissingEventPhases(const Catalog::Event &refEv,
                         CatalogPtr &refEvCatalog,
                         const CatalogCPtr &searchCatalog,
                         const NeighboursPtr &neighbours);

  typedef std::pair<std::string, Catalog::Phase::Type> MissingStationPhase;

  std::vector<MissingStationPhase>
  getMissingPhases(const Catalog::Event &refEv,
                   CatalogPtr &refEvCatalog,
                   const CatalogCPtr &searchCatalog) const;

  typedef std::pair<Catalog::Event, Catalog::Phase> PhasePeer;
  std::vector<PhasePeer> findPhasePeers(const Catalog::Station &station,
                                        const Catalog::Phase::Type &phaseType,
                                        const CatalogCPtr &searchCatalog,
                                        const NeighboursPtr &neighbours) const;

  Catalog::Phase
  createThoreticalPhase(const Catalog::Station &station,
                        const Catalog::Phase::Type &phaseType,
                        const Catalog::Event &refEv,
                        const std::vector<HypoDD::PhasePeer> &peers,
                        double phaseVelocity);

  XCorrCache buildXCorrCache(CatalogPtr &catalog,
                             const std::list<NeighboursPtr> &neighbourCats,
                             bool computeTheoreticalPhases);

  void buildXcorrDiffTTimePairs(CatalogPtr &catalog,
                                const NeighboursPtr &neighbours,
                                const Catalog::Event &refEv,
                                XCorrCache &xcorr);

  void fixPhases(CatalogPtr &catalog,
                 const Catalog::Event &refEv,
                 XCorrCache &xcorr);

  bool xcorrPhases(const Catalog::Event &event1,
                   const Catalog::Phase &phase1,
                   Waveform::LoaderPtr ph1Cache,
                   const Catalog::Event &event2,
                   const Catalog::Phase &phase2,
                   Waveform::LoaderPtr ph2Cache,
                   double &coeffOut,
                   double &lagOut);

  bool _xcorrPhases(const Catalog::Event &event1,
                    const Catalog::Phase &phase1,
                    Waveform::LoaderPtr ph1Cache,
                    const Catalog::Event &event2,
                    const Catalog::Phase &phase2,
                    Waveform::LoaderPtr ph2Cache,
                    double &coeffOut,
                    double &lagOut);

  Core::TimeWindow xcorrTimeWindowLong(const Catalog::Phase &phase) const;

  Core::TimeWindow xcorrTimeWindowShort(const Catalog::Phase &phase) const;

  GenericRecordCPtr getWaveform(const Core::TimeWindow &tw,
                                const Catalog::Event &ev,
                                const Catalog::Phase &ph,
                                Waveform::LoaderPtr wfLoader);

  void resetCounters();
  void printCounters() const;
  void updateCounters() const
  {
    updateCounters(_wfDiskCache, _wfSnrFilter, _wfMemCache);
  }
  void updateCounters(Waveform::LoaderPtr diskCache,
                      Waveform::SnrFilteredLoaderPtr snrFilter,
                      Waveform::LoaderPtr memCache) const;

private:
  bool _workingDirCleanup = true;
  std::string _workingDir;
  std::string _cacheDir;
  std::string _tmpCacheDir;
  std::string _wfDebugDir;

  CatalogCPtr _srcCat;
  CatalogCPtr _bgCat;

  const Config _cfg;

  bool _useCatalogWaveformDiskCache = true;
  bool _waveformCacheAll            = false;
  bool _waveformDebug               = false;

  bool _useArtificialPhases = true;

  TravelTimeTableInterfacePtr _ttt;

  Waveform::DiskCachedLoaderPtr _wfDiskCache;
  Waveform::SnrFilteredLoaderPtr _wfSnrFilter;
  Waveform::MemCachedLoaderPtr _wfMemCache;

  std::unordered_set<std::string> _unloadableWfs;

  struct
  {
    unsigned xcorr_performed;
    unsigned xcorr_performed_theo;
    unsigned xcorr_performed_s;
    unsigned xcorr_performed_s_theo;
    unsigned xcorr_good_cc;
    unsigned xcorr_good_cc_theo;
    unsigned xcorr_good_cc_s;
    unsigned xcorr_good_cc_s_theo;
    unsigned wf_downloaded;
    unsigned wf_no_avail;
    unsigned wf_disk_cached;
    unsigned wf_snr_low;
  } mutable _counters;

  // For waveforms that are cached to disk store at least DISK_TRACE_MIN_LEN
  // secs of data (centered at pick time)
  // This is to avoid re-downloading waveforms every time scrtdd is restarted
  // with a minimum change of the xcorr configuration, which happens when
  // the user is experimenting with the configuration optiobns.
  // This is a little overhead for the disk space but saves lot of precious user
  // time
  static constexpr const double DISK_TRACE_MIN_LEN = 10;
};

} // namespace HDD
} // namespace Seiscomp

#endif

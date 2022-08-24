/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as          *
 * published by the Free Software Foundation, either version 3 of the      *
 * License, or (at your option) any later version.                         *
 *                                                                         *
 * This software is distributed in the hope that it will be useful,        *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#ifndef __HDD_DD_H__
#define __HDD_DD_H__

#include "catalog.h"
#include "clustering.h"
#include "solver.h"
#include "timewindow.h"
#include "ttt.h"
#include "waveform.h"
#include "xcorrcache.h"

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace HDD {

struct Config
{
  std::vector<std::string> validPphases = {"Pg", "P", "Px"};
  std::vector<std::string> validSphases = {"Sg", "S", "Sx"};

  std::vector<std::pair<std::string, std::string>> compatibleChannels;

  // For waveforms that are cached to disk, store at least `diskTraceMinLen`
  // secs of data (centered at pick time).
  // This is to avoid re-downloading waveforms at every cross-correlation
  // configuration change
  double diskTraceMinLen = 10;

  struct XCorr
  {
    double minCoef;     // min cross-correlatation coefficient required (0-1)
    double startOffset; // secs
    double endOffset;   // secs
    double maxDelay;    // secs
    std::vector<std::string> components; // priority list of components to use
  };
  std::map<Catalog::Phase::Type, struct XCorr> xcorr = {
      {Catalog::Phase::Type::P, {0.50, -0.50, 0.50, 0.50, {"Z"}}},
      {Catalog::Phase::Type::S, {0.50, -0.50, 0.75, 0.50, {"H"}}}};

  struct
  {
    std::string filterStr = "ITAPER(1)>>BW_HLP(2,1,20)"; // "" -> no filtering
    double resampleFreq   = 0;                           // 0 -> no resampling
    // Extra waveform to load before and after the required window. This extra
    // len is useful to initialize the filter and to discard potential filter
    // artifacts at the beginning and end of trace
    double extraTraceLen = 1; // seconds
  } wfFilter;

  struct
  {
    double minSnr      = 2; // 0 -> no SNR check
    double noiseStart  = -3.0;
    double noiseEnd    = -0.350;
    double signalStart = -0.350;
    double signalEnd   = 1;
  } snr;
};

struct ClusteringOptions
{
  // min weight of phases required (0-1)
  double minWeight = 0;
  // min hypocenter-station to interevent distance ratio required
  double minEStoIEratio = 0;
  // min hypocenter-station distance required
  double minESdist = 0;
  // max hypocenter-station distance allowed
  double maxESdist = -1; // -1 -> disable
  // min neighbors required
  unsigned minNumNeigh = 1;
  // max neighbors allowed (furthest events are discarded)
  unsigned maxNumNeigh = 0; // 0 -> disable
  // min differential times per event pair required (Including P+S)
  unsigned minDTperEvt = 1;
  // max differential times per event pair required (Including P+S)
  unsigned maxDTperEvt = 0; // 0 -> disable
  // From Waldhauser 2009: to assure a spatially homogeneous subsampling,
  // reference events are selected within each of five concentric, vertically
  // longated ellipsoidal layers of increasing thickness. Each layer has 8
  // quadrants.
  unsigned numEllipsoids  = 5;
  double maxEllipsoidSize = 10; // km

  // cross-correlation observations specific
  double xcorrMaxEvStaDist = -1; // max event to station distance -1 -> disable
  double xcorrMaxInterEvDist = -1; // max inter-event distance -1 -> disable
  // use cross-correlation to detect phase for stations without one
  // (single-event only)
  bool xcorrDetectMissingPhases = false;
};

struct SolverOptions
{
  std::string type             = "LSMR"; // LSMR or LSQR
  bool L2normalization         = true;
  unsigned solverIterations    = 0; // 0 -> auto
  unsigned algoIterations      = 20;
  double absLocConstraintStart = 0; // 0 -> disable absolute location constraint
  double absLocConstraintEnd   = 0; // 0 -> disable absolute location constraint
  double dampingFactorStart    = 0.3;        // 0 -> disable damping factor
  double dampingFactorEnd      = 0.3;        // 0 -> disable damping factor
  double downWeightingByResidualStart = 10.; // 0 -> disbale downweighting
  double downWeightingByResidualEnd   = 3.;  // 0 -> disbale downweighting
  bool usePickUncertainty             = false;
  double absTTDiffObsWeight           = 0.5;
  double xcorrObsWeight               = 1.0;

  // Air-quakes are events whose depth shift above the range of the velocity
  // model (typically 0) during the inversion
  enum class AQ_ACTION
  {
    NONE,
    RESET,
    RESET_DEPTH
  };
  struct
  {
    // meters, threshold above which an event is considered an air-quake
    double elevationThreshold = 0;
    AQ_ACTION action          = AQ_ACTION::NONE;
  } airQuakes;
};

class DD
{

public:
  DD(const Catalog &catalog,
     const Config &cfg,
     std::unique_ptr<TravelTimeTable> ttt,
     std::unique_ptr<Waveform::Proxy> wf =
         std::unique_ptr<Waveform::Proxy>(new Waveform::NoWaveformProxy()));
  ~DD() = default;

  DD(const DD &other) = delete;
  DD &operator=(const DD &other) = delete;

  const Catalog &getCatalog() const { return _srcCat; }

  // Save relocation data into working directory: e.g. input catalog,
  // output catalog, found clusters, logs
  void disableSaveProcessing();
  void enableSaveProcessing(const std::string &workingDir);
  std::string saveProcessingDir() const { return _workingDir; }
  bool saveProcessing() const { return _saveProcessing; }

  // Enable/disable the usage of disk cache for the background
  // catalog waveforms
  void disableCatalogWaveformDiskCache();
  void enableCatalogWaveformDiskCache(const std::string &cacheDir);
  std::string catalogWaveformDiskCacheDir() const { return _cacheDir; }
  bool useCatalogWaveformDiskCache() const
  {
    return _useCatalogWaveformDiskCache;
  }

  // Enable/disable the usage of disk cache for temporary waveforms
  // (e.g. single-event). Normaly only background catalog waveforms
  // are cached to disk, but his might come in handy when testing or
  // developing and multiple relocations are attempted on the same
  // single-events
  void disableAllWaveformDiskCache();
  void enableAllWaveformDiskCache(const std::string &tmpCacheDir);
  std::string allWaveformDiskCacheDir() const { return _tmpCacheDir; }
  bool useAllWaveformCache() const { return _waveformCacheAll; }

  // preload all background catalog waveforms: store them on disk cache
  // (if enabled and not already there) then cache them in memory
  // already processed, ready for cross-correlation
  void preloadWaveforms();

  // free waveforms memory
  void unloadWaveforms();

  // save to disk the background catalog waveforms after being
  // resampled and filtered (for debugging)
  void dumpWaveforms(const std::string &basePath = "");

  // Find clusters in the background catalogs and return them
  std::list<Catalog> findClusters(const ClusteringOptions &clustOpt);

  // Multi-event relocation of the background catalog
  std::unique_ptr<Catalog>
  relocateMultiEvents(const ClusteringOptions &clustOpt,
                      const SolverOptions &solverOpt,
                      XCorrCache &precomputed);

  std::unique_ptr<Catalog>
  relocateMultiEvents(const ClusteringOptions &clustOpt,
                      const SolverOptions &solverOpt)
  {
    XCorrCache empty;
    return relocateMultiEvents(clustOpt, solverOpt, empty);
  }

  // Single-event relocation against background catalog
  std::unique_ptr<Catalog>
  relocateSingleEvent(const Catalog &singleEvent,
                      const ClusteringOptions &clustOpt1,
                      const ClusteringOptions &clustOpt2,
                      const SolverOptions &solverOpt);

  struct XCorrEvalStats
  {
    std::vector<unsigned> skipped;
    std::vector<unsigned> performed;
    std::vector<double> coeff;
    std::vector<double> lag;
    void summarize(unsigned &skipped,
                   unsigned &performed,
                   double &meanCoeff,
                   double &meanCoeffAbsDev,
                   double &medianCoeff,
                   double &medianCoeffAbsDev,
                   double &meanLag,
                   double &meanLagAbsDev,
                   double &medianLag,
                   double &medianLagAbsDev) const;
  };

  using evalXcorrCallback = std::function<void(
      const XCorrEvalStats &pTotStats,
      const XCorrEvalStats &sTotStats,
      const std::map<std::string, XCorrEvalStats> &pStatsByStation,
      const std::map<std::string, XCorrEvalStats> &sStatsByStation,
      const std::map<unsigned, XCorrEvalStats> &pStatsByStaDistance,
      const std::map<unsigned, XCorrEvalStats> &sStatsByStaDistance,
      const std::map<unsigned, XCorrEvalStats> &pStatsByInterEvDistance,
      const std::map<unsigned, XCorrEvalStats> &sStatsByInterEvDistance,
      double interEvDistStep,
      double staDistStep,
      double completionPercent)>;

  // Compute statistics on cross-correlatio results
  void evalXCorr(const ClusteringOptions &clustOpt,
                 const evalXcorrCallback &cb,
                 XCorrCache &precomputed,
                 const double interEvDistStep = 0.1, // km
                 const double staDistStep     = 3,   // km
                 bool theoretical             = false);

  // Compute statistics on cross-correlatio results
  void evalXCorr(const ClusteringOptions &clustOpt,
                 const evalXcorrCallback &cb,
                 const double interEvDistStep = 0.1, // km
                 const double staDistStep     = 3,   // km
                 bool theoretical             = false)
  {
    XCorrCache empty;
    evalXCorr(clustOpt, cb, empty, interEvDistStep, staDistStep, theoretical);
  }

  //
  // Static method
  //
  static std::string relocationReport(const Catalog &relocatedEv);

  static void xcorr(const Trace &tr1,
                    const Trace &tr2,
                    double maxDelay,
                    double &delayOut,
                    double &coeffOut);

private:
  void createWaveformCache();
  void
  replaceWaveformCacheLoader(const std::shared_ptr<Waveform::Loader> &baseLdr);

  std::string generateWorkingSubDir(const std::string &prefix) const;
  std::string generateWorkingSubDir(const Catalog::Event &ev) const;

  std::unique_ptr<Catalog>
  relocateEventSingleStep(const Catalog &bgCat,
                          const Catalog &evToRelocateCat,
                          const std::string &workingDir,
                          const ClusteringOptions &clustOpt,
                          const SolverOptions &solverOpt,
                          bool doXcorr);

  std::unique_ptr<Catalog>
  relocate(const Catalog &catalog,
           const std::unordered_map<unsigned, std::unique_ptr<Neighbours>>
               &neighCluster,
           const SolverOptions &solverOpt,
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
      double travelTimeResidual;
      double takeOfAngleAzim;
      double takeOfAngleDip;
      double velocityAtSrc;
      bool computeEvChanges;
    };
    bool add(HDD::TravelTimeTable &ttt,
             const Catalog::Event &event,
             const Catalog::Station &station,
             const Catalog::Phase &phase,
             bool computeEvChanges);
    const Entry &
    get(unsigned eventId, const std::string stationId, char phaseType) const;
    void addToSolver(Solver &solver) const;

  private:
    std::unordered_map<std::string, Entry> _entries;
  };

  void addObservations(Solver &solver,
                       double absTTDiffObsWeight,
                       double xcorrObsWeight,
                       const Catalog &catalog,
                       const Neighbours &neighbours,
                       bool keepNeighboursFixed,
                       bool usePickUncertainty,
                       const XCorrCache &xcorr,
                       ObservationParams &obsparams) const;

  std::unique_ptr<Catalog> updateRelocatedEvents(
      const Solver &solver,
      const Catalog &catalog,
      const SolverOptions &solverOpt,
      const std::unordered_map<unsigned, std::unique_ptr<Neighbours>>
          &neighCluster,
      ObservationParams &obsparams,
      double pickWeightScaler,
      std::unordered_map<unsigned, std::unique_ptr<Neighbours>>
          &finalNeighCluster) const;

  std::unique_ptr<Catalog> updateRelocatedEventsFinalStats(
      const Catalog &startingCatalog,
      const Catalog &finalCatalog,
      const std::unordered_map<unsigned, std::unique_ptr<Neighbours>>
          &neighCluster) const;

  void addMissingEventPhases(const Catalog::Event &refEv,
                             Catalog &refEvCatalog,
                             const Catalog &searchCatalog,
                             const Neighbours &neighbours);

  std::vector<Catalog::Phase>
  findMissingEventPhases(const Catalog::Event &refEv,
                         Catalog &refEvCatalog,
                         const Catalog &searchCatalog,
                         const Neighbours &neighbours);

  typedef std::pair<std::string, Catalog::Phase::Type> MissingStationPhase;

  std::vector<MissingStationPhase>
  getMissingPhases(const Catalog::Event &refEv,
                   Catalog &refEvCatalog,
                   const Catalog &searchCatalog) const;

  typedef std::pair<Catalog::Event, Catalog::Phase> PhasePeer;
  std::vector<PhasePeer> findPhasePeers(const Catalog::Station &station,
                                        const Catalog::Phase::Type &phaseType,
                                        const Catalog &searchCatalog,
                                        const Neighbours &neighbours) const;

  Catalog::Phase createThoreticalPhase(const Catalog::Station &station,
                                       const Catalog::Phase::Type &phaseType,
                                       const Catalog::Event &refEv,
                                       const std::vector<DD::PhasePeer> &peers,
                                       double phaseVelocity);

  XCorrCache buildXCorrCache(
      Catalog &catalog,
      const std::unordered_map<unsigned, std::unique_ptr<Neighbours>>
          &neighCluster,
      bool computeTheoreticalPhases,
      double xcorrMaxEvStaDist      = -1,
      double xcorrMaxInterEvDist    = -1,
      const XCorrCache &precomputed = XCorrCache());

  void buildXcorrDiffTTimePairs(Catalog &catalog,
                                const Neighbours &neighbours,
                                const Catalog::Event &refEv,
                                double xcorrMaxEvStaDist,   // -1 to disable
                                double xcorrMaxInterEvDist, // -1 to disable
                                XCorrCache &xcorr);

  void fixPhases(Catalog &catalog,
                 const Catalog::Event &refEv,
                 XCorrCache &xcorr) const;

  bool xcorrPhases(const Catalog::Event &event1,
                   const Catalog::Phase &phase1,
                   Waveform::Processor &ph1Cache,
                   const Catalog::Event &event2,
                   const Catalog::Phase &phase2,
                   Waveform::Processor &ph2Cache,
                   double &coeffOut,
                   double &lagOut,
                   std::string &componentOut);

  bool xcorrPhasesOneComponent(const Catalog::Event &event1,
                               const Catalog::Phase &phase1,
                               Waveform::Processor &ph1Cache,
                               const Catalog::Event &event2,
                               const Catalog::Phase &phase2,
                               Waveform::Processor &ph2Cache,
                               const std::string &component,
                               double &coeffOut,
                               double &lagOut);

  TimeWindow xcorrTimeWindowLong(const Catalog::Phase &phase) const;

  TimeWindow xcorrTimeWindowShort(const Catalog::Phase &phase) const;

  std::shared_ptr<const Trace> getWaveform(Waveform::Processor &wfLoader,
                                           const TimeWindow &tw,
                                           const Catalog::Event &ev,
                                           const Catalog::Phase &ph,
                                           const std::string &component);

  std::shared_ptr<Waveform::Processor>
  preloadNonCatalogWaveforms(Catalog &catalog,
                             const Neighbours &neighbours,
                             const Catalog::Event &refEv,
                             double xcorrMaxEvStaDist,
                             double xcorrMaxInterEvDist);

  void logXCorrSummary(const XCorrCache &xcorr);

  const std::vector<std::string>
  xcorrComponents(const Catalog::Phase &phase) const;

private:
  const Config _cfg;

  const Catalog _srcCat;
  const Catalog _bgCat;

  std::unique_ptr<TravelTimeTable> _ttt;
  std::shared_ptr<Waveform::Proxy> _wf;

  bool _saveProcessing = true;

  std::string _workingDir;
  std::string _cacheDir;
  std::string _tmpCacheDir;
  bool _useCatalogWaveformDiskCache = true;
  bool _waveformCacheAll            = false;

  struct
  {
    std::shared_ptr<Waveform::Loader> loader;
    std::shared_ptr<Waveform::DiskCachedLoader> diskCache;
    std::shared_ptr<Waveform::ExtraLenLoader> extraLen;
    std::shared_ptr<Waveform::SnrFilterPrc> snrFilter;
    std::shared_ptr<Waveform::MemCachedProc> memCache;
  } _wfAccess;
};

} // namespace HDD

#endif

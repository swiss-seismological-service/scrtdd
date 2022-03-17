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

  std::string recordStreamURL; // where to fetch waveforms from

  struct XCorr
  {
    double minCoef;     // min cross-correlatation coefficient required (0-1)
    double startOffset; // secs
    double endOffset;   // secs
    double maxDelay;    // secs
    std::vector<std::string> components; // priority list of components to use
  };
  std::map<Catalog::Phase::Type, struct XCorr> xcorr = {
      {Catalog::Phase::Type::P, {0.50, -0.50, 0.50, 0.350, {"Z"}}},
      {Catalog::Phase::Type::S, {0.50, -0.50, 0.75, 0.350, {"L"}}}};

  struct
  {
    std::string filterStr = "ITAPER(1)>>BW_HLP(2,1,20)"; // "" -> no filtering
    double resampleFreq   = 400;                         // 0 -> no resampling
  } wfFilter;

  struct
  {
    double minSnr      = 2; // 0 -> no SNR check
    double noiseStart  = -3.0;
    double noiseEnd    = -0.350;
    double signalStart = -0.350;
    double signalEnd   = 0.350;
  } snr;
};

struct ClusteringOptions
{
  // min weight of phases required (0-1)
  double minWeight = 0;
  // min hypocenter-station to interevent distance ration required
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
  double downWeightingByResidualEnd   = 6.;  // 0 -> disbale downweighting
  bool usePickUncertainty             = false;
  double absTTDiffObsWeight           = 0.5;
  double xcorrObsWeight               = 1.0;
};

class DD
{

public:
  DD(const Catalog &catalog,
     const Config &cfg,
     const std::string &workingDir,
     std::unique_ptr<HDD::TravelTimeTable> ttt);
  ~DD() = default;

  DD(const DD &other) = delete;
  DD &operator=(const DD &other) = delete;

  void preloadWaveforms();

  void unloadWaveforms() { createWaveformCache(); }

  void dumpWaveforms(const std::string &basePath = "");

  void dumpClusters(const ClusteringOptions &clustOpt,
                    const std::string &basePath = "");

  const Catalog &getCatalog() const { return _srcCat; }

  std::unique_ptr<Catalog>
  relocateMultiEvents(const ClusteringOptions &clustOpt,
                      const SolverOptions &solverOpt);
  std::unique_ptr<Catalog>
  relocateSingleEvent(const Catalog &singleEvent,
                      const ClusteringOptions &clustOpt1,
                      const ClusteringOptions &clustOpt2,
                      const SolverOptions &solverOpt);
  void evalXCorr(const ClusteringOptions &clustOpt, bool theoretical);

  void setSaveProcessing(bool dump) { _saveProcessing = dump; }
  bool saveProcessing() const { return _saveProcessing; }

  void setUseCatalogWaveformDiskCache(bool cache);
  bool useCatalogWaveformDiskCache() const
  {
    return _useCatalogWaveformDiskCache;
  }

  void setWaveformCacheAll(bool all) { _waveformCacheAll = all; }
  bool waveformCacheAll() const { return _waveformCacheAll; }

  void setUseArtificialPhases(bool use) { _useArtificialPhases = use; }
  bool useArtificialPhases() const { return _useArtificialPhases; }

  static std::string relocationReport(const Catalog &relocatedEv);

  static bool xcorr(const Trace &tr1,
                    const Trace &tr2,
                    double maxDelay,
                    bool qualityCheck,
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
                          bool doXcorr,
                          bool computeTheoreticalPhases);

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
      double xcorrMaxEvStaDist   = -1,
      double xcorrMaxInterEvDist = -1);

  void buildXcorrDiffTTimePairs(Catalog &catalog,
                                const Neighbours &neighbours,
                                const Catalog::Event &refEv,
                                double xcorrMaxEvStaDist,   // -1 to disable
                                double xcorrMaxInterEvDist, // -1 to disable
                                XCorrCache &xcorr);

  void
  fixPhases(Catalog &catalog, const Catalog::Event &refEv, XCorrCache &xcorr);

  bool xcorrPhases(const Catalog::Event &event1,
                   const Catalog::Phase &phase1,
                   Waveform::Processor &ph1Cache,
                   const Catalog::Event &event2,
                   const Catalog::Phase &phase2,
                   Waveform::Processor &ph2Cache,
                   double &coeffOut,
                   double &lagOut);

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

  void printCounters() const;

private:
  const Config _cfg;

  const Catalog _srcCat;
  const Catalog _bgCat;

  std::unique_ptr<HDD::TravelTimeTable> _ttt;

  bool _saveProcessing = true;

  const std::string _workingDir;
  const std::string _cacheDir;
  const std::string _tmpCacheDir;
  bool _useCatalogWaveformDiskCache = true;
  bool _waveformCacheAll            = false;
  bool _useArtificialPhases         = true;

  struct
  {
    std::shared_ptr<Waveform::Loader> loader;
    std::shared_ptr<Waveform::DiskCachedLoader> diskCache;
    std::shared_ptr<Waveform::ExtraLenLoader> extraLen;
    std::shared_ptr<Waveform::SnrFilterPrc> snrFilter;
    std::shared_ptr<Waveform::MemCachedProc> memCache;
  } _wfAccess;

  struct
  {
    unsigned xcorr_skipped;
    unsigned xcorr_performed;
    unsigned xcorr_performed_theo;
    unsigned xcorr_performed_s;
    unsigned xcorr_performed_s_theo;
    unsigned xcorr_good_cc;
    unsigned xcorr_good_cc_theo;
    unsigned xcorr_good_cc_s;
    unsigned xcorr_good_cc_s_theo;
    void reset() { *this = {0}; }
  } mutable _counters;

  // For waveforms that are cached to disk, store at least `DISK_TRACE_MIN_LEN`
  // secs of data (centered at pick time).
  // This is to avoid re-downloading waveforms every time the application is
  // restarted with a minimum change of the cross-correlation configuration,
  // which happens when the user is experimenting with the configuration
  // options.
  // Note that this approach requires slightly more disk space, but saves lot of
  // precious user time.
  static constexpr double DISK_TRACE_MIN_LEN = 10;
};

} // namespace HDD

#endif

/***************************************************************************
 * MIT License                                                             *
 *                                                                         *
 * Copyright (C) by ETHZ/SED                                               *
 *                                                                         *
 * Permission is hereby granted, free of charge, to any person obtaining a *
 * copy of this software and associated documentation files (the           *
 * “Software”), to deal in the Software without restriction, including     *
 * without limitation the rights to use, copy, modify, merge, publish,     *
 * distribute, sublicense, and/or sell copies of the Software, and to      *
 * permit persons to whom the Software is furnished to do so, subject to   *
 * the following conditions:                                               *
 *                                                                         *
 * The above copyright notice and this permission notice shall be          *
 * included in all copies or substantial portions of the Software.         *
 *                                                                         *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,         *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      *
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  *
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY    *
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,    *
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE       *
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                  *
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
  std::vector<std::string> validPphases      = {"Pg", "P", "Px"};
  std::vector<std::string> validSphases      = {"Sg", "S", "Sx"};
  std::vector<double> pickUncertaintyClasses = {0.000, 0.025, 0.050,
                                                0.100, 0.200, 0.400};
  bool PSTableOnly                           = true;
  struct
  {
    std::string filterStr = ""; // "" -> no filtering
    double resampleFreq   = 0;  // 0 -> no resampling
    // Extra waveform to load before and after the required window. This extra
    // len is useful to initialize the filter and to discard potential filter
    // artifacts at the beginning and end of trace
    double extraTraceLen = 1; // seconds
  } wfFilter;
};

struct ClusteringOptions
{
  // min hypocenter-station to interevent distance ratio required
  double minEStoIEratio = 5;
  // min hypocenter-station distance required
  double minESdist = 0;
  // max hypocenter-station distance allowed
  double maxESdist = -1; // -1 -> disable
  // min neighbors required
  unsigned minNumNeigh = 8;
  // max neighbors allowed (furthest events are discarded)
  unsigned maxNumNeigh = 40; // 0 -> disable
  // min phases per event pair required
  unsigned minNumPhases = 4;
  // max phases per event pair used
  unsigned maxNumPhases = 0; // 0 -> disable
  // max max neighbour distsance [km]
  double maxNeighbourDist = 5;
  // From Waldhauser 2009: to assure a spatially homogeneous subsampling,
  // reference events are selected within each of five concentric, vertically
  // longated ellipsoidal layers of increasing thickness. Each layer has 8
  // quadrants.
  unsigned numEllipsoids = 5;
};

struct XcorrOptions
{
  bool enable           = false; // perform or not cross-correlation
  double minEvStaDist   = 0;     // min event to station distance
  double maxEvStaDist   = -1;    // max event to station distance -1 -> disable
  double maxInterEvDist = -1;    // max inter-event distance -1 -> disable
  struct XCorr
  {
    double minCoef;     // min cross-correlatation coefficient required (0-1)
    double startOffset; // window start: secs before pick
    double endOffset;   // window end: secs after pick
    double winScaling;  // Window scaling coefficient
                        // WinLen = (endOff-startOff) + TravelTime * winScaling
    double maxDelay;    // secs
    std::vector<std::string> components; // priority list of components to use
  };
  std::map<Catalog::Phase::Type, struct XCorr> phase = {
      {Catalog::Phase::Type::P, {0.70, -0.50, 0.50, 0.02, 0.50, {"Z"}}},
      {Catalog::Phase::Type::S, {0.70, -0.50, 1.00, 0.04, 0.50, {"H"}}}};
};

struct SolverOptions
{
  std::string type                    = "LSMR"; // LSMR or LSQR
  bool L2normalization                = true;
  unsigned solverIterations           = 0; // 0 -> auto
  unsigned algoIterations             = 20;
  double absLocConstraintStart        = 0.3;  // 0 -> disable
  double absLocConstraintEnd          = 0.3;  // 0 -> disable
  double dampingFactorStart           = 0.01; // 0 -> disable
  double dampingFactorEnd             = 0.01; // 0 -> disable
  double downWeightingByResidualStart = 10.;  // 0 -> disable
  double downWeightingByResidualEnd   = 3.;   // 0 -> disable
  bool usePickUncertainties           = false;
  double xcorrWeightScaler            = 1.5;
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

  DD(const DD &other)            = delete;
  DD &operator=(const DD &other) = delete;

  const Catalog &getCatalog() const { return _bgCat; }

  //
  // Enable/disable the usage of disk cache for the background
  // catalog waveforms. Store at least `diskTraceMinLen`
  // secs of data (centered at pick time).
  //
  void enableCatalogWaveformDiskCache(const std::string &cacheDir,
                                      double diskTraceMinLen = 10.);
  void disableCatalogWaveformDiskCache();
  // Preload all background catalog waveforms: store them on disk cache
  // (if cache enabled)
  void preloadCatalogWaveformDiskCache(const XcorrOptions &xcorrOpt);

  std::string catalogWaveformDiskCacheDir() const { return _wfAccess.cacheDir; }
  bool useCatalogWaveformDiskCache() const { return _wfAccess.useCache; }
  double catalogWaveformDiskCacheTraceMinLen() const
  {
    return _wfAccess.diskTraceMinLen;
  }

  // free waveforms memory
  void unloadWaveforms();

  // save to disk the background catalog waveforms after being
  // resampled and filtered (for debugging)
  void dumpWaveforms(const XcorrOptions &xcorrOpt,
                     const std::string &basePath = "");

  // Find clusters in the background catalogs and return them
  std::list<std::unordered_map<unsigned, Neighbours>>
  findClusters(const ClusteringOptions &clustOpt);

  // Multi-event relocation of the background catalog
  Catalog relocateMultiEvents(
      std::list<std::unordered_map<unsigned, Neighbours>> &clusters,
      XCorrCache &xcorrData,
      const ClusteringOptions &clustOpt,
      const XcorrOptions &xcorrOpt,
      const SolverOptions &solverOpt,
      bool saveProcessing           = false,
      std::string processingDataDir = "");

  // Single-event relocation against background catalog
  Catalog relocateSingleEvent(const Catalog &singleEvent,
                              bool isManual,
                              const ClusteringOptions &clustOpt1,
                              const ClusteringOptions &clustOpt2,
                              const XcorrOptions &xcorrOpt,
                              const SolverOptions &solverOpt,
                              bool saveProcessing           = false,
                              std::string processingDataDir = "");
  //
  // Static method
  //
  static void xcorr(const Trace &tr1,
                    const Trace &tr2,
                    double maxDelay,
                    double &lagOut,
                    double &coeffOut);

private:
  void initWaveformAccess();
  void replaceWaveformLoader(const std::shared_ptr<Waveform::Loader> &baseLdr);

  std::string generateWorkingSubDir(const std::string &prefix) const;
  std::string generateWorkingSubDir(const Catalog::Event &ev) const;

  Catalog relocateEventSingleStep(const Catalog &bgCat,
                                  const Catalog &evToRelocateCat,
                                  const ClusteringOptions &clustOpt,
                                  const XcorrOptions &xcorrOpt,
                                  const SolverOptions &solverOpt,
                                  bool saveProcessing           = false,
                                  std::string processingDataDir = "");

  Catalog relocate(const Catalog &catalog,
                   const std::unordered_map<unsigned, Neighbours> &cluster,
                   const XcorrOptions &xcorrOpt,
                   const SolverOptions &solverOpt,
                   bool keepNeighboursFixed,
                   const XCorrCache &xcorr,
                   std::vector<Solver::DoubleDifference> startDDs,
                   std::vector<Solver::DoubleDifference> finalDDs) const;

  bool addObservations(Solver &solver,
                       const Catalog &catalog,
                       const Neighbours &neighbours,
                       bool keepNeighboursFixed,
                       bool usePickUncertainties,
                       double xcorrWeightScaler,
                       const XcorrOptions &xcorrOpt,
                       const XCorrCache &xcorr) const;

  bool addObservationParams(Solver &solver,
                            TravelTimeTable &ttt,
                            const Catalog::Event &event,
                            const Catalog::Station &station,
                            const Catalog::Phase &phase,
                            bool computeEvChanges) const;

  Catalog
  computeEventResiduals(const Solver &solver,
                        const Catalog &catalog,
                        const SolverOptions &solverOpt,
                        const std::unordered_map<unsigned, Neighbours> &cluster,
                        bool isFirstIteration) const;

  Catalog
  updateRelocatedEvents(const Solver &solver,
                        const Catalog &catalog,
                        const SolverOptions &solverOpt,
                        const std::unordered_map<unsigned, Neighbours> &cluster,
                        bool isFirstIteration) const;

  XCorrCache
  buildXCorrCache(const Catalog &catalog,
                  const std::unordered_map<unsigned, Neighbours> &cluster,
                  const XcorrOptions &xcorrOpt,
                  bool cacheRefEvWf,
                  const XCorrCache &precomputed = XCorrCache());

  void buildXcorrDiffTTimePairs(const Catalog &catalog,
                                const Neighbours &neighbours,
                                const Catalog::Event &refEv,
                                bool cacheRefEvWf,
                                const XcorrOptions &xcorrOpt,
                                const XCorrCache &precomputed,
                                XCorrCache &xcorr);

  bool xcorrPhases(const XcorrOptions &xcorrOpt,
                   const Catalog::Event &event1,
                   const Catalog::Phase &phase1,
                   Waveform::Processor &ph1Cache,
                   const Catalog::Event &event2,
                   const Catalog::Phase &phase2,
                   Waveform::Processor &ph2Cache,
                   double &coeffOut,
                   double &lagOut,
                   std::string &componentOut);

  bool xcorrPhasesOneComponent(const XcorrOptions &xcorrOpt,
                               const Catalog::Event &event1,
                               const Catalog::Phase &phase1,
                               Waveform::Processor &ph1Cache,
                               const Catalog::Event &event2,
                               const Catalog::Phase &phase2,
                               Waveform::Processor &ph2Cache,
                               const std::string &component,
                               double &coeffOut,
                               double &lagOut);

  TimeWindow xcorrTimeWindowLong(const XcorrOptions &xcorrOpt,
                                 const Catalog::Event &event,
                                 const Catalog::Phase &phase) const;

  TimeWindow xcorrTimeWindowShort(const XcorrOptions &xcorrOpt,
                                  const Catalog::Event &event,
                                  const Catalog::Phase &phase) const;

  std::shared_ptr<const Trace> getWaveform(Waveform::Processor &wfLoader,
                                           const TimeWindow &tw,
                                           const Catalog::Event &ev,
                                           const Catalog::Phase &ph,
                                           const std::string &component);

  std::shared_ptr<Waveform::Processor>
  preloadNonCatalogWaveforms(const Catalog &catalog,
                             const Neighbours &neighbours,
                             const Catalog::Event &refEv,
                             const XcorrOptions &xcorrOpt);

  void logXCorrSummary(const std::unordered_map<unsigned, Neighbours> &cluster,
                       const XcorrOptions &xcorrOpt,
                       const XCorrCache &xcorr);

  const std::vector<std::string>
  xcorrComponents(const XcorrOptions &xcorrOpt,
                  const Catalog::Phase &phase) const;

private:
  const Config _cfg;
  const Catalog _bgCat;
  const EventTree _evTree;
  std::unique_ptr<TravelTimeTable> _ttt;
  std::shared_ptr<Waveform::Proxy> _proxy;

  struct
  {
    bool useCache = false;
    std::string cacheDir;
    double diskTraceMinLen = 10;
    std::shared_ptr<Waveform::Loader> loader;
    std::shared_ptr<Waveform::DiskCachedLoader> diskCache;
    std::shared_ptr<Waveform::ExtraLenLoader> extraLen;
    std::shared_ptr<Waveform::BasicProcessor> basicProc;
    std::shared_ptr<Waveform::MemCachedProc> memCache;
  } _wfAccess;
};

} // namespace HDD

#endif

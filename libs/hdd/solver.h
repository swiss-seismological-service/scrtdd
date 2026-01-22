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

#ifndef __HDD_SOLVER_H__
#define __HDD_SOLVER_H__

#include "index.h"

#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace HDD {

/*
 * Store data for a double-difference problem as described in Waldhauser &
 * Ellsworth 2000 paper:
 *
 *      W G m = d W;
 *
 * Where G contains the partial derivatives of the travel times with respect to
 * event/station location and 1 in the column  corresponding to the origin
 * time correction term
 * m is a vector containing the changes in hypocentral parameters we wish to
 * determine for each event (delta x, delta y, delta z and delta travel time)
 * d is the data vector containing the double-differences
 * W is a diagonal matrix to weight each equation.
 *
 * This class also contains additional equations for constraining the shift of
 * earthquakes according to travel time residuals.
 *
 * We take advantage of the sparsness of G matrix, so G is not a full matrix.
 */
struct DDSystem
{
  // weight of each row of G matrix
  double *W = nullptr;
  // The G matrix stores data in a compressed format since it is a sparse
  // matrix: 3 partial derivatives for each event/station pair + tt (dx,dy,dz,1)
  double (*G)[4] = nullptr;
  // changes for each event hypocentral parameters we wish to determine
  // (x,y,z,t)
  double *m = nullptr;
  // double differences + optional travel time constraints
  double *d = nullptr;
  // L2 norm scaler for each G column
  double *L2NScaler = nullptr;
  // map of 2 event identifiers for each observation (index -1 means no
  // parameters)
  int *evByObs[2] = {nullptr, nullptr};
  // map of station identifiers for each observation
  unsigned *phStaByObs = nullptr;

  // number of observations
  unsigned nObs = 0;
  // number of events
  unsigned nEvts = 0;
  // number of stations
  unsigned nPhStas = 0;
  // number of obtional travel time constraints
  unsigned nTTconstraints = 0;

  unsigned numColsG           = 0;
  unsigned numRowsG           = 0;
  unsigned numCompressedColsG = 0;
  unsigned numCompressedRowsG = 0;

  DDSystem(unsigned _nObs,
           unsigned _nEvts,
           unsigned _nPhStas,
           unsigned _nTTconstraints = 0)
      : nObs(_nObs), nEvts(_nEvts), nPhStas(_nPhStas),
        nTTconstraints(_nTTconstraints), numColsG(nEvts * 4),
        numRowsG(nObs + nTTconstraints), numCompressedColsG(4),
        numCompressedRowsG(nEvts * nPhStas)
  {
    _W          = std::vector<double>(numRowsG);
    _G          = std::vector<double[4]>(numCompressedRowsG);
    _m          = std::vector<double>(numColsG);
    _d          = std::vector<double>(numRowsG);
    _L2NScaler  = std::vector<double>(numColsG);
    _evByObs[0] = std::vector<int>(numRowsG);
    _evByObs[1] = std::vector<int>(numRowsG);
    _phStaByObs = std::vector<unsigned>(numRowsG);

    W          = _W.data();
    G          = _G.data();
    m          = _m.data();
    d          = _d.data();
    L2NScaler  = _L2NScaler.data();
    evByObs[0] = _evByObs[0].data();
    evByObs[1] = _evByObs[1].data();
    phStaByObs = _phStaByObs.data();
  }

  DDSystem() = default;

  DDSystem(const DDSystem &other)            = default;
  DDSystem &operator=(const DDSystem &other) = default;

  DDSystem(DDSystem &&other)            = default;
  DDSystem &operator=(DDSystem &&other) = default;

  ~DDSystem() = default;

private:
  std::vector<double> _W;
  std::vector<double[4]> _G;
  std::vector<double> _m;
  std::vector<double> _d;
  std::vector<double> _L2NScaler;
  std::vector<int> _evByObs[2];
  std::vector<unsigned> _phStaByObs;
};

/*
 * Solver for double difference problems.
 *
 * For details, see Waldhauser & Ellsworth 2000 paper.
 */
class Solver
{

public:
  Solver(std::string type) : _type(type) {}
  ~Solver() = default;

  Solver(const Solver &other)           = delete;
  Solver operator=(const Solver &other) = delete;

  void addEvent(unsigned evId, double evLat, double evLon, double evDepth);

  void addStation(const std::string &staId,
                  double staLat,
                  double staLon,
                  double staElevation);

  void addObservation(unsigned evId1,
                      unsigned evId2,
                      const std::string &staId,
                      const std::string &phase,
                      double timeDiff,
                      double aPrioriWeight,
                      bool xcorrUsed,
                      double xcorrCoeff);

  void addObservationParams(unsigned evId,
                            const std::string &staId,
                            const std::string &phase,
                            bool computeEvChanges,
                            double travelTime,
                            double travelTimeResidual,
                            double dx,
                            double dy,
                            double dz);

  bool
  getEvent(unsigned evId, double &evLat, double &evLon, double &evDepth) const;

  bool getStation(const std::string &staId,
                  double &staLat,
                  double &staLon,
                  double &staElevation) const;

  bool getObservationParams(unsigned evId,
                            const std::string &staId,
                            const std::string &phase,
                            bool &computeEvChanges,
                            double &travelTime,
                            double &travelTimeResidual,
                            double &dx,
                            double &dy,
                            double &dz) const;

  void prepare(double ttConstraint = 0, double residualDownWeight = 0);

  struct DoubleDifference
  {
    unsigned evId1;
    unsigned evId2;
    std::string staId;
    std::string phase;
    double weight;
    bool xcorrUsed;
    double xcorrCoeff;
    double observedTimeDiff;   // sec
    double computedTimeDiff;   // sec
    double doubleDifference;   // sec
    double interEventDistance; // km
  };
  std::vector<DoubleDifference> getDoubleDifferences() const;

  void solve(unsigned numIterations = 0,
             double dampingFactor   = 0,
             bool normalizeG        = true);

  bool getEventChanges(unsigned evId,
                       double &deltaX,           // km
                       double &deltaY,           // km
                       double &deltaZ,           // km
                       double &deltaTime) const; // sec

  bool isEventPhaseUsed(unsigned evId,
                        const std::string &staId,
                        const std::string &phase) const;

private:
  std::multimap<double, unsigned> computeInterEventDistance();

  std::vector<double>
  computeResidualWeights(const std::vector<double> &residuals,
                         const double alpha) const;

  template <typename T>
  void _solve(unsigned numIterations,
              double dampingFactor,
              bool normalizeG,
              std::set<unsigned> rejectStoppingReasons);

  void loadSolutions();

private:
  IdToIndex<unsigned> _eventIdConverter;
  IdToIndex<std::string> _staIdConverter;
  IdToIndex<std::string> _phStaIdConverter;
  IdToIndex<std::string> _obsIdConverter;

  struct Observation
  {
    // user provided
    unsigned ev1Idx;
    unsigned ev2Idx;
    unsigned phStaIdx;
    double timeDiff;
    double aPrioriWeight;
    bool xcorrUsed;
    double xcorrCoeff;
    // computed
    double computedTimeDiff;
    double doubleDifference;
    double weight;
    double interEventDistance;
  };
  std::unordered_map<unsigned, Observation> _observations; // key = obsIdx

  struct EventParams
  {
    double lat, lon, depth;
  };
  std::unordered_map<unsigned, EventParams> _eventParams; // key = evIdx

  struct StationParams
  {
    double lat, lon, elevation;
  };
  std::unordered_map<unsigned, StationParams> _stationParams; // key = staIdx

  struct ObservationParams
  {
    bool computeEvChanges;
    double travelTime;         // sec
    double travelTimeResidual; // sec
    double dx;                 // sec/km
    double dy;                 // sec/km
    double dz;                 // sec/km
  };
  // key1=evIdx  key2=phStaIdx
  std::unordered_map<unsigned, std::unordered_map<unsigned, ObservationParams>>
      _obsParams;

  struct Stats
  {
    unsigned finalTotalObs  = 0;
    double totalFinalWeight = 0;
  };
  // key1=evIdx  key2=phStaIdx
  std::unordered_map<unsigned, std::unordered_map<unsigned, Stats>> _stats;

  struct EventDeltas
  {
    double time; // sec
    double x;    // km
    double y;    // km
    double z;    // km
  };
  std::unordered_map<unsigned, EventDeltas> _eventDeltas; // key = evIdx

  DDSystem _dd;
  std::string _type;
};

} // namespace HDD

#endif

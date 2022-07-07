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
 * event location and origin times.
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

  // number of observations
  const unsigned nObs;
  // number of events
  const unsigned nEvts;
  // number of stations
  const unsigned nPhStas;
  // number of obtional travel time constraints
  const unsigned nTTconstraints;
  // weight of each row of G matrix
  double *W;
  // The G matrix stores data in a dense format since it is a sparse matrix:
  // 3 partial derivatives for each event/station pair + tt (dx,dy,dz,1)
  double (*G)[4];
  // changes for each event hypocentral parameters we wish to determine
  // (x,y,z,t)
  double *m;
  // double differences + optional travel time constraints
  double *d;
  // L2 norm scaler for each G column
  double *L2NScaler;
  // map of 2 event identifiers for each observation (index -1 means no
  // parameters)
  int *evByObs[2];
  // map of station identifiers for each observation
  unsigned *phStaByObs;

  const unsigned numColsG;
  const unsigned numRowsG;

  DDSystem(unsigned _nObs,
           unsigned _nEvts,
           unsigned _nPhStas,
           unsigned _nTTconstraints = 0)
      : nObs(_nObs), nEvts(_nEvts), nPhStas(_nPhStas),
        nTTconstraints(_nTTconstraints), numColsG(nEvts * 4),
        numRowsG(_nObs + _nTTconstraints)
  {
    W          = new double[numRowsG];
    G          = new double[nEvts * nPhStas][4];
    m          = new double[numColsG];
    d          = new double[numRowsG];
    L2NScaler  = new double[numColsG];
    evByObs[0] = new int[numRowsG];
    evByObs[1] = new int[numRowsG];
    phStaByObs = new unsigned[numRowsG];
  }

  ~DDSystem()
  {
    delete[] phStaByObs;
    delete[] evByObs[0];
    delete[] evByObs[1];
    delete[] L2NScaler;
    delete[] d;
    delete[] m;
    delete[] G;
    delete[] W;
  }

  DDSystem(const DDSystem &other) = delete;
  DDSystem operator=(const DDSystem &other) = delete;
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

  Solver(const Solver &other) = delete;
  Solver operator=(const Solver &other) = delete;

  void addObservation(unsigned evId1,
                      unsigned evId2,
                      const std::string &staId,
                      char phase,
                      double diffTime,
                      double aPrioriWeight,
                      bool isXcorr);

  void addObservationParams(unsigned evId,
                            const std::string &staId,
                            char phase,
                            double evLat,
                            double evLon,
                            double evDepth,
                            double staLat,
                            double staLon,
                            double staElevation,
                            bool computeEvChanges,
                            double travelTime,
                            double travelTimeResidual,
                            double takeOffAngleAzim,
                            double takeOffAngleDip,
                            double velocityAtSrc);

  void solve(unsigned numIterations    = 0,
             double ttConstraint       = 0,
             double dampingFactor      = 0,
             double residualDownWeight = 0,
             bool normalizeG           = true);

  bool getEventChanges(unsigned evId,
                       double &deltaLat,
                       double &deltaLon,
                       double &deltaDepth,
                       double &deltaTT) const;

  bool getObservationParamsChanges(unsigned evId,
                                   const std::string &staId,
                                   char phase,
                                   unsigned &startingTTObs,
                                   unsigned &startingCCObs,
                                   unsigned &finalTotalObs,
                                   double &meanAPrioriWeight,
                                   double &meanFinalWeight,
                                   double &meanObsResidual,
                                   std::set<unsigned> &evIds) const;

private:
  void computePartialDerivatives();

  std::multimap<double, unsigned> computeInterEventDistance() const;

  std::vector<double>
  computeResidualWeights(const std::vector<double> &residuals,
                         const double alpha) const;

  void prepareDDSystem(double ttConstraint,
                       double dampingFactor,
                       double residualDownWeight);

  template <typename T>
  void _solve(unsigned numIterations,
              double ttConstraint,
              double dampingFactor,
              double residualDownWeight,
              bool normalizeG);

  void loadSolutions();

private:
  IdToIndex<unsigned> _eventIdConverter;
  IdToIndex<std::string> _phStaIdConverter;
  IdToIndex<std::string> _obsIdConverter;

  struct Observation
  {
    unsigned ev1Idx;
    unsigned ev2Idx;
    unsigned phStaIdx;
    double observedDiffTime;
    double aPrioriWeight;
    bool isXcorr;
  };
  std::unordered_map<unsigned, Observation> _observations; // key = obsIdx

  struct EventParams
  {
    double lat, lon, depth;
    double x, y, z; // km
  };
  std::unordered_map<unsigned, EventParams> _eventParams; // key = evIdx

  struct StationParams
  {
    double lat, lon, elevation;
    double x, y, z; // km
  };
  std::unordered_map<unsigned, StationParams> _stationParams; // key = phStaIdx

  struct ObservationParams
  {
    bool computeEvChanges;
    double travelTime;
    double travelTimeResidual;
    double takeOffAngleAzim;
    double takeOffAngleDip;
    double velocityAtSrc;
    double dx;
    double dy;
    double dz;
  };
  // key1=evIdx  key2=phStaIdx
  std::unordered_map<unsigned, std::unordered_map<unsigned, ObservationParams>>
      _obsParams;

  struct ParamStats
  {
    unsigned startingTTObs    = 0;
    unsigned startingCCObs    = 0;
    unsigned finalTotalObs    = 0;
    double totalAPrioriWeight = 0;
    double totalFinalWeight   = 0;
    double totalResiduals     = 0;
    std::set<unsigned> peerEvIds;
  };
  // key1=evIdx  key2=phStaIdx
  std::unordered_map<unsigned, std::unordered_map<unsigned, ParamStats>>
      _paramStats;

  struct
  {
    double lat, lon, depth;
  } _centroid;

  struct EventDeltas
  {
    double time;      // sec
    double depth;     // km
    double latitude;  // km
    double longitude; // km
  };
  std::unordered_map<unsigned, EventDeltas> _eventDeltas; // key = evIdx

  std::vector<double> _residuals;
  std::unique_ptr<DDSystem> _dd;
  std::string _type;
};

} // namespace HDD

#endif

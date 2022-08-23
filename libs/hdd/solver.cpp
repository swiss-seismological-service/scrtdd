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

#include "solver.h"
#include "log.h"
#include "lsmr.h"
#include "lsqr.h"
#include "utils.h"

#include <sstream>

using namespace std;
using HDD::degToRad;
using HDD::Exception;
using HDD::radToDeg;
using HDD::square;

namespace {

/**
 * Common DDSystem adapter for both LSQR and LSMR solvers
 * T can be `lsqrBase` or `lsmrBase`.
 */
template <typename T> class Adapter : public T
{
private:
  HDD::DDSystem &_dd; // doesn't own the DDSystem

public:
  Adapter(HDD::DDSystem &dd) : _dd(dd) {}
  virtual ~Adapter() = default;

  /*
   * Scale G by normalizing the L2-norm of each column as suggested
   * by LSQR and LSMR solvers.
   */
  void L2normalize()
  {
    std::fill_n(_dd.L2NScaler, _dd.numColsG, 0.);

    for (unsigned int ob = 0; ob < _dd.numRowsG; ob++)
    {
      const double obsW = _dd.W[ob];
      if (obsW == 0.) continue;

      const unsigned phStaIdx =
          _dd.phStaByObs[ob]; // station for this observation

      const int evIdx1 = _dd.evByObs[0][ob]; // event 1 for this observation
      if (evIdx1 >= 0)
      {
        const unsigned idxG     = evIdx1 * _dd.nPhStas + phStaIdx;
        const unsigned evOffset = evIdx1 * 4;
        _dd.L2NScaler[evOffset + 0] += square(_dd.G[idxG][0] * obsW);
        _dd.L2NScaler[evOffset + 1] += square(_dd.G[idxG][1] * obsW);
        _dd.L2NScaler[evOffset + 2] += square(_dd.G[idxG][2] * obsW);
        _dd.L2NScaler[evOffset + 3] += square(_dd.G[idxG][3] * obsW);
      }

      const int evIdx2 = _dd.evByObs[1][ob]; // event 2 for this observation
      if (evIdx2 >= 0)
      {
        const unsigned idxG     = evIdx2 * _dd.nPhStas + phStaIdx;
        const unsigned evOffset = evIdx2 * 4;
        _dd.L2NScaler[evOffset + 0] += square(_dd.G[idxG][0] * obsW);
        _dd.L2NScaler[evOffset + 1] += square(_dd.G[idxG][1] * obsW);
        _dd.L2NScaler[evOffset + 2] += square(_dd.G[idxG][2] * obsW);
        _dd.L2NScaler[evOffset + 3] += square(_dd.G[idxG][3] * obsW);
      }
    }

    for (unsigned col = 0; col < _dd.numColsG; col++)
    {
      _dd.L2NScaler[col] = 1. / std::sqrt(_dd.L2NScaler[col]);
    }
  }

  /*
   * Rescale m back to the initial scaling.
   */
  void L2DeNormalize()
  {
    for (unsigned evOffset = 0; evOffset < _dd.numColsG; evOffset += 4)
    {
      _dd.m[evOffset + 0] *= _dd.L2NScaler[evOffset + 0];
      _dd.m[evOffset + 1] *= _dd.L2NScaler[evOffset + 1];
      _dd.m[evOffset + 2] *= _dd.L2NScaler[evOffset + 2];
      _dd.m[evOffset + 3] *= _dd.L2NScaler[evOffset + 3];
    }
  }

  /**
   * Required by `lsqrBase` and `lsmrBase`:
   *
   * computes y = y + A*x without altering x,
   * where A is a matrix of dimensions A[m][n].
   * The size of the vector x is n.
   * The size of the vector y is m.
   */
  void Aprod1(unsigned int m, unsigned int n, const double *x, double *y) const
  {
    if (m != _dd.numRowsG || n != _dd.numColsG)
    {
      string msg =
          HDD::strf("Solver: Internal logic error (m=%u n=%u but G=%ux%u)", m,
                    n, _dd.numRowsG, _dd.numColsG);
      throw Exception(msg);
    }

    for (unsigned int ob = 0; ob < _dd.numRowsG; ob++)
    {
      if (_dd.W[ob] == 0.) continue;

      const unsigned phStaIdx =
          _dd.phStaByObs[ob]; // station for this observation
      double sum = 0;

      const int evIdx1 = _dd.evByObs[0][ob]; // event 1 for this observation
      if (evIdx1 >= 0)
      {
        const unsigned idxG     = evIdx1 * _dd.nPhStas + phStaIdx;
        const unsigned evOffset = evIdx1 * 4;
        sum += _dd.G[idxG][0] * _dd.L2NScaler[evOffset + 0] * x[evOffset + 0];
        sum += _dd.G[idxG][1] * _dd.L2NScaler[evOffset + 1] * x[evOffset + 1];
        sum += _dd.G[idxG][2] * _dd.L2NScaler[evOffset + 2] * x[evOffset + 2];
        sum += _dd.G[idxG][3] * _dd.L2NScaler[evOffset + 3] * x[evOffset + 3];
      }

      const int evIdx2 = _dd.evByObs[1][ob]; // event 2 for this observation
      if (evIdx2 >= 0)
      {
        const unsigned idxG     = evIdx2 * _dd.nPhStas + phStaIdx;
        const unsigned evOffset = evIdx2 * 4;
        sum -= _dd.G[idxG][0] * _dd.L2NScaler[evOffset + 0] * x[evOffset + 0];
        sum -= _dd.G[idxG][1] * _dd.L2NScaler[evOffset + 1] * x[evOffset + 1];
        sum -= _dd.G[idxG][2] * _dd.L2NScaler[evOffset + 2] * x[evOffset + 2];
        sum -= _dd.G[idxG][3] * _dd.L2NScaler[evOffset + 3] * x[evOffset + 3];
      }

      y[ob] += _dd.W[ob] * sum;
    }
  }

  /**
   * Required by `lsqrBase` and `lsmrBase`:
   *
   * computes x = x + A'*y without altering y,
   * where A is a matrix of dimensions A[m][n].
   * The size of the vector x is n.
   * The size of the vector y is m.
   */
  void Aprod2(unsigned int m, unsigned int n, double *x, const double *y) const
  {
    if (m != _dd.numRowsG || n != _dd.numColsG)
    {
      string msg =
          HDD::strf("Solver: Internal logic error (m=%u n=%u but G=%ux%u)", m,
                    n, _dd.numRowsG, _dd.numColsG);
      throw Exception(msg);
    }

    for (unsigned int ob = 0; ob < _dd.numRowsG; ob++)
    {
      const double wY = y[ob] * _dd.W[ob];
      if (wY == 0.) continue;

      const unsigned phStaIdx =
          _dd.phStaByObs[ob]; // station for this observation

      const int evIdx1 = _dd.evByObs[0][ob]; // event 1 for this observation
      if (evIdx1 >= 0)
      {
        const unsigned idxG     = evIdx1 * _dd.nPhStas + phStaIdx;
        const unsigned evOffset = evIdx1 * 4;
        x[evOffset + 0] += _dd.G[idxG][0] * _dd.L2NScaler[evOffset + 0] * wY;
        x[evOffset + 1] += _dd.G[idxG][1] * _dd.L2NScaler[evOffset + 1] * wY;
        x[evOffset + 2] += _dd.G[idxG][2] * _dd.L2NScaler[evOffset + 2] * wY;
        x[evOffset + 3] += _dd.G[idxG][3] * _dd.L2NScaler[evOffset + 3] * wY;
      }

      const int evIdx2 = _dd.evByObs[1][ob]; // event 2 for this observation
      if (evIdx2 >= 0)
      {
        const unsigned idxG     = evIdx2 * _dd.nPhStas + phStaIdx;
        const unsigned evOffset = evIdx2 * 4;
        x[evOffset + 0] -= _dd.G[idxG][0] * _dd.L2NScaler[evOffset + 0] * wY;
        x[evOffset + 1] -= _dd.G[idxG][1] * _dd.L2NScaler[evOffset + 1] * wY;
        x[evOffset + 2] -= _dd.G[idxG][2] * _dd.L2NScaler[evOffset + 2] * wY;
        x[evOffset + 3] -= _dd.G[idxG][3] * _dd.L2NScaler[evOffset + 3] * wY;
      }
    }
  }
};

} // namespace

namespace HDD {

void Solver::addObservation(unsigned evId1,
                            unsigned evId2,
                            const std::string &staId,
                            char phase,
                            double observedDiffTime,
                            double aPrioriWeight,
                            bool isXcorr)
{
  string phStaId    = string(1, phase) + "@" + staId;
  string obsId      = to_string(evId1) + "+" + to_string(evId2) + "_" + phStaId;
  unsigned evIdx1   = _eventIdConverter.convert(evId1);
  unsigned evIdx2   = _eventIdConverter.convert(evId2);
  unsigned phStaIdx = _phStaIdConverter.convert(phStaId);
  unsigned obsIdx   = _obsIdConverter.convert(obsId);
  _observations[obsIdx] = Observation{evIdx1,           evIdx2,        phStaIdx,
                                      observedDiffTime, aPrioriWeight, isXcorr};
}

void Solver::addObservationParams(unsigned evId,
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
                                  double velocityAtSrc)
{
  string phStaId      = string(1, phase) + "@" + staId;
  int evIdx           = _eventIdConverter.convert(evId);
  unsigned phStaIdx   = _phStaIdConverter.convert(phStaId);
  _eventParams[evIdx] = EventParams{evLat, evLon, evDepth, 0, 0, 0};
  _stationParams[phStaIdx] =
      StationParams{staLat, staLon, staElevation, 0, 0, 0};
  _obsParams[evIdx][phStaIdx] = ObservationParams{computeEvChanges,
                                                  travelTime,
                                                  travelTimeResidual,
                                                  takeOffAngleAzim,
                                                  takeOffAngleDip,
                                                  velocityAtSrc,
                                                  0,
                                                  0,
                                                  0};
}

bool Solver::getEventChanges(unsigned evId,
                             double &deltaLat,
                             double &deltaLon,
                             double &deltaDepth,
                             double &deltaTT) const
{
  unsigned evIdx;
  if (!_eventIdConverter.hasId(evId, evIdx)) return false;

  if (_eventDeltas.find(evIdx) == _eventDeltas.end()) return false;

  const EventDeltas &evDelta = _eventDeltas.at(evIdx);
  deltaLat                   = evDelta.latitude;
  deltaLon                   = evDelta.longitude;
  deltaDepth                 = evDelta.depth;
  deltaTT                    = evDelta.time;
  return true;
}

bool Solver::getObservationParamsChanges(unsigned evId,
                                         const std::string &staId,
                                         char phase,
                                         unsigned &startingTTObs,
                                         unsigned &startingCCObs,
                                         unsigned &finalTotalObs,
                                         double &meanAPrioriWeight,
                                         double &meanFinalWeight,
                                         double &meanObsResiduals,
                                         std::set<unsigned> &evIds) const
{
  unsigned evIdx;
  if (!_eventIdConverter.hasId(evId, evIdx)) return false;

  string phStaId = string(1, phase) + "@" + staId;
  unsigned phStaIdx;
  if (!_phStaIdConverter.hasId(phStaId, phStaIdx)) return false;

  const auto &it1 = _paramStats.find(evIdx);
  if (it1 == _paramStats.end()) return false;

  const auto &it2 = it1->second.find(phStaIdx);
  if (it2 == it1->second.end()) return false;

  const ParamStats &prmSts = it2->second;

  startingTTObs     = prmSts.startingTTObs;
  startingCCObs     = prmSts.startingCCObs;
  finalTotalObs     = prmSts.finalTotalObs;
  meanAPrioriWeight = 0;
  meanFinalWeight   = 0;
  meanObsResiduals  = 0;
  if ((startingTTObs + startingCCObs) > 0)
  {
    meanAPrioriWeight =
        prmSts.totalAPrioriWeight / (startingTTObs + startingCCObs);
  }
  if (finalTotalObs > 0)
  {
    meanFinalWeight  = prmSts.totalFinalWeight / finalTotalObs;
    meanObsResiduals = prmSts.totalResiduals / finalTotalObs;
  }
  evIds.insert(prmSts.peerEvIds.begin(), prmSts.peerEvIds.end());
  return true;
}

void Solver::loadSolutions()
{
  auto computeEventDelta = [this](unsigned evIdx,
                                  EventDeltas &evDelta) -> bool {
    const unsigned evOffset = evIdx * 4;

    evDelta.longitude = _dd->m[evOffset + 0];
    evDelta.latitude  = _dd->m[evOffset + 1];
    evDelta.depth     = _dd->m[evOffset + 2];
    evDelta.time      = _dd->m[evOffset + 3];

    if (!std::isfinite(evDelta.longitude) || !std::isfinite(evDelta.latitude) ||
        !std::isfinite(evDelta.depth) || !std::isfinite(evDelta.time))
    {
      return false;
    }

    return true;
  };

  //
  // Build a map of events that has at least one non-zero-weight observation
  // (i.e. discard events that lost all their observations due to
  // downweighting).
  //
  for (const auto &kv1 : _paramStats)
  {
    unsigned evIdx = kv1.first;
    bool allZero   = true;

    for (const auto &kv2 : kv1.second)
    {
      const ParamStats &pweight = kv2.second;
      if (pweight.totalFinalWeight > 0)
      {
        allZero = false;
        break;
      }
    }

    if (!allZero) _eventDeltas[evIdx] = {0};
  }

  //
  // Load change in event parameters for all events that have at least
  // one non-zero-weight observation.
  //
  for (auto it = _eventDeltas.begin(); it != _eventDeltas.end();)
  {
    const unsigned evIds = it->first;
    EventDeltas &evDelta = it->second;
    it                   = computeEventDelta(evIds, evDelta) ? std::next(it)
                                           : _eventDeltas.erase(it);
  }

  // free some memory
  _eventParams.clear();
  _residuals.clear();
  _dd = nullptr;
}

void Solver::computePartialDerivatives()
{
  //
  // Convert all events and station coordinates to an euclidean space in X,Y,Z
  // cartesian system centered around the cluster centroid, which is an
  // approximation of the original non-euclidean system (The approximation is
  // valid if the area covered by all events is small such that the Earth
  // curvature can be assumed flat.).
  // Next, compute the partial derivatives of the travel times with respect to
  // the new system.
  //
  _centroid = {0, 0, 0};
  vector<double> longitudes;
  for (const auto &kw : _eventParams)
  {
    _centroid.lat += kw.second.lat;
    _centroid.depth += kw.second.depth;
    longitudes.push_back(degToRad(kw.second.lon));
  }
  _centroid.lat /= _eventParams.size();
  _centroid.depth /= _eventParams.size();
  _centroid.lon = normalizeLon(radToDeg(computeCircularMean(longitudes)));

  auto convertCoord = [this](double lat, double lon, double depth, double &x,
                             double &y, double &z) {
    double distance, az;
    distance = computeDistance(this->_centroid.lat, this->_centroid.lon, 0, lat,
                               lon, 0, &az);
    az       = degToRad(az);
    x        = distance * std::sin(az);
    y        = distance * std::cos(az);
    z        = depth - _centroid.depth;
  };

  // convert events' coordinates
  for (auto &kv : _eventParams)
  {
    EventParams &evprm = kv.second;
    convertCoord(evprm.lat, evprm.lon, evprm.depth, evprm.x, evprm.y, evprm.z);
  }

  // convert stations' coordinates
  for (auto &kv : _stationParams)
  {
    StationParams &staprm = kv.second;
    convertCoord(staprm.lat, staprm.lon, -staprm.elevation / 1000., staprm.x,
                 staprm.y, staprm.z);
  }

  // compute derivatives
  for (auto &kv1 : _obsParams)
  {
    for (auto &kv2 : kv1.second)
    {
      ObservationParams &obprm = kv2.second;

      // dip angle:  0(down):180(up) -> -90(down):+90(up)
      const double dip = obprm.takeOffAngleDip - degToRad(90);
      // azimuth angle to backazimuth
      const double azi      = obprm.takeOffAngleAzim - degToRad(180);
      const double slowness = 1. / obprm.velocityAtSrc;

      obprm.dx = slowness * std::cos(dip) * std::sin(azi);
      obprm.dy = slowness * std::cos(dip) * std::cos(azi);
      obprm.dz = slowness * std::sin(dip);
    }
  }
}

multimap<double, unsigned> Solver::computeInterEventDistance() const
{
  if (_observations.size() < 1)
  {
    return multimap<double, unsigned>();
  }

  map<string, double> distCache;
  multimap<double, unsigned> dists;

  for (const auto &kw : _observations)
  {
    unsigned obIdx           = kw.first;
    const Observation &obsrv = kw.second;

    double interEvDistance;

    string key = obsrv.ev1Idx < obsrv.ev2Idx
                     ? to_string(obsrv.ev1Idx) + "-" + to_string(obsrv.ev2Idx)
                     : to_string(obsrv.ev2Idx) + "-" + to_string(obsrv.ev1Idx);

    auto it = distCache.find(key);
    if (it != distCache.end())
    {
      interEvDistance = it->second;
    }
    else
    {
      const EventParams &ev1Prm = _eventParams.at(obsrv.ev1Idx);
      const EventParams &ev2Prm = _eventParams.at(obsrv.ev2Idx);
      interEvDistance = computeDistance(ev1Prm.lat, ev1Prm.lon, ev1Prm.depth,
                                        ev2Prm.lat, ev2Prm.lon, ev2Prm.depth);
    }

    dists.emplace(interEvDistance, obIdx);
  }

  return dists;
}

/*
 * From Waldhauser's:
 *
 * W = max^2 ( 0, 1 - ( res / (alpha*resMAD/0.67449) )^2 )
 *
 */
vector<double> Solver::computeResidualWeights(const vector<double> &residuals,
                                              const double alpha) const
{
  if (residuals.size() < 1)
  {
    return vector<double>();
  }

  //
  // Find the median absolute deviation of residuals (MAD).
  //
  const double median = computeMedian(residuals);
  const double MAD    = computeMedianAbsoluteDeviation(residuals, median);

  logInfo("Solver: num DD %lu residual median %.3f [msec] "
          "MedianAbsoluteDeviation %.3f [msec]",
          _observations.size(), median * 1000, MAD * 1000);

  //
  // compute weights
  //
  vector<double> weights(residuals.size());

  const double MAD_gn = 0.67449; // MAD for gaussian noise
  for (unsigned i = 0; i < residuals.size(); i++)
  {
    double weight = 1. - square(residuals.at(i) / (alpha * MAD / MAD_gn));
    weight        = std::max(weight, 0.);
    weight        = square(weight);

    weights[i] = weight;
  }

  return weights;
}

void Solver::prepareDDSystem(double ttConstraint,
                             double dampingFactor,
                             double residualDownWeight)
{
  computePartialDerivatives();

  // Count how many travel time constraints we need in the DD system.
  unsigned ttconstraintNum = 0;
  if (ttConstraint > 0)
  {
    for (const auto &kv1 : _obsParams)
      for (const auto &kv2 : kv1.second)
        if (kv2.second.computeEvChanges) ttconstraintNum++;
  }

  // allocate DD system memory
  _dd.reset(new DDSystem(_observations.size(), _eventIdConverter.size(),
                         _phStaIdConverter.size(), ttconstraintNum));

  // initialize `m` and `L2NScaler`
  std::fill_n(_dd->m, _dd->numColsG, 0);
  std::fill_n(_dd->L2NScaler, _dd->numColsG, 1.);

  // initialize `G`
  for (const auto &kv1 : _obsParams)
  {
    unsigned evIdx = kv1.first;
    for (const auto &kv2 : kv1.second)
    {
      unsigned phStaIdx                          = kv2.first;
      const ObservationParams &obprm             = kv2.second;
      _dd->G[evIdx * _dd->nPhStas + phStaIdx][0] = obprm.dx;
      _dd->G[evIdx * _dd->nPhStas + phStaIdx][1] = obprm.dy;
      _dd->G[evIdx * _dd->nPhStas + phStaIdx][2] = obprm.dz;
      _dd->G[evIdx * _dd->nPhStas + phStaIdx][3] = 1.; // travel time
    }
  }

  // initialize: `W`, `d`, `evByObsi`, `phStaByObs` (`m` is zero initialized)
  for (auto &kw : _observations)
  {
    unsigned obIdx                  = kw.first;
    Observation &ob                 = kw.second;
    const ObservationParams &obprm1 = _obsParams.at(ob.ev1Idx).at(ob.phStaIdx);
    const ObservationParams &obprm2 = _obsParams.at(ob.ev2Idx).at(ob.phStaIdx);

    _dd->W[obIdx]          = ob.aPrioriWeight;
    _dd->evByObs[0][obIdx] = obprm1.computeEvChanges ? ob.ev1Idx : -1;
    _dd->evByObs[1][obIdx] = obprm2.computeEvChanges ? ob.ev2Idx : -1;
    _dd->phStaByObs[obIdx] = ob.phStaIdx;
    // compute double difference
    _dd->d[obIdx] =
        ob.observedDiffTime - (obprm1.travelTime - obprm2.travelTime);

    // apply weights to d
    _dd->d[obIdx] *= _dd->W[obIdx];

    // (bookkeeping) Keep track of the weights of observation parameters.
    if (obprm1.computeEvChanges)
    {
      ParamStats &prmSts = _paramStats[ob.ev1Idx][ob.phStaIdx];
      if (ob.isXcorr)
        prmSts.startingCCObs++;
      else
        prmSts.startingTTObs++;
      prmSts.totalAPrioriWeight += _dd->W[obIdx];
      prmSts.peerEvIds.insert(_eventIdConverter.fromIdx(ob.ev2Idx));
    }

    if (obprm2.computeEvChanges)
    {
      ParamStats &prmSts = _paramStats[ob.ev2Idx][ob.phStaIdx];
      if (ob.isXcorr)
        prmSts.startingCCObs++;
      else
        prmSts.startingTTObs++;
      prmSts.totalAPrioriWeight += _dd->W[obIdx];
      prmSts.peerEvIds.insert(_eventIdConverter.fromIdx(ob.ev1Idx));
    }
  }

  // downweight observations by residuals
  _residuals = vector<double>(_dd->d, _dd->d + _dd->nObs);
  if (residualDownWeight > 0)
  {
    vector<double> resWeights =
        computeResidualWeights(_residuals, residualDownWeight);
    for (unsigned obIdx = 0; obIdx < _dd->nObs; obIdx++)
    {
      _dd->W[obIdx] *= resWeights[obIdx];
      _dd->d[obIdx] *= resWeights[obIdx];
    }
  }

  //
  // (bookkeeping) Compute final stats for each observation.
  //
  for (unsigned int obIdx = 0; obIdx < _dd->nObs; obIdx++)
  {
    double observationWeight = _dd->W[obIdx];

    if (observationWeight == 0.) continue;

    const unsigned phStaIdx =
        _dd->phStaByObs[obIdx]; // station for this observation

    const int evIdx1 = _dd->evByObs[0][obIdx]; // event 1 for this observation
    if (evIdx1 >= 0)
    {
      ParamStats &prmSts = _paramStats.at(evIdx1).at(phStaIdx);
      prmSts.finalTotalObs++;
      prmSts.totalFinalWeight += observationWeight;
      prmSts.totalResiduals += _residuals.at(obIdx);
    }

    const int evIdx2 = _dd->evByObs[1][obIdx]; // event 2 for this observation
    if (evIdx2 >= 0)
    {
      ParamStats &prmSts = _paramStats.at(evIdx2).at(phStaIdx);
      prmSts.finalTotalObs++;
      prmSts.totalFinalWeight += observationWeight;
      prmSts.totalResiduals += _residuals.at(obIdx);
    }
  }

  // add travel time residual constraints after DD observations
  if (ttConstraint > 0)
  {
    unsigned ttconstraintIdx = _dd->nObs - 1;
    for (const auto &kv1 : _obsParams)
    {
      unsigned evIdx            = kv1.first;
      const auto &evObsParamMap = kv1.second;
      for (const auto &kv2 : evObsParamMap)
      {
        unsigned phStaIdx              = kv2.first;
        const ObservationParams &obprm = kv2.second;

        if (!obprm.computeEvChanges) continue;

        const auto &it1 = _paramStats.find(evIdx);
        if (it1 == _paramStats.end()) continue;

        const auto &it2 = it1->second.find(phStaIdx);
        if (it2 == it1->second.end()) continue;

        const ParamStats &prmSts = it2->second;

        if (++ttconstraintIdx >= _dd->numRowsG)
        {
          string msg = strf("Solver: internal logic error "
                            "(ttconstraintIdx=%u but _dd->numRowsG=%u)",
                            ttconstraintIdx, _dd->numRowsG);
          throw Exception(msg);
        }

        _dd->W[ttconstraintIdx] =
            ttConstraint *
            (prmSts.finalTotalObs != 0
                 ? (prmSts.totalFinalWeight / prmSts.finalTotalObs)
                 : 0);
        _dd->d[ttconstraintIdx] =
            -obprm.travelTimeResidual * _dd->W[ttconstraintIdx];
        _dd->evByObs[0][ttconstraintIdx] = evIdx;
        _dd->evByObs[1][ttconstraintIdx] = -1;
        _dd->phStaByObs[ttconstraintIdx] = phStaIdx;
      }
    }

    // In case _obsParams contains more entries than required by
    // _observations we fill the remaining entries with defaults
    while (++ttconstraintIdx < _dd->numRowsG)
    {
      _dd->W[ttconstraintIdx]          = 0;
      _dd->d[ttconstraintIdx]          = 0;
      _dd->evByObs[0][ttconstraintIdx] = -1;
      _dd->evByObs[1][ttconstraintIdx] = -1;
    }

    // just a safety belt
    if (ttconstraintNum != (ttconstraintIdx - _dd->nObs))
    {
      string msg = strf("Solver: internal logic error (ttconstraintNum=%u "
                        "but only added %u constraints)",
                        ttconstraintNum, (ttconstraintIdx + 1 - _dd->nObs));
      throw Exception(msg);
    }
  }

  // print residual by inter-event distance
  multimap<double, unsigned> obByDist = computeInterEventDistance();
  auto obByDistIt                     = obByDist.begin();
  while (obByDistIt != obByDist.end())
  {
    unsigned decileSize = (obByDist.size() / 10) + 1;
    vector<double> decileRes;
    decileRes.reserve(decileSize);

    double startingDist = obByDistIt->first;
    double finalDist    = obByDistIt->first;
    while (obByDistIt != obByDist.end() && decileRes.size() < decileSize)
    {
      unsigned obIdx = obByDistIt->second;
      decileRes.push_back(_residuals.at(obIdx));
      finalDist = obByDistIt->first;
      obByDistIt++;
    }

    const double median = computeMedian(decileRes);
    const double MAD    = computeMedianAbsoluteDeviation(decileRes, median);

    logInfo("Solver: Inter-event dist %.3f-%-.3f [km] num DD %lu residual "
            "median %6.3f [msec] MedianAbsoluteDeviation %6.3f [msec]",
            startingDist, finalDist, decileRes.size(), median * 1000,
            MAD * 1000);
  }

  // free some memory
  _observations.clear();
  _obsParams.clear();
  _stationParams.clear();
}

void Solver::solve(unsigned numIterations,
                   double ttConstraint,
                   double dampingFactor,
                   double residualDownWeight,
                   bool normalizeG)
{
  if (_observations.size() == 0)
  {
    throw Exception("Solver: no observations given");
  }

  if (_type == "LSQR")
  {
    _solve<lsqrBase>(numIterations, ttConstraint, dampingFactor,
                     residualDownWeight, normalizeG);
  }
  else if (_type == "LSMR")
  {
    _solve<lsmrBase>(numIterations, ttConstraint, dampingFactor,
                     residualDownWeight, normalizeG);
  }
  else
  {
    throw Exception(
        "Solver: invalid type, only LSQR and LSMR are valid methods");
  }
}

template <typename T>
void Solver::_solve(unsigned numIterations,
                    double ttConstraint,
                    double dampingFactor,
                    double residualDownWeight,
                    bool normalizeG)
{
  prepareDDSystem(ttConstraint, dampingFactor, residualDownWeight);

  Adapter<T> solver(*_dd); // keeps only a reference to _dd, doesn't copy it!!!
  if (normalizeG)
  {
    solver.L2normalize();
  }
  solver.SetDamp(dampingFactor);
  solver.SetMaximumNumberOfIterations(numIterations ? numIterations
                                                    : _dd->numColsG / 2);
  const double eps = std::numeric_limits<double>::epsilon();
  solver.SetEpsilon(eps);
  solver.SetToleranceA(1e-6); // we use [km] and [sec] in the DD system, so
  solver.SetToleranceB(1e-6); // this tolerance looks like enough (mm and usec)
  solver.SetUpperLimitOnConditional(1.0 / (10 * sqrt(eps)));

  std::ostringstream solverLogs;
  solver.SetOutputStream(solverLogs);

  solver.Solve(_dd->numRowsG, _dd->numColsG, _dd->d, _dd->m);

  logDebug("%s", solverLogs.str().c_str());

  logInfo("Stopped because %u : %s (used %u iterations)",
          solver.GetStoppingReason(), solver.GetStoppingReasonMessage().c_str(),
          solver.GetNumberOfIterationsPerformed());

  if (solver.GetStoppingReason() == 4)
  {
    _dd        = nullptr;
    string msg = strf("Solver: no solution found (%s)",
                      solver.GetStoppingReasonMessage().c_str());
    throw Exception(msg);
  }

  if (normalizeG)
  {
    solver.L2DeNormalize();
  }

  loadSolutions();

  if (_eventDeltas.empty())
  {
    throw Exception("Solver: no event has been relocated");
  }
}

} // namespace HDD

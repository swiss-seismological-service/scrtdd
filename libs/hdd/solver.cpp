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

#include "solver.h"
#include "log.h"
#include "lsmr.h"
#include "lsqr.h"
#include "utils.h"

#include <sstream>

using namespace std;
using namespace HDD::Logger;
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
  Adapter(HDD::DDSystem &dd) : _dd(dd)
  {
    std::fill_n(_dd.m, _dd.numColsG, 0);
    std::fill_n(_dd.L2NScaler, _dd.numColsG, 1.);

    for (unsigned int ob = 0; ob < _dd.numRowsG; ob++)
    {
      if (_dd.W[ob] != 0)
      {
        _dd.d[ob] *= _dd.W[ob];
      }
      else
      {
        _dd.d[ob] = 0;
      }
    }
  }

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
  void Aprod1(unsigned int m,
              unsigned int n,
              const double *x,
              double *y) const override
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
  void Aprod2(unsigned int m,
              unsigned int n,
              double *x,
              const double *y) const override
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
                            const std::string &phase,
                            double timeDiff,
                            double aPrioriWeight,
                            bool xcorrUsed,
                            double xcorrCoeff)
{
  string phStaId    = phase + "@" + staId;
  string obsId      = to_string(evId1) + "+" + to_string(evId2) + "_" + phStaId;
  unsigned evIdx1   = _eventIdConverter.convert(evId1);
  unsigned evIdx2   = _eventIdConverter.convert(evId2);
  unsigned phStaIdx = _phStaIdConverter.convert(phStaId);
  unsigned obsIdx   = _obsIdConverter.convert(obsId);
  _observations.insert(
      {obsIdx, Observation{evIdx1, evIdx2, phStaIdx, timeDiff, aPrioriWeight,
                           xcorrUsed, xcorrCoeff}});
}

void Solver::addEvent(unsigned evId, double evLat, double evLon, double evDepth)
{
  int evIdx = _eventIdConverter.convert(evId);
  _eventParams.insert({evIdx, EventParams{evLat, evLon, evDepth}});
}

void Solver::addStation(const std::string &staId,
                        double staLat,
                        double staLon,
                        double staElevation)
{
  unsigned staIdx = _staIdConverter.convert(staId);
  _stationParams.insert({staIdx, StationParams{staLat, staLon, staElevation}});
}

void Solver::addObservationParams(unsigned evId,
                                  const std::string &staId,
                                  const std::string &phase,
                                  bool computeEvChanges,
                                  double travelTime,
                                  double travelTimeResidual,
                                  double dx,
                                  double dy,
                                  double dz)
{
  string phStaId    = phase + "@" + staId;
  int evIdx         = _eventIdConverter.convert(evId);
  unsigned phStaIdx = _phStaIdConverter.convert(phStaId);
  _obsParams[evIdx].insert(
      {phStaIdx, ObservationParams{computeEvChanges, travelTime,
                                   travelTimeResidual, dx, dy, dz}});
}

bool Solver::getEvent(unsigned evId,
                      double &evLat,
                      double &evLon,
                      double &evDepth) const
{
  unsigned evIdx;
  if (!_eventIdConverter.hasId(evId, evIdx)) return false;

  const auto &it = _eventParams.find(evIdx);
  if (it == _eventParams.end()) return false;

  const EventParams &eprms = it->second;

  evLat   = eprms.lat;
  evLon   = eprms.lon;
  evDepth = eprms.depth;

  return true;
}

bool Solver::getStation(const std::string &staId,
                        double &staLat,
                        double &staLon,
                        double &staElevation) const
{
  unsigned staIdx;
  if (!_staIdConverter.hasId(staId, staIdx)) return false;

  const auto &it = _stationParams.find(staIdx);
  if (it == _stationParams.end()) return false;

  const StationParams &sprms = it->second;

  staLat       = sprms.lat;
  staLon       = sprms.lon;
  staElevation = sprms.elevation;

  return true;
}

bool Solver::getObservationParams(unsigned evId,
                                  const std::string &staId,
                                  const std::string &phase,
                                  bool &computeEvChanges,
                                  double &travelTime,
                                  double &travelTimeResidual,
                                  double &dx,
                                  double &dy,
                                  double &dz) const
{
  unsigned evIdx;
  if (!_eventIdConverter.hasId(evId, evIdx)) return false;

  string phStaId = phase + "@" + staId;
  unsigned phStaIdx;
  if (!_phStaIdConverter.hasId(phStaId, phStaIdx)) return false;

  const auto &it1 = _obsParams.find(evIdx);
  if (it1 == _obsParams.end()) return false;

  const auto &it2 = it1->second.find(phStaIdx);
  if (it2 == it1->second.end()) return false;

  const ObservationParams &oprms = it2->second;

  computeEvChanges   = oprms.computeEvChanges;
  travelTime         = oprms.travelTime;
  travelTimeResidual = oprms.travelTimeResidual;
  dx                 = oprms.dx;
  dy                 = oprms.dy;
  dz                 = oprms.dz;

  return true;
}

std::vector<Solver::DoubleDifference> Solver::getDoubleDifferences() const
{
  std::vector<DoubleDifference> results;
  results.reserve(_observations.size());
  for (auto &kw : _observations)
  {
    // unsigned obIdx        = kw.first;
    const Observation &ob = kw.second;

    // Get staId and phase
    string phStaId = _phStaIdConverter.fromIdx(ob.phStaIdx);
    size_t pos     = phStaId.find('@');
    if (pos == std::string::npos)
    {
      throw Exception("Solver: internal logic error (phStaId formmating)");
    }
    string phase = phStaId.substr(0, pos);
    string staId = phStaId.substr(pos + 1);

    DoubleDifference entry{_eventIdConverter.fromIdx(ob.ev1Idx),
                           _eventIdConverter.fromIdx(ob.ev2Idx),
                           staId,
                           phase,
                           ob.residualDownWeight,
                           ob.xcorrUsed,
                           ob.xcorrCoeff,
                           ob.timeDiff,
                           ob.computedTimeDiff,
                           ob.doubleDifference,
                           ob.interEventDistance};
    results.emplace_back(std::move(entry));
  }
  return results;
}

bool Solver::getEventChanges(unsigned evId,
                             double &deltaX,
                             double &deltaY,
                             double &deltaZ,
                             double &deltaTime) const
{
  unsigned evIdx;
  if (!_eventIdConverter.hasId(evId, evIdx)) return false;

  if (_eventDeltas.find(evIdx) == _eventDeltas.end()) return false;

  const EventDeltas &evDelta = _eventDeltas.at(evIdx);
  deltaX                     = evDelta.x;
  deltaY                     = evDelta.y;
  deltaZ                     = evDelta.z;
  deltaTime                  = evDelta.time;
  return true;
}

bool Solver::isEventPhaseUsed(unsigned evId,
                              const std::string &staId,
                              const std::string &phase) const
{
  unsigned evIdx;
  if (!_eventIdConverter.hasId(evId, evIdx)) return false;

  string phStaId = phase + "@" + staId;
  unsigned phStaIdx;
  if (!_phStaIdConverter.hasId(phStaId, phStaIdx)) return false;

  const auto &it1 = _stats.find(evIdx);
  if (it1 == _stats.end()) return false;

  const auto &it2 = it1->second.find(phStaIdx);
  if (it2 == it1->second.end()) return false;

  const Stats &stat = it2->second;

  return stat.finalTotalObs > 0;
}

void Solver::loadSolutions()
{
  auto fetchEventDelta = [this](unsigned evIdx, EventDeltas &evDelta) -> bool {
    const unsigned evOffset = evIdx * 4;

    evDelta.x    = _dd.m[evOffset + 0];
    evDelta.y    = _dd.m[evOffset + 1];
    evDelta.z    = _dd.m[evOffset + 2];
    evDelta.time = _dd.m[evOffset + 3];

    if (!std::isfinite(evDelta.x) || !std::isfinite(evDelta.y) ||
        !std::isfinite(evDelta.z) || !std::isfinite(evDelta.time))
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
  for (const auto &kv1 : _stats)
  {
    unsigned evIdx        = kv1.first;
    bool someObservations = false;

    for (const auto &kv2 : kv1.second)
    {
      const Stats &stat = kv2.second;
      if (stat.finalTotalObs > 0)
      {
        someObservations = true;
        break;
      }
    }

    if (someObservations) _eventDeltas[evIdx] = {};
  }

  //
  // Load change in event parameters for all events that have at least
  // one non-zero-weight observation.
  //
  for (auto it = _eventDeltas.begin(); it != _eventDeltas.end();)
  {
    const unsigned evIds = it->first;
    EventDeltas &evDelta = it->second;
    it                   = fetchEventDelta(evIds, evDelta) ? std::next(it)
                                                           : _eventDeltas.erase(it);
  }

  // free some memory
  _dd = DDSystem();
}

multimap<double, unsigned> Solver::computeInterEventDistance()
{
  if (_observations.size() < 1)
  {
    return multimap<double, unsigned>();
  }

  map<string, double> distCache;
  multimap<double, unsigned> dists;

  for (auto &kw : _observations)
  {
    unsigned obIdx     = kw.first;
    Observation &obsrv = kw.second;

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
    obsrv.interEventDistance = interEvDistance;
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

  logInfoF(
      "Total DD %zu, dd-residual [msec] median %.3f MedianAbsoluteDeviation "
      "%.3f",
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

void Solver::prepare(double ttConstraint, double residualDownWeight)
{
  //
  // Count how many travel time constraints we need in the DD system.
  //
  unsigned ttconstraintNum = 0;
  if (ttConstraint > 0)
  {
    for (const auto &kv1 : _obsParams)
      for (const auto &kv2 : kv1.second)
        if (kv2.second.computeEvChanges) ttconstraintNum++;
  }

  //
  // allocate DD system memory
  //
  _dd = DDSystem(_observations.size(), _eventIdConverter.size(),
                 _phStaIdConverter.size(), ttconstraintNum);

  //
  // initialize `G`
  //
  for (const auto &kv1 : _obsParams)
  {
    unsigned evIdx = kv1.first;
    for (const auto &kv2 : kv1.second)
    {
      unsigned phStaIdx                        = kv2.first;
      const ObservationParams &obprm           = kv2.second;
      _dd.G[evIdx * _dd.nPhStas + phStaIdx][0] = obprm.dx;
      _dd.G[evIdx * _dd.nPhStas + phStaIdx][1] = obprm.dy;
      _dd.G[evIdx * _dd.nPhStas + phStaIdx][2] = obprm.dz;
      _dd.G[evIdx * _dd.nPhStas + phStaIdx][3] = 1.; // travel time
    }
  }

  //
  // initialize: `W`, `d`, `evByObs`, `phStaByObs`
  //
  for (auto &kw : _observations)
  {
    unsigned obIdx                  = kw.first;
    Observation &ob                 = kw.second;
    const ObservationParams &obprm1 = _obsParams.at(ob.ev1Idx).at(ob.phStaIdx);
    const ObservationParams &obprm2 = _obsParams.at(ob.ev2Idx).at(ob.phStaIdx);

    // compute double difference
    ob.computedTimeDiff = obprm1.travelTime - obprm2.travelTime;
    ob.doubleDifference = ob.timeDiff - ob.computedTimeDiff;

    _dd.W[obIdx]          = ob.aPrioriWeight;
    _dd.evByObs[0][obIdx] = obprm1.computeEvChanges ? ob.ev1Idx : -1;
    _dd.evByObs[1][obIdx] = obprm2.computeEvChanges ? ob.ev2Idx : -1;
    _dd.phStaByObs[obIdx] = ob.phStaIdx;
    _dd.d[obIdx]          = ob.doubleDifference;
  }

  //
  // downweight observations by residuals
  //
  if (residualDownWeight > 0)
  {
    vector<double> residuals = vector<double>(_dd.d, _dd.d + _dd.nObs);
    vector<double> resWeights =
        computeResidualWeights(residuals, residualDownWeight);
    for (unsigned obIdx = 0; obIdx < _dd.nObs; obIdx++)
    {
      _dd.W[obIdx] *= resWeights[obIdx];
      _observations.at(obIdx).residualDownWeight = _dd.W[obIdx];
    }
  }

  //
  // (bookkeeping) Compute final stats for each observation.
  //
  for (unsigned int obIdx = 0; obIdx < _dd.nObs; obIdx++)
  {
    double observationWeight = _dd.W[obIdx];

    if (observationWeight == 0.) continue;

    const unsigned phStaIdx =
        _dd.phStaByObs[obIdx]; // station for this observation

    const int evIdx1 = _dd.evByObs[0][obIdx]; // event 1 for this observation
    if (evIdx1 >= 0)
    {
      Stats &stat = _stats[evIdx1][phStaIdx];
      stat.finalTotalObs++;
      stat.totalFinalWeight += observationWeight;
    }

    const int evIdx2 = _dd.evByObs[1][obIdx]; // event 2 for this observation
    if (evIdx2 >= 0)
    {
      Stats &stat = _stats[evIdx2][phStaIdx];
      stat.finalTotalObs++;
      stat.totalFinalWeight += observationWeight;
    }
  }

  //
  // add travel time residual constraints after DD observations
  //
  if (ttConstraint > 0)
  {
    unsigned ttconstraintIdx = _dd.nObs - 1;
    for (const auto &kv1 : _obsParams)
    {
      unsigned evIdx            = kv1.first;
      const auto &evObsParamMap = kv1.second;
      for (const auto &kv2 : evObsParamMap)
      {
        unsigned phStaIdx              = kv2.first;
        const ObservationParams &obprm = kv2.second;

        if (!obprm.computeEvChanges) continue;

        const auto &it1 = _stats.find(evIdx);
        if (it1 == _stats.end()) continue;

        const auto &it2 = it1->second.find(phStaIdx);
        if (it2 == it1->second.end()) continue;

        const Stats &stat = it2->second;

        if (++ttconstraintIdx >= _dd.numRowsG)
        {
          string msg = strf("Solver: internal logic error "
                            "(ttconstraintIdx=%u but _dd.numRowsG=%u)",
                            ttconstraintIdx, _dd.numRowsG);
          throw Exception(msg);
        }

        _dd.W[ttconstraintIdx] =
            ttConstraint * (stat.finalTotalObs > 0
                                ? (stat.totalFinalWeight / stat.finalTotalObs)
                                : 0);
        _dd.d[ttconstraintIdx]          = -obprm.travelTimeResidual;
        _dd.evByObs[0][ttconstraintIdx] = evIdx;
        _dd.evByObs[1][ttconstraintIdx] = -1;
        _dd.phStaByObs[ttconstraintIdx] = phStaIdx;
      }
    }

    // In case _obsParams contains more entries than required by
    // _observations we fill the remaining entries with defaults
    while (++ttconstraintIdx < _dd.numRowsG)
    {
      _dd.W[ttconstraintIdx]          = 0;
      _dd.d[ttconstraintIdx]          = 0;
      _dd.evByObs[0][ttconstraintIdx] = -1;
      _dd.evByObs[1][ttconstraintIdx] = -1;
    }

    // just a safety belt
    if (ttconstraintNum != (ttconstraintIdx - _dd.nObs))
    {
      string msg = strf("Solver: internal logic error (ttconstraintNum=%u "
                        "but only added %u constraints)",
                        ttconstraintNum, (ttconstraintIdx + 1 - _dd.nObs));
      throw Exception(msg);
    }
  }

  //
  // print residual by inter-event distance
  //
  multimap<double, unsigned> obByDist = computeInterEventDistance();
  auto obByDistIt                     = obByDist.begin();
  vector<string> xcorrLines;
  while (obByDistIt != obByDist.end())
  {
    unsigned decileSize = (obByDist.size() / 10) + 1;
    vector<double> decileRes, decileCoeff;
    decileRes.reserve(decileSize);
    decileCoeff.reserve(decileSize);

    double startingDist = obByDistIt->first;
    double finalDist    = obByDistIt->first;
    while (obByDistIt != obByDist.end() && decileRes.size() < decileSize)
    {
      unsigned obIdx = obByDistIt->second;
      decileRes.push_back(_dd.d[obIdx]);
      if (_observations.at(obIdx).xcorrCoeff > 0)
      {
        decileCoeff.push_back(_observations.at(obIdx).xcorrCoeff);
      }
      finalDist = obByDistIt->first;
      obByDistIt++;
    }

    double min, max, q1, q2, q3;
    compute5numberSummary(decileRes, min, max, q1, q2, q3);
    logInfoF("Event dist %.4f-%-.4f [km], num DD %zu dd-residual "
             "[msec] 1st quartile %6.3f median %6.3f 3rd quartile %6.3f",
             startingDist, finalDist, decileRes.size(), q1 * 1000, q2 * 1000,
             q3 * 1000);
    // when xcorr is not used there no entries
    if (decileCoeff.size() > 0)
    {
      compute5numberSummary(decileCoeff, min, max, q1, q2, q3);
      string line = strf(
          "Event dist %.4f-%-.4f [km], num CC %zu corr-coeff "
          "min %.2f 1st quartile %.2f median %.2f 3rd quartile %.2f max %.2f",
          startingDist, finalDist, decileCoeff.size(), min, q1, q2, q3, max);
      xcorrLines.push_back(std::move(line));
    }
  }

  // Print Xcorr 5 number summary
  for (const string &line : xcorrLines)
  {
    logInfo(line);
  }
}

void Solver::solve(unsigned numIterations,
                   double dampingFactor,
                   bool normalizeG)
{
  if (_observations.size() == 0)
  {
    throw Exception("Solver: no observations given");
  }

  // free some memory
  _observations.clear();
  _eventParams.clear();
  _stationParams.clear();
  _obsParams.clear();

  if (_type == "LSQR")
  {
    _solve<lsqrBase>(numIterations, dampingFactor, normalizeG,
                     std::set<unsigned>{4});
  }
  else if (_type == "LSMR")
  {
    _solve<lsmrBase>(numIterations, dampingFactor, normalizeG,
                     std::set<unsigned>{3, 6});
  }
  else
  {
    throw Exception(
        "Solver: invalid type, only LSQR and LSMR are valid methods");
  }
}

template <typename T>
void Solver::_solve(unsigned numIterations,
                    double dampingFactor,
                    bool normalizeG,
                    std::set<unsigned> rejectStoppingReasons)
{
  Adapter<T> solver(_dd); // keeps only a reference to _dd, doesn't copy it!!!
  if (normalizeG)
  {
    solver.L2normalize();
  }
  solver.SetDamp(dampingFactor);
  solver.SetMaximumNumberOfIterations(numIterations ? numIterations
                                                    : _dd.numColsG / 2);
  const double eps = std::numeric_limits<double>::epsilon();
  solver.SetEpsilon(eps);
  solver.SetToleranceA(1e-6); // we use [km] and [sec] in the DD system, so
  solver.SetToleranceB(1e-6); // this tolerance looks like enough (mm and usec)
  solver.SetUpperLimitOnConditional(1.0 / (10 * sqrt(eps)));

  std::ostringstream solverLogs;
  solver.SetOutputStream(solverLogs);

  solver.Solve(_dd.numRowsG, _dd.numColsG, _dd.d, _dd.m);

  logDebugF("%s", solverLogs.str().c_str());

  logInfoF("Stopped because %u : %s (used %u iterations)",
           solver.GetStoppingReason(),
           solver.GetStoppingReasonMessage().c_str(),
           solver.GetNumberOfIterationsPerformed());

  if (rejectStoppingReasons.count(solver.GetStoppingReason()) != 0)
  {
    _dd        = DDSystem(); // free memory
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

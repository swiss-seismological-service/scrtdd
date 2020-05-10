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

#include "solver.h"
#include "utils.h"

#include <stdexcept>
#include <sstream>
#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/math.h>
#include <seiscomp3/core/strings.h>

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp3/logging/log.h>

using namespace std;
using Seiscomp::Core::stringify;

namespace {

/**
 * Common DDSystem adapter for both LSQR and LSMR solvers
 * T can be lsqrBase or lsmrBase
 */
template <class T>
class Adapter : public T
{

public:

    Adapter() { }
    virtual ~Adapter() { }

    void setDDSytem(const Seiscomp::HDD::DDSystemPtr& dd)
    {
        _dd = dd;
    }

    /*
     * Scale G by normalizing the L2-norm of each column as suggested
     * by LSQR and LSMR solvers
     */
    void L2normalize()
    {
        std::fill_n(_dd->L2NScaler, _dd->numColsG, 0.);

        for ( unsigned int ob = 0; ob < _dd->nObs; ob++ )
        {
            const double obsW = _dd->W[ob];
            if ( obsW == 0. )
                continue;

            const unsigned phStaIdx = _dd->phStaByObs[ob]; // station for this observation

            const int evIdx1 = _dd->evByObs[ob][0]; // event 1 for this observation
            if ( evIdx1 >= 0 )
            {
                const unsigned idxG = evIdx1 * _dd->nPhStas + phStaIdx;
                const unsigned evOffset = evIdx1 * 4;
                _dd->L2NScaler[evOffset+0] += std::pow(_dd->G[idxG][0] * obsW, 2);
                _dd->L2NScaler[evOffset+1] += std::pow(_dd->G[idxG][1] * obsW, 2);
                _dd->L2NScaler[evOffset+2] += std::pow(_dd->G[idxG][2] * obsW, 2);
                _dd->L2NScaler[evOffset+3] += std::pow(_dd->G[idxG][3] * obsW, 2);
            }

            const int evIdx2 = _dd->evByObs[ob][1]; // event 2 for this observation
            if ( evIdx2 >= 0 )
            {
                const unsigned idxG = evIdx2 * _dd->nPhStas + phStaIdx;
                const unsigned evOffset = evIdx2 * 4;
                _dd->L2NScaler[evOffset+0] += std::pow(_dd->G[idxG][0] * obsW, 2);
                _dd->L2NScaler[evOffset+1] += std::pow(_dd->G[idxG][1] * obsW, 2);
                _dd->L2NScaler[evOffset+2] += std::pow(_dd->G[idxG][2] * obsW, 2);
                _dd->L2NScaler[evOffset+3] += std::pow(_dd->G[idxG][3] * obsW, 2);
            }
        }

        double const* meanShiftWeight = &_dd->W[_dd->nObs];
        if ( meanShiftWeight[0] != 0 || meanShiftWeight[1] != 0 ||
             meanShiftWeight[2] != 0 || meanShiftWeight[3] != 0 )
        {
            for (unsigned evOffset = 0; evOffset < _dd->numColsG; evOffset += 4 )
            {
                _dd->L2NScaler[evOffset+0] += std::pow(meanShiftWeight[0], 2);
                _dd->L2NScaler[evOffset+1] += std::pow(meanShiftWeight[1], 2);
                _dd->L2NScaler[evOffset+2] += std::pow(meanShiftWeight[2], 2);
                _dd->L2NScaler[evOffset+3] += std::pow(meanShiftWeight[3], 2); 
            }
        }

        for ( unsigned col = 0; col < _dd->numColsG; col++)
        {
            _dd->L2NScaler[col] = 1. / std::sqrt(_dd->L2NScaler[col]);
        }
    }

    /*
     * Rescale m back to the initial scaling
     */
    void L2DeNormalize()
    {
        for (unsigned evOffset = 0; evOffset < _dd->numColsG; evOffset += 4 )
        {
            _dd->m[evOffset+0] *= _dd->L2NScaler[evOffset+0];
            _dd->m[evOffset+1] *= _dd->L2NScaler[evOffset+1];
            _dd->m[evOffset+2] *= _dd->L2NScaler[evOffset+2];
            _dd->m[evOffset+3] *= _dd->L2NScaler[evOffset+3];
        }
    }

    /**
     * Required by lsqrBase and lsmrBase:
     *
     * computes y = y + A*x without altering x,
     * where A is a matrix of dimensions A[m][n].
     * The size of the vector x is n.
     * The size of the vector y is m.
     */
    void Aprod1(unsigned int m, unsigned int n, const double * x, double * y ) const
    {
        if ( m != _dd->numRowsG || n != _dd->numColsG )
        {
            string msg = stringify("Solver: Internal logic error (m=%u n=%u but G=%ux%u)",
                                   m, n, _dd->numRowsG, _dd->numColsG);
            throw std::runtime_error(msg.c_str());
        }

        for ( unsigned int ob = 0; ob < _dd->nObs; ob++ )
        {
            if ( _dd->W[ob] == 0. )
                continue;

            const unsigned phStaIdx = _dd->phStaByObs[ob]; // station for this observation
            double sum = 0;

            const int evIdx1 = _dd->evByObs[ob][0]; // event 1 for this observation
            if ( evIdx1 >= 0 )
            {
                const unsigned idxG = evIdx1 * _dd->nPhStas + phStaIdx;
                const unsigned evOffset = evIdx1 * 4;
                sum += _dd->G[idxG][0] * _dd->L2NScaler[evOffset+0] * x[evOffset+0];
                sum += _dd->G[idxG][1] * _dd->L2NScaler[evOffset+1] * x[evOffset+1];
                sum += _dd->G[idxG][2] * _dd->L2NScaler[evOffset+2] * x[evOffset+2];
                sum += _dd->G[idxG][3] * _dd->L2NScaler[evOffset+3] * x[evOffset+3];
            }

            const int evIdx2 = _dd->evByObs[ob][1]; // event 2 for this observation
            if ( evIdx2 >= 0 )
            {
                const unsigned idxG = evIdx2 * _dd->nPhStas + phStaIdx;
                const unsigned evOffset = evIdx2 * 4;
                sum -= _dd->G[idxG][0] * _dd->L2NScaler[evOffset+0] * x[evOffset+0];
                sum -= _dd->G[idxG][1] * _dd->L2NScaler[evOffset+1] * x[evOffset+1];
                sum -= _dd->G[idxG][2] * _dd->L2NScaler[evOffset+2] * x[evOffset+2];
                sum -= _dd->G[idxG][3] * _dd->L2NScaler[evOffset+3] * x[evOffset+3];
            }

            y[ob] += _dd->W[ob] * sum;
        }

        double *meanShiftWeight = &_dd->W[_dd->nObs];
        if ( meanShiftWeight[0] != 0 || meanShiftWeight[1] != 0 ||
             meanShiftWeight[2] != 0 || meanShiftWeight[3] != 0 )
        {
            double meanShift[4] = {0};
            for (unsigned evOffset = 0; evOffset < _dd->numColsG; evOffset += 4 )
            {
                meanShift[0] += x[evOffset+0] * _dd->L2NScaler[evOffset+0];
                meanShift[1] += x[evOffset+1] * _dd->L2NScaler[evOffset+1];
                meanShift[2] += x[evOffset+2] * _dd->L2NScaler[evOffset+2];
                meanShift[3] += x[evOffset+3] * _dd->L2NScaler[evOffset+3]; 
            }
            y[_dd->nObs+0] += meanShift[0] * meanShiftWeight[0];
            y[_dd->nObs+1] += meanShift[1] * meanShiftWeight[1];
            y[_dd->nObs+2] += meanShift[2] * meanShiftWeight[2];
            y[_dd->nObs+3] += meanShift[3] * meanShiftWeight[3];
        }
    }

    /**
     * Required by lsqrBase and lsmrBase:
     *
     * computes x = x + A'*y without altering y,
     * where A is a matrix of dimensions A[m][n].
     * The size of the vector x is n.
     * The size of the vector y is m.
     */
    void Aprod2(unsigned int m, unsigned int n, double * x, const double * y ) const
    {
        if ( m != _dd->numRowsG || n != _dd->numColsG )
        {
            string msg = stringify("Solver: Internal logic error (m=%u n=%u but G=%ux%u)",
                                   m, n, _dd->numRowsG, _dd->numColsG);
            throw std::runtime_error(msg.c_str());
        }

        for ( unsigned int ob = 0; ob < _dd->nObs; ob++ )
        {
            const double wY = y[ob] * _dd->W[ob];
            if ( wY == 0. )
                continue;

            const unsigned phStaIdx = _dd->phStaByObs[ob]; // station for this observation

            const int evIdx1 = _dd->evByObs[ob][0]; // event 1 for this observation
            if ( evIdx1 >= 0 )
            {
                const unsigned idxG = evIdx1 * _dd->nPhStas + phStaIdx;
                const unsigned evOffset = evIdx1 * 4;
                x[evOffset+0] += _dd->G[idxG][0] * _dd->L2NScaler[evOffset+0] * wY;
                x[evOffset+1] += _dd->G[idxG][1] * _dd->L2NScaler[evOffset+1] * wY;
                x[evOffset+2] += _dd->G[idxG][2] * _dd->L2NScaler[evOffset+2] * wY;
                x[evOffset+3] += _dd->G[idxG][3] * _dd->L2NScaler[evOffset+3] * wY;
            }

            const int evIdx2 = _dd->evByObs[ob][1]; // event 2 for this observation
            if ( evIdx2 >= 0 )
            {
                const unsigned idxG = evIdx2 * _dd->nPhStas + phStaIdx;
                const unsigned evOffset = evIdx2 * 4;
                x[evOffset+0] -= _dd->G[idxG][0] * _dd->L2NScaler[evOffset+0] * wY;
                x[evOffset+1] -= _dd->G[idxG][1] * _dd->L2NScaler[evOffset+1] * wY;
                x[evOffset+2] -= _dd->G[idxG][2] * _dd->L2NScaler[evOffset+2] * wY;
                x[evOffset+3] -= _dd->G[idxG][3] * _dd->L2NScaler[evOffset+3] * wY;
            }
        }

        double *meanShiftWeight = &_dd->W[_dd->nObs];
        if ( meanShiftWeight[0] != 0 || meanShiftWeight[1] != 0 ||
             meanShiftWeight[2] != 0 || meanShiftWeight[3] != 0 )
        {
            for (unsigned evOffset = 0; evOffset < _dd->numColsG; evOffset += 4 )
            {
                x[evOffset+0] += meanShiftWeight[0] * y[_dd->nObs+0] * _dd->L2NScaler[evOffset+0];
                x[evOffset+1] += meanShiftWeight[1] * y[_dd->nObs+1] * _dd->L2NScaler[evOffset+1];
                x[evOffset+2] += meanShiftWeight[2] * y[_dd->nObs+2] * _dd->L2NScaler[evOffset+2];
                x[evOffset+3] += meanShiftWeight[3] * y[_dd->nObs+3] * _dd->L2NScaler[evOffset+3];
            }
        }
    }

private:

    Seiscomp::HDD::DDSystemPtr _dd;
};


}


namespace Seiscomp {
namespace HDD {


void
Solver::addObservation(unsigned evId1, unsigned evId2, const std::string& staId, char phase,
                       double observedDiffTime, double aPrioriWeight,
                       bool computeEv1Changes, bool computeEv2Changes, bool isXcorr)
{
    string phStaId = string(1,phase) + "@" + staId;
    string obsId = to_string(evId1) + "+" + to_string(evId2) + "_" + phStaId;
    unsigned evIdx1 = _eventIdConverter.convert(evId1);
    unsigned evIdx2 = _eventIdConverter.convert(evId2);
    unsigned phStaIdx = _phStaIdConverter.convert(phStaId);
    unsigned obsIdx = _obsIdConverter.convert(obsId);
    _observations[obsIdx] = Observation( {evIdx1, evIdx2, phStaIdx, computeEv1Changes,
            computeEv2Changes, observedDiffTime, aPrioriWeight, isXcorr});
}


void
Solver::addObservationParams(unsigned evId, const std::string& staId, char phase,
                             double evLat, double evLon, double evDepth,
                             double staLat, double staLon, double staElevation,
                             double travelTime)
{
    string phStaId = string(1,phase) + "@" + staId;
    int evIdx  = _eventIdConverter.convert(evId);
    unsigned phStaIdx = _phStaIdConverter.convert(phStaId);
    _eventParams[evIdx] = EventParams( {evLat, evLon, evDepth, 0, 0, 0} );
    _stationParams[phStaIdx] = StationParams( {staLat, staLon, staElevation, 0, 0, 0} );
    _obsParams[evIdx][phStaIdx] = ObservationParams( {travelTime, 0, 0, 0, 0} );
}


bool
Solver::getEventChanges(unsigned evId, double &deltaLat, double &deltaLon, double &deltaDepth, double &deltaTT) const
{
    if ( ! _eventIdConverter.hasId(evId) )
        return false;

    unsigned evIdx = _eventIdConverter.toIdx(evId);

    if ( _eventDeltas.find(evIdx) == _eventDeltas.end() )
        return false;

    const EventDeltas& evDelta = _eventDeltas.at(evIdx);
    deltaLat = evDelta.deltaLat;
    deltaLon = evDelta.deltaLon;
    deltaDepth = evDelta.deltaDepth;
    deltaTT = evDelta.deltaTT;
    return true;
}


bool
Solver::getObservationParamsChanges(unsigned evId, const std::string& staId, char phase,
                                    unsigned &startingObservations, 
                                    unsigned &startingXcorrObservations,
                                    unsigned &totalFinalObservations, 
                                    double &meanAPrioriWeight,
                                    double &meanFinalWeight) const
{
    if ( ! _eventIdConverter.hasId(evId) )
        return false;

    string phStaId = string(1,phase) + "@" + staId;
    if ( ! _phStaIdConverter.hasId(phStaId) )
        return false;

    int evIdx  = _eventIdConverter.toIdx(evId);
    unsigned phStaIdx = _phStaIdConverter.toIdx(phStaId);

    const auto& it1 = _paramStats.find(evIdx);
    if ( it1 == _paramStats.end() )
        return false;

    const auto& it2 = it1->second.find(phStaIdx);
    if ( it2 == it1->second.end() )
        return false;

    const ParamStats& prmSts = it2->second;

    startingObservations = prmSts.startingObservations;
    startingXcorrObservations = prmSts.startingXcorrObservations;
    totalFinalObservations    = prmSts.totalFinalObservations;
    meanAPrioriWeight = 0;
    if ( (startingObservations + startingXcorrObservations) > 0 )
        meanAPrioriWeight = prmSts.totalAPrioriWeight / (startingObservations + startingXcorrObservations);
    meanFinalWeight = totalFinalObservations ? (prmSts.totalFinalWeight / totalFinalObservations) : 0;
    return true;
}


void
Solver::loadSolutions()
{
    auto computeEventDelta = [this](unsigned evIdx, EventDeltas& evDelta)
    {
        const EventParams& evprm = _eventParams.at(evIdx);
        const unsigned evOffset = evIdx * 4;

        double deltaX      = _dd->m[evOffset+0];
        double deltaY      = _dd->m[evOffset+1];
        evDelta.deltaDepth = _dd->m[evOffset+2];
        evDelta.deltaTT    = _dd->m[evOffset+3];

        double newX = evprm.x + deltaX;
        double newY = evprm.y + deltaY;

        // compute distance and azimuth of evId to centroid (0,0,0)
        double hdist = std::sqrt( std::pow(newX,2) + std::pow(newY, 2) );
        hdist = Math::Geo::km2deg(hdist); // distance to degree

        double azimuth  = std::atan2(newX, newY);
        azimuth = rad2deg(azimuth);

        // Computes the coordinates (lat, lon) of the point which
        // is at a degree azimuth of 'azi' and a distance of 'dist' as seen
        // from the centroid (lat0, lon0)
        double newLat, newLon;
        Math::Geo::delandaz2coord(hdist, azimuth, _centroid.lat, _centroid.lon, &newLat, &newLon);

        evDelta.deltaLat = newLat - evprm.lat;
        evDelta.deltaLon = newLon - evprm.lon;
    };

    //
    // Compute final weghts for each ObservationParams
    // Note: we could have done this even before solving the system
    //       but here is more convenient because we might eventually
    //       add more information depending on the solution
    //
    for ( unsigned int ob = 0; ob < _dd->nObs; ob++ )
    {
        double observationWeight = _dd->W[ob];

        if ( observationWeight == 0. ) continue;

        const unsigned phStaIdx = _dd->phStaByObs[ob]; // station for this observation

        const int evIdx1 = _dd->evByObs[ob][0]; // event 1 for this observation
        if ( evIdx1 >= 0 )
        {
            ParamStats& prmSts = _paramStats.at(evIdx1).at(phStaIdx);
            prmSts.totalFinalObservations++;
            prmSts.totalFinalWeight += observationWeight;
        }

        const int evIdx2 = _dd->evByObs[ob][1]; // event 2 for this observation
        if ( evIdx2 >= 0 )
        {
            ParamStats& prmSts = _paramStats.at(evIdx2).at(phStaIdx);
            prmSts.totalFinalObservations++;
            prmSts.totalFinalWeight += observationWeight;
        }
    }

    //
    // Now build a map of events that have at least one observation whose weight is non zero
    // (i.e. discard events that lost all their observations due to downweighting )
    //
    for ( const auto& kv1 : _paramStats )
    {
        unsigned evIdx = kv1.first;
        bool allZero = true;

        for ( const auto& kv2 : kv1.second )
        {
            const ParamStats& pweight = kv2.second;
            if ( pweight.totalFinalWeight > 0 )
            {
                allZero = false;
                break;
            }
        }

        if ( ! allZero ) _eventDeltas[evIdx] = { 0 };
    }

    //
    // Load change in event parameters for all events that have at least
    // one non-zero-weight observation
    //
    for ( auto& kw : _eventDeltas )
    {
        const unsigned evIds = kw.first;
        EventDeltas& evDelta = kw.second;
        computeEventDelta( evIds, evDelta );
    }

    // free some memory
    _eventParams.clear();
    _dd = nullptr;
}


void
Solver::computePartialDerivatives()
{
    //
    // Fist convert all events and staions coordinates to an euclidean space
    // in X,Y,Z cartesian system centered around the cluster centroid, which
    // is an approximation of the original non-euclidean system.
    // The approximation works when the area covered by all events is small
    // and the Earth curvature can be assumed flat
    // Then compute the partial derivatives of the travel times with respect to
    // the new system
    //
    _centroid = {0,0,0};
    for ( const auto& kw: _eventParams )
    {
        _centroid.lat   += kw.second.lat;
        _centroid.lon   += kw.second.lon;
        _centroid.depth += kw.second.depth;
    }
    _centroid.lat   /= _eventParams.size();
    _centroid.lon   /= _eventParams.size();
    _centroid.depth /= _eventParams.size();

    auto convertCoord = [this](double lat, double lon, double depth,
                           double& x, double& y, double& z)
    {
        double distance, az;
        distance = computeDistance(this->_centroid.lat, this->_centroid.lon, 0, lat, lon, 0, &az);
        az = deg2rad(az);
        x = distance * std::sin(az);
        y = distance * std::cos(az);
        z = depth - _centroid.depth;
    };

    //
    // convert events coordinates
    //
    for ( auto& kv : _eventParams )
    {
        EventParams& evprm = kv.second;
        convertCoord(evprm.lat, evprm.lon, evprm.depth, evprm.x,  evprm.y, evprm.z);
    }

    // convert stations coordinates
    for ( auto& kv : _stationParams )
    {
        StationParams& staprm = kv.second;
        convertCoord(staprm.lat, staprm.lon, -staprm.elevation/1000., staprm.x,  staprm.y, staprm.z);
    }

    //
    // compute derivatives
    //
    for ( auto& kv1 : _obsParams )
    {
        unsigned evIdx = kv1.first;

        for ( auto& kv2 : kv1.second )
        {
            unsigned phStaIdx = kv2.first;
            ObservationParams& obsprm = kv2.second;
            const EventParams& evprm = _eventParams.at(evIdx);
            const StationParams& staprm = _stationParams.at(phStaIdx);

            double distance = computeDistance(evprm.lat, evprm.lon, evprm.depth,
                                               staprm.lat, staprm.lon, -staprm.elevation/1000.);

            double angle   = std::atan2( evprm.y - staprm.y, evprm.x - staprm.x);
            double takeOff = std::atan2( evprm.z - staprm.z, evprm.x - staprm.x);

            obsprm.slowness = obsprm.travelTime / distance;
            obsprm.dx = obsprm.slowness * std::cos(angle);
            obsprm.dy = obsprm.slowness * std::sin(angle);
            obsprm.dz = obsprm.slowness * std::sin(takeOff);
        }
    }
}


multimap<double,unsigned>
Solver::computeInterEventDistance() const
{
    if ( _observations.size() < 1 )
    {
        return multimap<double,unsigned>();
    }

    map<string, double> distCache;
    multimap<double,unsigned> dists;

    for ( const auto& kw: _observations )
    {
        unsigned obIdx = kw.first;
        const Observation& obsrv = kw.second;

        double interEvDistance;

        string key =  obsrv.ev1Idx < obsrv.ev2Idx
                   ? to_string(obsrv.ev1Idx) + "-" + to_string(obsrv.ev2Idx)
                   : to_string(obsrv.ev2Idx) + "-" + to_string(obsrv.ev1Idx);

        auto it = distCache.find(key);
        if ( it != distCache.end() )
        {
            interEvDistance = it->second;
        }
        else
        {
            const EventParams& ev1Prm = _eventParams.at(obsrv.ev1Idx);
            const EventParams& ev2Prm = _eventParams.at(obsrv.ev2Idx);
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
vector<double>
Solver::computeResidualWeights(const vector<double>& residuals, const double alpha) const
{
    if ( residuals.size() < 1 )
    {
        return vector<double>();
    }

    //
    // Find the median absolute deviation of residuals (MAD)
    // 
    const double median = computeMedian(residuals);
    const double MAD = computeMedianAbsoluteDeviation(residuals, median);

    SEISCOMP_INFO("Solver: #observations %lu residual median %.1f [msec] MedianAbsoluteDeviation %.1f [msec]",
                  _observations.size(), median*1000, MAD*1000);

    //
    // Compute weights
    //
    vector<double> weights( residuals.size() );
 
    const double MAD_gn = 0.67449; // MAD for gaussian noise
    for ( unsigned i = 0; i < residuals.size(); i++ )
    {
        double weight = 1. - std::pow( residuals.at(i)/(alpha*MAD/MAD_gn), 2);
        weight = std::max(weight, 0.);
        weight = std::pow(weight, 2);

        weights[i] = weight;
    }

    return weights;
}


void
Solver::prepareDDSystem(array<double,4> meanShiftConstraint, double residualDownWeight)
{
    computePartialDerivatives();

    _dd = DDSystemPtr(new DDSystem(_observations.size(),  _eventIdConverter.size(), _phStaIdConverter.size()) );

    // Init m and L2NScaler
    std::fill_n(_dd->m, _dd->numColsG, 0);
    std::fill_n(_dd->L2NScaler, _dd->numColsG, 1.);

    // initialize G
    for ( const auto& kv1 : _obsParams )
    {
        unsigned evIdx = kv1.first;
        for ( const auto& kv2 : kv1.second )
        {
            unsigned phStaIdx = kv2.first;
            const ObservationParams& obsprm = kv2.second;
            _dd->G[evIdx * _dd->nPhStas + phStaIdx][0] = obsprm.dx;
            _dd->G[evIdx * _dd->nPhStas + phStaIdx][1] = obsprm.dy;
            _dd->G[evIdx * _dd->nPhStas + phStaIdx][2] = obsprm.dz;
            _dd->G[evIdx * _dd->nPhStas + phStaIdx][3] = 1.; // travel time
        }
    }

    // initialize: W, d, evByObsi, phStaByObs
    // note: m is zero initialized
    for ( auto& kw: _observations )
    {
        unsigned obIdx = kw.first;
        Observation& obsrv = kw.second;

        _dd->W[obIdx] = obsrv.aPrioriWeight;
        _dd->evByObs[obIdx][0] = obsrv.computeEv1Changes ? obsrv.ev1Idx : -1;
        _dd->evByObs[obIdx][1] = obsrv.computeEv2Changes ? obsrv.ev2Idx : -1;
        _dd->phStaByObs[obIdx] = obsrv.phStaIdx;

        // compute double difference
        const ObservationParams& obsprm1 = _obsParams.at(obsrv.ev1Idx).at(obsrv.phStaIdx);
        const ObservationParams& obsprm2 = _obsParams.at(obsrv.ev2Idx).at(obsrv.phStaIdx); 
        _dd->d[obIdx] = obsrv.observedDiffTime - (obsprm1.travelTime - obsprm2.travelTime);

        // apply weights to d
        _dd->d[obIdx] *= _dd->W[obIdx];

        // keep track of the wights for these obsparms
        if ( obsrv.computeEv1Changes )
        {
            ParamStats& prmSts = _paramStats[obsrv.ev1Idx][obsrv.phStaIdx];
            if ( obsrv.isXcorr ) prmSts.startingXcorrObservations++;
            else prmSts.startingObservations++;
            prmSts.totalAPrioriWeight += _dd->W[obIdx];
        }

        if ( obsrv.computeEv2Changes )
        {
            ParamStats& prmSts = _paramStats[obsrv.ev2Idx][obsrv.phStaIdx];
            if ( obsrv.isXcorr ) prmSts.startingXcorrObservations++;
            else prmSts.startingObservations++;
            prmSts.totalAPrioriWeight += _dd->W[obIdx]; 
        }
    }

    // Init remaining 4 equations for cluster zero mean shift and their weights
    _dd->d[_dd->nObs+0] = 0;
    _dd->d[_dd->nObs+1] = 0;
    _dd->d[_dd->nObs+2] = 0;
    _dd->d[_dd->nObs+3] = 0;
    _dd->W[_dd->nObs+0] = meanShiftConstraint[0];
    _dd->W[_dd->nObs+1] = meanShiftConstraint[1];
    _dd->W[_dd->nObs+2] = meanShiftConstraint[2];
    _dd->W[_dd->nObs+3] = meanShiftConstraint[3];

    // downweight observations by residuals
    vector<double> residuals(_dd->d, _dd->d+_dd->nObs);
    if ( residualDownWeight > 0 )
    {
        vector<double> resWeights = computeResidualWeights(residuals, residualDownWeight);
        for ( unsigned obIdx = 0; obIdx < _dd->nObs; obIdx++ )
        {
            _dd->W[obIdx] *= resWeights[obIdx];
            _dd->d[obIdx] *= resWeights[obIdx]; 
        }
    }

    // Print residual by inter-event distance information
    multimap<double,unsigned> obByDist = computeInterEventDistance();
    auto obByDistIt = obByDist.begin();
    while ( obByDistIt != obByDist.end() )
    {
        unsigned decileSize = (obByDist.size() / 10) + 1;
        vector<double> decileRes;
        decileRes.reserve( decileSize );

        double startingDist = obByDistIt->first;
        double finalDist    = obByDistIt->first;
        while ( obByDistIt != obByDist.end() && decileRes.size() < decileSize )
        {
            unsigned obIdx = obByDistIt->second;
            decileRes.push_back( residuals.at(obIdx) );
            finalDist    = obByDistIt->first;
            obByDistIt++;
        }

        const double median = computeMedian(decileRes);
        const double MAD = computeMedianAbsoluteDeviation(decileRes, median);

        SEISCOMP_INFO("Solver: Inter-event dist %.2f-%-.2f [km] #observations %lu residual median %4.1f [msec] MedianAbsoluteDeviation %4.1f [msec]", 
                startingDist, finalDist, decileRes.size(), median*1000, MAD*1000); 

    }

    // free some memory 
    _observations.clear();
    _obsParams.clear();
    _stationParams.clear();
}


void 
Solver::solve(unsigned numIterations,
              double dampingFactor,
              double residualDownWeight,
              double meanLonShiftConstraint,
              double meanLatShiftConstraint,
              double meanDepthShiftConstraint,
              double meanTTShiftConstraint,
              bool normalizeG)
{
    if ( _observations.size() == 0 )
    {
        throw runtime_error("Solver: no observations given");
    }

    array<double,4> meanShiftConstraint = {
        meanLonShiftConstraint,
        meanLatShiftConstraint,
        meanDepthShiftConstraint,
        meanTTShiftConstraint,
    };

    if ( _type == "LSQR" )
    {
        _solve<lsqrBase>(numIterations, dampingFactor, residualDownWeight,
                         meanShiftConstraint, normalizeG);
    }
    else if ( _type == "LSMR" )
    {
        _solve<lsmrBase>(numIterations, dampingFactor, residualDownWeight,
                         meanShiftConstraint, normalizeG);
    }
    else
    {
        throw runtime_error("Solver: invalid type, only LSQR and LSMR are valid methods");
    }
}


template <class T>
void Solver::_solve(unsigned numIterations,
                    double dampingFactor,
                    double residualDownWeight,
                    array<double,4> meanShiftConstraint,
                    bool normalizeG)
{
    prepareDDSystem(meanShiftConstraint, residualDownWeight);

    Adapter<T> solver;
    solver.setDDSytem(_dd);
    if ( normalizeG )
    {
        solver.L2normalize();
    }

    solver.SetDamp(dampingFactor);
    solver.SetMaximumNumberOfIterations(numIterations ? numIterations : _dd->numColsG/2);

    const double eps = 1e-15;
    solver.SetEpsilon( eps );
    solver.SetToleranceA( 1e-16 );
    solver.SetToleranceB( 1e-16 );
    solver.SetUpperLimitOnConditional( 1.0 / ( 10 * sqrt( eps ) ) );

    //std::ostringstream solverLogs;
    //solver.SetOutputStream( solverLogs );

    solver.Solve(_dd->numRowsG, _dd->numColsG, _dd->d, _dd->m );

    //SEISCOMP_DEBUG("%s", solverLogs.str().c_str() );

    SEISCOMP_INFO("Stopped because %u : %s (used %u Iterations)", solver.GetStoppingReason(),
          solver.GetStoppingReasonMessage().c_str(), solver.GetNumberOfIterationsPerformed());

    if ( solver.GetStoppingReason() == 4 )
    {
        _dd = nullptr;
        string msg = stringify("Solver: no solution found (%s)", solver.GetStoppingReasonMessage().c_str() );
        throw runtime_error(msg.c_str());
    }

    if ( normalizeG )
    {
        solver.L2DeNormalize();
    }

    loadSolutions();

    if ( _eventDeltas.empty() )
    {
        throw runtime_error("Solver: no event has been relocated");
    } 
}


} // HDD
} // Seiscomp

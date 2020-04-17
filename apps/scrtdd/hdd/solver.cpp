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
 * Common adapter for both LSQR and LSMR solvers:
 *
 *      W*G*m = W*d -> A*x = b
 *
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
        if ( m != _dd->numRowsG )
            throw std::runtime_error("Solver: Internal logic error (m != numRowsG)");
        if ( n != (_dd->numColsG) )
            throw std::runtime_error("Solver: Internal logic error (n != numColsG");

        for ( unsigned int ob = 0; ob < _dd->nObs; ob++ )
        {
            const unsigned phStaIdx = _dd->phStaByObs[ob]; // station for this observation

            const int evIdx1 = _dd->evByObs[ob][0]; // event 1 for this observation
            if ( evIdx1 >= 0 )
            {
                y[ob] += _dd->G[evIdx1 * _dd->nPhStas + phStaIdx][0] * x[evIdx1 * 4 + 0];
                y[ob] += _dd->G[evIdx1 * _dd->nPhStas + phStaIdx][1] * x[evIdx1 * 4 + 1];
                y[ob] += _dd->G[evIdx1 * _dd->nPhStas + phStaIdx][2] * x[evIdx1 * 4 + 2];
                y[ob] += _dd->G[evIdx1 * _dd->nPhStas + phStaIdx][3] * x[evIdx1 * 4 + 3];
            }

            const int evIdx2 = _dd->evByObs[ob][1]; // event 2 for this observation
            if ( evIdx2 >= 0 )
            {
                y[ob] -= _dd->G[evIdx2 * _dd->nPhStas + phStaIdx][0] * x[evIdx2 * 4 + 0];
                y[ob] -= _dd->G[evIdx2 * _dd->nPhStas + phStaIdx][1] * x[evIdx2 * 4 + 1];
                y[ob] -= _dd->G[evIdx2 * _dd->nPhStas + phStaIdx][2] * x[evIdx2 * 4 + 2];
                y[ob] -= _dd->G[evIdx2 * _dd->nPhStas + phStaIdx][3] * x[evIdx2 * 4 + 3];
            } 
        }

        // do not waste computing if meanshift minimization is not used
        double *meanShiftWeight = &_dd->W[_dd->nObs];
        if ( meanShiftWeight[0] != 0 || meanShiftWeight[1] != 0 ||
             meanShiftWeight[2] != 0 || meanShiftWeight[3] != 0 )
        {
            double meanShift[4] = {0};
            for (unsigned evIdx = 0; evIdx < _dd->nEvts; evIdx++ )
            {
                meanShift[0] += x[evIdx * 4 + 0];
                meanShift[1] += x[evIdx * 4 + 1];
                meanShift[2] += x[evIdx * 4 + 2];
                meanShift[3] += x[evIdx * 4 + 3]; 
            }
            y[_dd->nObs + 0] += meanShift[0] * meanShiftWeight[0];
            y[_dd->nObs + 1] += meanShift[1] * meanShiftWeight[1];
            y[_dd->nObs + 2] += meanShift[2] * meanShiftWeight[2];
            y[_dd->nObs + 3] += meanShift[3] * meanShiftWeight[3];
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
        if ( m != _dd->numRowsG )
            throw std::runtime_error("Solver: Internal logic error (m != numRowsG)");
        if ( n != (_dd->numColsG) )
            throw std::runtime_error("Solver: Internal logic error (n != numColsG"); 

        for ( unsigned int ob = 0; ob < _dd->nObs; ob++ )
        {
            const unsigned phStaIdx = _dd->phStaByObs[ob]; // station for this observation

            const int evIdx1 = _dd->evByObs[ob][0]; // event 1 for this observation
            if ( evIdx1 >= 0 )
            {
                x[evIdx1 * 4 + 0] += _dd->G[evIdx1 * _dd->nPhStas + phStaIdx][0] * y[ob];
                x[evIdx1 * 4 + 1] += _dd->G[evIdx1 * _dd->nPhStas + phStaIdx][1] * y[ob];
                x[evIdx1 * 4 + 2] += _dd->G[evIdx1 * _dd->nPhStas + phStaIdx][2] * y[ob];
                x[evIdx1 * 4 + 3] += _dd->G[evIdx1 * _dd->nPhStas + phStaIdx][3] * y[ob];
            }

            const int evIdx2 = _dd->evByObs[ob][1]; // event 2 for this observation
            if ( evIdx2 >= 0 )
            {
                x[evIdx2 * 4 + 0] -= _dd->G[evIdx2 * _dd->nPhStas + phStaIdx][0] * y[ob];
                x[evIdx2 * 4 + 1] -= _dd->G[evIdx2 * _dd->nPhStas + phStaIdx][1] * y[ob];
                x[evIdx2 * 4 + 2] -= _dd->G[evIdx2 * _dd->nPhStas + phStaIdx][2] * y[ob];
                x[evIdx2 * 4 + 3] -= _dd->G[evIdx2 * _dd->nPhStas + phStaIdx][3] * y[ob];
            }
        }

        // do not waste computing if meanshift minimization is not used
        double *meanShiftWeight = &_dd->W[_dd->nObs];
        if ( meanShiftWeight[0] != 0 || meanShiftWeight[1] != 0 ||
             meanShiftWeight[2] != 0 || meanShiftWeight[3] != 0 )
        {
            for (unsigned evIdx = 0; evIdx < _dd->nEvts; evIdx++ )
            {
                x[evIdx * 4 + 0] += meanShiftWeight[0] * y[_dd->nObs + 0];
                x[evIdx * 4 + 1] += meanShiftWeight[1] * y[_dd->nObs + 1];
                x[evIdx * 4 + 2] += meanShiftWeight[2] * y[_dd->nObs + 2];
                x[evIdx * 4 + 3] += meanShiftWeight[3] * y[_dd->nObs + 3];
            }
        } 
    }

private:

    Seiscomp::HDD::DDSystemPtr _dd;
}; 


}


namespace Seiscomp {
namespace HDD {


/*
 * Compute distance in km between two points and optionally
 * azimuth and backazimuth
 */
double
Solver::computeDistance(double lat1, double lon1, double depth1,
                        double lat2, double lon2, double depth2,
                        double *azimuth, double *backAzimuth)
{
    double Hdist, az, baz;
    Math::Geo::delazi(lat1, lon1, lat2, lon2, &Hdist, &az, &baz);
    Hdist = Math::Geo::deg2km(Hdist);

    if (azimuth) *azimuth = az;
    if (backAzimuth) *backAzimuth = baz;

    if ( depth1 == depth2 )
        return Hdist;

    // this is an approximation that works when the distance is small
    // and the Earth curvature can be assumed flat 
    double Vdist = abs(depth1 - depth2);
    return std::sqrt( std::pow(Hdist,2) + std::pow(Vdist,2) );
}

 
Solver::Solver(std::string type)
    : _type(type)
{
}


void
Solver::addObservation(unsigned evId1, unsigned evId2, const std::string& staId, char phase,
                       double observedDiffTime, double weight, bool ev2Fixed)
{
    int evIdx1 = _eventIdConverter.convert(evId1);
    int evIdx2 = _eventIdConverter.convert(evId2);
    unsigned phStaIdx = _phStaIdConverter.convert(string(1,phase) + "@" + staId);
    _observations.push_back( Observation( {evIdx1, evIdx2, phStaIdx, observedDiffTime, weight, ev2Fixed} ) );
}


void
Solver::addObservationParams(unsigned evId, const std::string& staId, char phase,
                             double evLat, double evLon, double evDepth,
                             double staLat, double staLon, double staElevation,
                             double travelTime)
{
    int evIdx  = _eventIdConverter.convert(evId);
    unsigned phStaIdx = _phStaIdConverter.convert(string(1,phase) + "@" + staId);
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
    const EventParams& evprm = _eventParams.at(evIdx);

    double deltaX = _dd->m[evIdx*4];
    double deltaY = _dd->m[evIdx*4+1];
    deltaDepth    = -_dd->m[evIdx*4+2]; // make depth positive
    deltaTT       = _dd->m[evIdx*4+3];

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

    deltaLat = newLat - evprm.lat;
    deltaLon = newLon - evprm.lon;

    return true;
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
        z = _centroid.depth - depth;
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


void
Solver::prepareDDSystem(bool useObservationWeghts, double meanShiftWeight)
{
    computePartialDerivatives();

    _dd = DDSystemPtr(new DDSystem(_observations.size(),  _eventIdConverter.size(), _phStaIdConverter.size()) );

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
    int ob = 0;
    for ( const Observation& obsrv : _observations )
    {
        _dd->W[ob] = obsrv.weight;
        _dd->evByObs[ob][0] = obsrv.ev1Idx;
        _dd->evByObs[ob][1] = obsrv.ev2Fixed ? -1 : obsrv.ev2Idx;
        _dd->phStaByObs[ob] = obsrv.phStaIdx;

        // compute double difference
        const ObservationParams& obsprm1 = _obsParams.at(obsrv.ev1Idx).at(obsrv.phStaIdx);
        const ObservationParams& obsprm2 = _obsParams.at(obsrv.ev2Idx).at(obsrv.phStaIdx); 
        _dd->d[ob] = obsrv.observedDiffTime - (obsprm1.travelTime - obsprm2.travelTime);
        ob++;
    }

    // cluster zero mean shift and its weights
    _dd->d[_dd->nObs + 0] = 0;
    _dd->d[_dd->nObs + 1] = 0;
    _dd->d[_dd->nObs + 2] = 0;
    _dd->d[_dd->nObs + 3] = 0;
    _dd->W[_dd->nObs + 0] = meanShiftWeight;
    _dd->W[_dd->nObs + 1] = meanShiftWeight;
    _dd->W[_dd->nObs + 2] = meanShiftWeight;
    _dd->W[_dd->nObs + 3] = meanShiftWeight;

    // perform observation weighting
    if ( useObservationWeghts )
       _dd->applyWeights();

    // free memory 
    _observations.clear();
    _obsParams.clear();
    _stationParams.clear();
}


void Solver::solve(bool useObservationWeghts, double dampingFactor,
                   double meanShiftWeight, unsigned numIterations)
{
    if ( _type == "LSQR" )
    {
        _solve<lsqrBase>(useObservationWeghts, dampingFactor, meanShiftWeight, numIterations);
    }
    else if ( _type == "LSMR" )
    {
        _solve<lsmrBase>(useObservationWeghts, dampingFactor, meanShiftWeight, numIterations);
    }
    else
    {
        throw runtime_error("Solver: invalid type, only LSQR and LSMR are valid");
    }
}


template <class T>
void Solver::_solve(bool useObservationWeghts, double dampingFactor,
                    double meanShiftWeight, unsigned numIterations)
{
    prepareDDSystem(useObservationWeghts, meanShiftWeight);

    Adapter<T> solver;
    solver.setDDSytem(_dd);

    solver.SetDamp(dampingFactor);
    solver.SetMaximumNumberOfIterations(numIterations);

    const double eps = 1e-15;
    solver.SetEpsilon( eps );
    solver.SetToleranceA( 1e-16 );
    solver.SetToleranceB( 1e-16 );
    solver.SetUpperLimitOnConditional( 1.0 / ( 10 * sqrt( eps ) ) );

    std::ostringstream solverLogs;
    solver.SetOutputStream( solverLogs );

    solver.Solve(_dd->numRowsG, _dd->numColsG, _dd->d, _dd->m );

    SEISCOMP_INFO("%s", solverLogs.str().c_str() );
    SEISCOMP_INFO("Stopped because %u : %s", solver.GetStoppingReason(), solver.GetStoppingReasonMessage().c_str());
    SEISCOMP_INFO("Used %u Iterations", solver.GetNumberOfIterationsPerformed());
    SEISCOMP_INFO("Frobenius norm estimation of Abar = %.4f", solver.GetFrobeniusNormEstimateOfAbar());
    SEISCOMP_INFO("Condition number estimation of Abar = %.4f", solver.GetConditionNumberEstimateOfAbar());
    SEISCOMP_INFO("Estimate of final value of norm(rbar) = %.4f", solver.GetFinalEstimateOfNormRbar());
    SEISCOMP_INFO("Estimate of final value of norm of residuals = %.7f", solver.GetFinalEstimateOfNormOfResiduals());
    SEISCOMP_INFO("Estimate of norm of final solution = %.7f", solver.GetFinalEstimateOfNormOfX());
}


} // HDD
} // Seiscomp
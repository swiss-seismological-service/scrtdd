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
        if ( m != _dd->nObs ) throw std::runtime_error("Solver: Internal logic error (m != nObs)");
        if ( n != (_dd->nEvts*4) ) throw std::runtime_error("Solver: Internal logic error (n != nEvts)");

        for ( unsigned int ob = 0; ob < _dd->nObs; ob++ )
        {
            const unsigned phStaIdx = _dd->phStaByObs[ob]; // station for this observation
            double sum = 0.0;
            for ( unsigned int ev = 0; ev < 2; ev++ )
            {
                const int evIdx = _dd->evByObs[ob][ev]; // event for this observation

                if ( evIdx < 0 ) // this part of the observation is constant
                    continue;

                for ( unsigned int param = 0; param < 4; param++ )
                {
                    sum += _dd->G[evIdx * _dd->nPhStas + phStaIdx][param] * x[evIdx*4+param];
                }
            }
            y[ob] += sum;
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
        if ( m != _dd->nObs ) throw std::runtime_error("Solver: Internal logic error (m != nObs)");
        if ( n != (_dd->nEvts*4) ) throw std::runtime_error("Solver: Internal logic error (n != nEvts)");

        for ( unsigned int ob = 0; ob < _dd->nObs; ob++ )
        {
            const unsigned phStaIdx = _dd->phStaByObs[ob]; // station for this observation
            for ( unsigned int ev = 0; ev < 2; ev++ )
            {
                const int evIdx = _dd->evByObs[ob][ev]; //  event for this observation

                if ( evIdx < 0 ) // this part of the observation is constant
                    continue; 

                for ( unsigned int param = 0; param < 4; param++ )
                {
                    x[evIdx*4+param] += _dd->G[evIdx * _dd->nPhStas + phStaIdx][param] * y[ob];
                }
            }
        }
    }

private:

    Seiscomp::HDD::DDSystemPtr _dd;
}; 


}


namespace Seiscomp {
namespace HDD {


Solver::Solver(std::string type)
    : _type(type)
{
}


void
Solver::addObservation(unsigned evId1, unsigned evId2, const std::string& staId, char phase,
                       double doubleDiff, double weight)
{
    int evIdx1 = _eventIdConverter.convert(evId1);
    int evIdx2 = _eventIdConverter.convert(evId2);
    unsigned phStaIdx = _phStaIdConverter.convert(string(1,phase) + "@" + staId);
    _observations.push_back( Observation( {evIdx1, evIdx2, phStaIdx, doubleDiff, weight} ) );
}


void
Solver::addObservation(unsigned evId, const std::string& staId, char phase,
                       double doubleDiff, double weight)
{
    int evIdx = _eventIdConverter.convert(evId);
    unsigned phStaIdx = _phStaIdConverter.convert(string(1,phase) + "@" + staId);
    _observations.push_back( Observation( {evIdx, -1 , phStaIdx, doubleDiff, weight} ) );
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
    _partialDeriv[evIdx][phStaIdx] = PartialDerivatives( {travelTime, 0, 0, 0, 0} );
}


void
Solver::getEventChanges(unsigned evId, double &deltaLat, double &deltaLon, double &deltaDepth, double &deltaTT) const
{
    unsigned evIdx = _eventIdConverter.toIdx(evId);

    double deltaX = _dd->m[evIdx*4];
    double deltaY = _dd->m[evIdx*4+1];
    deltaDepth    = _dd->m[evIdx*4+2]; // no conversion required
    deltaTT       = _dd->m[evIdx*4+3]; // no conversion required

    double newX = _eventParams.at(evIdx).x + deltaX;
    double newY = _eventParams.at(evIdx).y + deltaY;

    // compute distance and azimuth of evId to centroid
    double distance = std::sqrt( std::pow(newX,2) + std::pow(newY, 2) );
    double azimuth  = std::atan2(newY, newX);

    // Computes the coordinates (lat, lon) of the point which
    // is at an azimuth of 'azi' and a distance of 'dist' as seen
    // from the centroid (lat0, lon0)
    Math::Geo::delandaz2coord(distance, azimuth, _centroid.lat, _centroid.lon,
                              &deltaLat, &deltaLon);
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
        double distance, az, baz;
        Math::Geo::delazi(lat, lon, this->_centroid.lat, this->_centroid.lon, &distance, &az, &baz);
        distance = Math::Geo::deg2km(distance); // distance to km
        az = deg2rad(az);
        x = distance * std::cos(az);
        y = distance * std::sin(az);
        z = depth - _centroid.depth;
    };

    for ( auto& kv1 : _partialDeriv )
    {
        unsigned evIdx = kv1.first;

        for ( auto& kv2 : kv1.second )
        {
            unsigned phStaIdx = kv2.first;
            PartialDerivatives& deriv = kv2.second;

            // convert event and station coordinates
            convertCoord(_eventParams.at(evIdx).lat, _eventParams.at(evIdx).lon,  _eventParams.at(evIdx).depth,
                         _eventParams.at(evIdx).x,  _eventParams.at(evIdx).y, _eventParams.at(evIdx).z);
            convertCoord(_stationParams.at(phStaIdx).lat, _stationParams.at(phStaIdx).lon,  -_stationParams.at(phStaIdx).elevation/1000.,
                         _stationParams.at(phStaIdx).x,  _stationParams.at(phStaIdx).y, _stationParams.at(phStaIdx).z);

            // compute derivatives
            double distance, az, baz;
            Math::Geo::delazi(_eventParams.at(evIdx).lat, _eventParams.at(evIdx).lon,
                              _stationParams.at(phStaIdx).lat, _stationParams.at(phStaIdx).lon,
                              &distance, &az, &baz);
            distance = Math::Geo::deg2km(distance); // distance to km
            deriv.slowness = deriv.travelTime / distance;

            double azimuth = std::atan2( _stationParams.at(phStaIdx).y - _eventParams.at(evIdx).y,
                                         _stationParams.at(phStaIdx).x - _eventParams.at(evIdx).x);
            double takeOff = std::atan2( _stationParams.at(phStaIdx).z - _eventParams.at(evIdx).z,
                                         _stationParams.at(phStaIdx).x - _eventParams.at(evIdx).x);

            deriv.dx = deriv.slowness * std::cos(azimuth);
            deriv.dy = deriv.slowness * std::sin(azimuth);
            deriv.dz = deriv.slowness * std::sin(takeOff);
        }
    }
}


void
Solver::prepareDDSystem()
{
    // compute partial derivatives
    computePartialDerivatives();

    _dd = DDSystemPtr(new DDSystem(_observations.size(),  _eventIdConverter.size(), _phStaIdConverter.size()) );

    // initialize G
    for ( const auto& kv1 : _partialDeriv )
    {
        unsigned evIdx = kv1.first;
        for ( const auto& kv2 : kv1.second )
        {
            unsigned phStaIdx = kv2.first;
            const PartialDerivatives& deriv = kv2.second;
            _dd->G[evIdx * _dd->nPhStas + phStaIdx][0] = deriv.dx;
            _dd->G[evIdx * _dd->nPhStas + phStaIdx][1] = deriv.dy;
            _dd->G[evIdx * _dd->nPhStas + phStaIdx][2] = deriv.dz;
            _dd->G[evIdx * _dd->nPhStas + phStaIdx][3] = 1.; // travel time
        }
    }

    // initialize: W, d, evByObsi, phStaByObs
    // note: m is zero initialized
    int ob = 0;
    for ( const Observation& obsrv : _observations )
    {
        // force throwing exception if no parameter was provided
        if ( obsrv.ev1Idx >= 0) _partialDeriv.at(obsrv.ev1Idx).at(obsrv.phStaIdx);
        if ( obsrv.ev2Idx >= 0) _partialDeriv.at(obsrv.ev2Idx).at(obsrv.phStaIdx);

        _dd->W[ob] = obsrv.weight;
        _dd->evByObs[ob][0] = obsrv.ev1Idx;
        _dd->evByObs[ob][1] = obsrv.ev2Idx;
        _dd->phStaByObs[ob] = obsrv.phStaIdx;
        _dd->d[ob] = obsrv.doubleDifference;
        ob++;
    }

    _dd->precomputeWeighting();

    // free memory 
    _observations.clear();
    _partialDeriv.clear();
    _stationParams.clear();
}


void Solver::solve()
{
    if ( _type == "LSQR" )
    {
        _solve<lsqrBase>();
    }
    else if ( _type == "LSMR" )
    {
        _solve<lsmrBase>();
    }
    else
    {
        throw runtime_error("Sovler: invalid type, only LSQR and LSMR are valid");
    }
}


template <class T>
void Solver::_solve()
{
    prepareDDSystem();

    Adapter<T> solver;
    solver.setDDSytem(_dd);

    const double eps = 1e-15;
    solver.SetEpsilon( eps );
    solver.SetDamp( 0.0 );
    solver.SetMaximumNumberOfIterations( 20 );
    solver.SetToleranceA( 1e-16 );
    solver.SetToleranceB( 1e-16 );
    solver.SetUpperLimitOnConditional( 1.0 / ( 10 * sqrt( eps ) ) );

    std::ostringstream solverLogs;
    solver.SetOutputStream( solverLogs );

//    double se[_dd->nEvts*4];
//    solver.SetStandardErrorEstimates( se );
//    solver.SetStandardErrorEstimatesFlag( true );

    solver.Solve(_dd->nObs, _dd->nEvts*4, _dd->d, _dd->m );

    SEISCOMP_INFO("%s", solverLogs.str().c_str() );
    SEISCOMP_INFO("Stopped because %u : %s", solver.GetStoppingReason(), solver.GetStoppingReasonMessage().c_str());
    SEISCOMP_INFO("Used %u Iterations", solver.GetNumberOfIterationsPerformed());
    SEISCOMP_INFO("Frobenius norm estimation of Abar = %.4f", solver.GetFrobeniusNormEstimateOfAbar());
    SEISCOMP_INFO("Condition number estimation of Abar = %.4f", solver.GetConditionNumberEstimateOfAbar());
    SEISCOMP_INFO("Estimate of final value of norm(rbar) = %.4f", solver.GetFinalEstimateOfNormRbar());
    SEISCOMP_INFO("Estimate of final value of norm of residuals = %.4f", solver.GetFinalEstimateOfNormOfResiduals());
    SEISCOMP_INFO("Estimate of norm of final solution = %.4f", solver.GetFinalEstimateOfNormOfX());
}


} // HDD
} // Seiscomp

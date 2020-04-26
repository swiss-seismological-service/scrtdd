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

#ifndef __RTDD_APPLICATIONS_SOLVER_H__
#define __RTDD_APPLICATIONS_SOLVER_H__

#include "lsqr.h"
#include "lsmr.h"

#include <seiscomp3/core/baseobject.h>
#include <unordered_map>
#include <vector>

namespace Seiscomp {
namespace HDD {


/*
 * Store data for a double-difference problem as described in Waldhaused &
 * Ellsworth 2000 paper:
 *
 *      W*G*m = d*W;
 *
 * Where G contains the partial derivatives of the travel times with respect to
 * event location and origin times.
 * m is a vector containing the changes in hypocentral parameters we wish to
 * determine for each event ( delta x, delta y, delta z and delta travel time)
 * d is the data vector containing the double-differences
 * W is a diagonal matrix to weight each equation.
 *
 * This class also contains 4 additional equations for constraining the mean
 * shift of all earthquakes during relocation.
 *
 * We take advantage of the sparsness of G matrix, so G is not a full matrix
 */
struct DDSystem : public Core::BaseObject {

    // number of observations
    const unsigned nObs;
    // number of events
    const unsigned nEvts;
    // number of stations
    const unsigned nPhStas;
    // W[nObs+4]: weight of each observation + cluster mean shift constraints (x,y,z,time) 
    double *W;
    // G[nEvts*nPhStas][4]: 3 partial derivatives for each event/station pair + tt (dx,dy,dz,1)
    double (*G)[4];
    // m[nEvts*4]: changes for each event hypocentral parameters we wish to determine (x,y,z,t)
    double (*m);
    // d[nObs+4]: double differences, one for each observation + mean shift constraints (x,y,z,time)
    double *d;
    // L2NScaler[nEvts*4]: L2 norm scaler for each G column
    double *L2NScaler;
    // evByObs[nObs][2]: map of 2 event idx for each observation (index -1 means no parameters)
    int (*evByObs)[2];
    // phStaByObs[nObs]: map of station idx for each observation
    unsigned *phStaByObs;

    const unsigned numRowsG;
    const unsigned numColsG;

    DDSystem(unsigned _nObs, unsigned _nEvts, unsigned _nPhStas)
        : nObs(_nObs), nEvts(_nEvts), nPhStas(_nPhStas), numRowsG(nObs+4), numColsG(nEvts*4)
    {
        W = new double[numRowsG];
        G = new double[nEvts*nPhStas][4];
        m = new double[numColsG];
        d = new double[numRowsG];
        L2NScaler = new double[numColsG];
        evByObs = new int[nObs][2];
        phStaByObs = new unsigned[nObs];
    }

    virtual ~DDSystem()
    {
        delete[] phStaByObs;
        delete[] evByObs;
        delete[] L2NScaler;
        delete[] d;
        delete[] m;
        delete[] G;
        delete[] W; 
    }

private:
    DDSystem( const DDSystem& other ) = delete;
    DDSystem operator=( const DDSystem& other ) = delete;
};

DEFINE_SMARTPOINTER(DDSystem);

/*
 * Solver for double difference problems.
 *
 * For details see Waldhauser & Ellsworth 2000 paper
 */
class Solver : public Core::BaseObject
{

public:
    Solver(std::string type) : _type(type) {}
    virtual ~Solver() {}

    void reset() { *this = Solver(_type); }

    void addObservation(unsigned evId1, unsigned evId2, const std::string& staId, char phase,
                        double diffTime, double aPrioriWeight,
                        bool computeEv1Changes, bool computeEv2Changes);

    void addObservationParams(unsigned evId, const std::string& staId, char phase,
                              double evLat, double evLon, double evDepth,
                              double staLat, double staLon, double staElevation,
                              double travelTime);

    void solve(unsigned numIterations=0, double dampingFactor=0,
               double meanShiftConstrainWeight=0, double residualDownWeight=0,
               bool normalizeG=true);

    bool getEventChanges(unsigned evId, double &deltaLat, double &deltaLon,
                         double &deltaDepth, double &deltaTT) const;

    static double computeDistance(double lat1, double lon1, double depth1,
                                  double lat2, double lon2, double depth2,
                                  double *azimuth = nullptr, double *backAzimuth = nullptr);
private:

    void computePartialDerivatives();

    std::vector<double> computeResidualWeights(std::vector<double> residuals, const double alpha);

    void prepareDDSystem(double meanShiftWeigh, double residualDownWeight);

    template <class T>
    void _solve(unsigned numIterations, double dampingFactor,
                double meanShiftConstrainWeight, double residualDownWeight,
                bool normalizeG );

private:

    /*
     *  Convert some hashable id of type T (e.g. string) to an alternative 
     *  rappresentaion as a sequentially growing integer starting from 0 
     *  (suitable for array index)
     */
    template <class T>
    class IdToIndex
    {
      public:
        unsigned convert(const T& id)
        {
            if ( _to.find(id) == _to.end() )
            {
                unsigned newIdx = _currentIdx++;
                _to[id] = newIdx;
                _from[newIdx] = id;
            }
            return _to.at(id);
        }

        unsigned toIdx(const T& id) const { return _to.at(id); }
        T fromIdx(unsigned idx) const { return _from.at(idx); }

        bool hasIdx(unsigned idx) const { return _from.find(idx) != _from.end(); } 
        bool hasId(const T& id) const { return _to.find(id) != _to.end(); } 

        unsigned size() { return _to.size(); }

      private:
        unsigned _currentIdx = 0;
        std::unordered_map<T,unsigned> _to;
        std::unordered_map<unsigned,T> _from;
    };
    IdToIndex<unsigned> _eventIdConverter;
    IdToIndex<std::string> _phStaIdConverter;
    IdToIndex<std::string> _obsIdConverter;

    struct Observation {
        unsigned ev1Idx;
        unsigned ev2Idx;
        unsigned phStaIdx;
        bool computeEv1Changes;
        bool computeEv2Changes;
        double observedDiffTime;
        double aPrioriWeight;
    };
    std::unordered_map<unsigned,Observation> _observations; // key = obsIdx

    struct EventParams {
        double lat, lon, depth;
        double x, y, z; // km
    };
    std::unordered_map<unsigned,EventParams> _eventParams; // key = evIdx

    struct StationParams {
        double lat, lon, elevation;
        double x, y, z; // km
    };
    std::unordered_map<unsigned,StationParams> _stationParams;  // key = phStaIdx

    struct ObservationParams {
        double travelTime;
        double slowness;
        double dx;
        double dy;
        double dz;
    };
    // key1=evIdx  key2=phStaIdx
    std::unordered_map<unsigned, std::unordered_map<unsigned,ObservationParams>> _obsParams;

    struct {
        double lat, lon, depth;
    } _centroid;

    DDSystemPtr _dd;
    std::string _type;
};

DEFINE_SMARTPOINTER(Solver);


}
}

#endif

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

#ifndef __RTDD_APPLICATIONS_XCORRCACHE_H__
#define __RTDD_APPLICATIONS_XCORRCACHE_H__

#include "catalog.h"

#include <unordered_map>

namespace Seiscomp {
namespace HDD {

class XCorrCache {

public:

    struct Entry {

        double mean_coeff;
        double mean_lag;
        double min_lag;
        double max_lag;
        unsigned ccCount;

        typedef struct {
            double coeff, lag, lowerUncertainty, upperUncertainty;
        } PeerInfo;
        std::unordered_map<unsigned,const PeerInfo> peers;

        std::string peersStr; // debug

        void update(const Catalog::Event& event, const Catalog::Phase& phase,
                    double coeff, double lag)
        {
            PeerInfo pi= {coeff, lag, phase.lowerUncertainty, phase.upperUncertainty};
            peers.insert( std::pair<unsigned,const PeerInfo>(event.id, pi) );
            peersStr   += std::string(event) + " ";
        }

        void computeStats()
        {
            ccCount = peers.size();
            mean_coeff = 0;
            mean_lag   = 0;
            min_lag    = 0;
            max_lag    = 0;
            for ( auto& pair : peers )
            {
                const PeerInfo& data = pair.second;
                mean_coeff += std::abs(data.coeff);
                mean_lag   += data.lag;
                min_lag    += data.lag - data.lowerUncertainty;
                max_lag    += data.lag + data.upperUncertainty; 
            }
            mean_coeff /= ccCount;
            mean_lag   /= ccCount;
            min_lag    /= ccCount;
            max_lag    /= ccCount; 
        }
    };

    Entry& getForUpdate(unsigned evId, const std::string& stationId,
                                  const Catalog::Phase::Type& type)
    {
        std::string key = make_key(evId, stationId, type);
        return resultsByPhase[key];
    }

    void remove(unsigned evId, const std::string& stationId,
                const Catalog::Phase::Type& type)
    {
        std::string key = make_key(evId, stationId, type);
        resultsByPhase.erase(key);
    }

    void computeStats()
    {
        for ( auto& pair : resultsByPhase )  pair.second.computeStats();
    }

    bool has(unsigned evId, const std::string& stationId, const Catalog::Phase::Type& type ) const
    {
        std::string key = make_key(evId, stationId, type);
        return resultsByPhase.count(key) != 0;
    }

    const Entry& get(unsigned evId, const std::string& stationId,
                     const Catalog::Phase::Type& type ) const
    {
        std::string key = make_key(evId, stationId, type);
        return resultsByPhase.at(key);
    }

    bool has(unsigned evId1, unsigned evId2, const std::string& stationId,
             const Catalog::Phase::Type& type ) const
    {
        return has(evId1, stationId, type) && (get(evId1, stationId, type).peers.count(evId2) != 0);
    }

    const Entry::PeerInfo& get(unsigned evId1, unsigned evId2,
                               const std::string& stationId,
                               const Catalog::Phase::Type& type ) const
    {
        return get(evId1, stationId, type).peers.at(evId2);
    } 

private:
    static std::string make_key(unsigned evId, const std::string& stationId,
                                const Catalog::Phase::Type& type )
    {
        return std::to_string(evId) + "." + stationId + "." + static_cast<char>(type);
    }

    // cache of computed xcorr
    std::unordered_map<std::string, Entry> resultsByPhase;

};


}
}

#endif

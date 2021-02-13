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

#ifndef __HDD_UTILS_H__
#define __HDD_UTILS_H__

#include "catalog.h"
#include <random>
#include <seiscomp3/core/strings.h>
#include <vector>

namespace Seiscomp {
namespace HDD {

template <class T> inline T square(T x) { return x * x; }

double computeDistance(double lat1,
                       double lon1,
                       double depth1,
                       double lat2,
                       double lon2,
                       double depth2,
                       double *azimuth     = nullptr,
                       double *backAzimuth = nullptr);

double computeDistance(double lat1,
                       double lon1,
                       double lat2,
                       double lon2,
                       double *azimuth     = nullptr,
                       double *backAzimuth = nullptr);

double computeDistance(const Catalog::Event &ev1,
                       const Catalog::Event &ev2,
                       double *azimuth     = nullptr,
                       double *backAzimuth = nullptr);

double computeDistance(const Catalog::Event &event,
                       const Catalog::Station &station,
                       double *azimuth     = nullptr,
                       double *backAzimuth = nullptr);

double computeMedian(const std::vector<double> &values);

double computeMedianAbsoluteDeviation(const std::vector<double> &values,
                                      const double median);

double computeMean(const std::vector<double> &values);

double computeMeanAbsoluteDeviation(const std::vector<double> &values,
                                    const double mean);

class Randomer
{

public:
  Randomer(size_t min, size_t max, unsigned int seed = std::random_device{}())
      : gen_{seed}, dist_{min, max}
  {}

  // if you want predictable numbers
  void setSeed(unsigned int seed) { gen_.seed(seed); }

  size_t next() { return dist_(gen_); }

private:
  // random seed by default
  std::mt19937 gen_;
  std::uniform_int_distribution<size_t> dist_;
};

/*
 *  Convert some hashable id of type T (e.g. `std::string`) to an alternative
 *  representation i.e. a sequentially growing integer starting from 0
 *  (suitable for array index).
 */
template <class T> class IdToIndex
{
public:
  unsigned convert(const T &id)
  {
    unsigned idx;
    if (hasId(id, idx)) return idx;
    unsigned newIdx = _currentIdx++;
    _to[id]         = newIdx;
    _from[newIdx]   = id;
    return newIdx;
  }

  unsigned toIdx(const T &id) const { return _to.at(id); }
  T fromIdx(unsigned idx) const { return _from.at(idx); }

  bool hasIdx(unsigned idx) const { return _from.find(idx) != _from.end(); }
  bool hasId(const T &id) const { return _to.find(id) != _to.end(); }

  bool hasIdx(unsigned idx, T &id) const
  {
    const auto &iter = _from.find(idx);
    if (iter == _from.end()) return false;
    id = iter->second;
    return true;
  }

  bool hasId(const T &id, unsigned &idx) const
  {
    const auto &iter = _to.find(id);
    if (iter == _to.end()) return false;
    idx = iter->second;
    return true;
  }

  unsigned size() const { return _to.size(); }

private:
  unsigned _currentIdx = 0;
  std::unordered_map<T, unsigned> _to;
  std::unordered_map<unsigned, T> _from;
};

} // namespace HDD
} // namespace Seiscomp

#endif

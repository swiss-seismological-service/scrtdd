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

#ifndef __HDD_UTILS_H__
#define __HDD_UTILS_H__

#include "catalog.h"
#include <initializer_list>
#include <random>
#include <regex>
#include <stdexcept>
#include <vector>

#include <seiscomp3/datamodel/databasequery.h>

using namespace Seiscomp;

namespace HDD {

template <typename T> inline T square(T x) { return x * x; }

template <typename... Args> std::string strf(Args &&... args)
{
  return Seiscomp::Core::stringify(std::forward<Args>(args)...);
}

std::vector<std::string> splitString(const std::string &str,
                                     const std::regex &regex);

DataModel::SensorLocation *findSensorLocation(const std::string &networkCode,
                                              const std::string &stationCode,
                                              const std::string &locationCode,
                                              const Core::Time &atTime);

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

class Exception : public std::runtime_error
{
public:
  Exception(const std::string &message) : std::runtime_error(message) {}
};

class Logger
{
public:
  static Logger &getInstance()
  {
    static Logger instance; // Guaranteed to be destroyed.
                            // Instantiated on first use.
    return instance;
  }
  Logger(const Logger &) = delete;
  void operator=(const Logger &) = delete;

  enum class Level
  {
    debug,
    info,
    warning,
    error,
  };

  template <typename... Args> void log(Level l, Args &&... args)
  {
    _log(l, strf(std::forward<Args>(args)...));
  }

  void logToFile(const std::string &logFile, const std::vector<Level> &levels);

private:
  Logger() = default;
  void _log(Level, const std::string &);
};

template <typename... Args> void logDebug(Args &&... args)
{
  Logger::getInstance().log(Logger::Level::debug, std::forward<Args>(args)...);
}

template <typename... Args> void logInfo(Args &&... args)
{
  Logger::getInstance().log(Logger::Level::info, std::forward<Args>(args)...);
}

template <typename... Args> void logWarning(Args &&... args)
{
  Logger::getInstance().log(Logger::Level::warning,
                            std::forward<Args>(args)...);
}

template <typename... Args> void logError(Args &&... args)
{
  Logger::getInstance().log(Logger::Level::error, std::forward<Args>(args)...);
}

class UniformRandomer
{

public:
  UniformRandomer(size_t min,
                  size_t max,
                  unsigned int seed = std::random_device{}())
      : _gen(seed), _dist(min, max)
  {}

  // if you want predictable numbers
  void setSeed(unsigned int seed) { _gen.seed(seed); }

  size_t next() { return _dist(_gen); }

private:
  // random seed by default
  std::mt19937 _gen;
  std::uniform_int_distribution<size_t> _dist;
};

class NormalRandomer
{

public:
  NormalRandomer(double mean,
                 double stdDev,
                 unsigned int seed = std::random_device{}())
      : _gen(seed), _dist(mean, stdDev)
  {}

  // if you want predictable numbers
  void setSeed(unsigned int seed) { _gen.seed(seed); }

  double next() { return _dist(_gen); }

private:
  // random seed by default
  std::mt19937 _gen;
  std::normal_distribution<double> _dist;
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

#endif

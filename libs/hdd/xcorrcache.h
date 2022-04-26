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

#ifndef __HDD_XCORRCACHE_H__
#define __HDD_XCORRCACHE_H__

#include "catalog.h"

#include <functional>
#include <string>
#include <unordered_map>

namespace HDD {

class XCorrCache
{
public:
  struct Entry
  {
    bool valid;
    double coeff, lag;
    std::string component;
  };

  XCorrCache()  = default;
  ~XCorrCache() = default;

  XCorrCache(const XCorrCache &other) = default;
  XCorrCache &operator=(const XCorrCache &other) = default;

  XCorrCache(XCorrCache &&other) = default;
  XCorrCache &operator=(XCorrCache &&other) = default;

  void add(unsigned evId1,
           unsigned evId2,
           const std::string &stationId,
           const Catalog::Phase::Type &type,
           bool valid,
           double coeff,
           double lag,
           const std::string component)
  {
    Entry e;
    e.valid                                 = valid;
    e.component                             = component;
    e.coeff                                 = coeff;
    e.lag                                   = lag;
    _entries[evId1][stationId][type][evId2] = e;
    e.lag                                   = -lag;
    _entries[evId2][stationId][type][evId1] = e;
  }

  void remove(unsigned evId1,
              unsigned evId2,
              const std::string &stationId,
              const Catalog::Phase::Type &type)
  {
    try
    {
      _entries.at(evId1).at(stationId).at(type).erase(evId2);
    }
    catch (...)
    {}
    try
    {
      _entries.at(evId2).at(stationId).at(type).erase(evId1);
    }
    catch (...)
    {}
  }

  void remove(unsigned evId1,
              const std::string &stationId,
              const Catalog::Phase::Type &type)
  {
    if (has(evId1, stationId, type))
    {
      std::unordered_map<unsigned, Entry> copy = get(evId1, stationId, type);
      for (const auto &kv : copy)
      {
        remove(evId1, kv.first, stationId, type);
      }
    }
  }

  bool has(unsigned evId1,
           unsigned evId2,
           const std::string &stationId,
           const Catalog::Phase::Type &type) const
  {
    try
    {
      _entries.at(evId1).at(stationId).at(type).at(evId2);
      return true;
    }
    catch (...)
    {
      return false;
    }
  }

  bool has(unsigned evId1,
           const std::string &stationId,
           const Catalog::Phase::Type &type) const
  {
    try
    {
      _entries.at(evId1).at(stationId).at(type);
      return true;
    }
    catch (...)
    {
      return false;
    }
  }

  const Entry &get(unsigned evId1,
                   unsigned evId2,
                   const std::string &stationId,
                   const Catalog::Phase::Type &type) const
  {
    return _entries.at(evId1).at(stationId).at(type).at(evId2);
  }

  const std::unordered_map<unsigned, Entry> &
  get(unsigned evId1,
      const std::string &stationId,
      const Catalog::Phase::Type &type) const
  {
    return _entries.at(evId1).at(stationId).at(type);
  }

  using ForEachCallback = std::function<void(unsigned ev1,
                                             unsigned ev2,
                                             const std::string &stationId,
                                             const Catalog::Phase::Type &type,
                                             const Entry &e)>;

  void forEach(const ForEachCallback &c) const
  {
    for (const auto &kv1 : _entries)
      for (const auto &kv2 : kv1.second)
        for (const auto &kv3 : kv2.second)
          for (const auto &kv4 : kv3.second)
          {
            c(kv1.first, kv4.first, kv2.first, kv3.first, kv4.second);
          }
  }

  void forEach(unsigned evId1, const ForEachCallback &c) const
  {
    try
    {
      for (const auto &kv2 : _entries.at(evId1))
        for (const auto &kv3 : kv2.second)
          for (const auto &kv4 : kv3.second)
          {
            c(evId1, kv4.first, kv2.first, kv3.first, kv4.second);
          }
    }
    catch (...)
    {}
  }

  void forEach(unsigned evId1,
               const std::string &stationId,
               const ForEachCallback &c) const
  {
    try
    {
      for (const auto &kv3 : _entries.at(evId1).at(stationId))
        for (const auto &kv4 : kv3.second)
        {
          c(evId1, kv4.first, stationId, kv3.first, kv4.second);
        }
    }
    catch (...)
    {}
  }

  void forEach(unsigned evId1,
               const std::string &stationId,
               const Catalog::Phase::Type &type,
               const ForEachCallback &c) const
  {
    try
    {
      for (const auto &kv4 : _entries.at(evId1).at(stationId).at(type))
      {
        c(evId1, kv4.first, stationId, type, kv4.second);
      }
    }
    catch (...)
    {}
  }

private:
  std::unordered_map<
      unsigned, // indexed by event id
      std::unordered_map<
          std::string,                             // indexed by station id
          std::unordered_map<Catalog::Phase::Type, // indexed by phase type
                             std::unordered_map<unsigned, // indexed by event id
                                                Entry>>>>
      _entries;
};

} // namespace HDD

#endif

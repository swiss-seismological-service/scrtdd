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
  };

  bool empty() const { return _entries.empty(); }

  void add(unsigned evId1,
           unsigned evId2,
           const std::string &stationId,
           const std::string &phase,
           bool valid,
           double coeff,
           double lag)
  {
    Entry e;
    e.valid                                  = valid;
    e.coeff                                  = coeff;
    e.lag                                    = lag;
    _entries[evId1][stationId][phase][evId2] = e;
    e.lag                                    = -lag;
    _entries[evId2][stationId][phase][evId1] = e;
  }

  void add(const XCorrCache &xcorr)
  {
    auto callback =
        [this](unsigned ev1, unsigned ev2, const std::string &stationId,
               const std::string &phase, const XCorrCache::Entry &e) {
          add(ev1, ev2, stationId, phase, e.valid, e.coeff, e.lag);
        };
    xcorr.forEach(callback);
  }

  void remove(unsigned evId1,
              unsigned evId2,
              const std::string &stationId,
              const std::string &phase)
  {
    try
    {
      _entries.at(evId1).at(stationId).at(phase).erase(evId2);
    }
    catch (...)
    {}
    try
    {
      _entries.at(evId2).at(stationId).at(phase).erase(evId1);
    }
    catch (...)
    {}
  }

  void
  remove(unsigned evId1, const std::string &stationId, const std::string &phase)
  {
    if (has(evId1, stationId, phase))
    {
      std::unordered_map<unsigned, Entry> copy = get(evId1, stationId, phase);
      for (const auto &kv : copy)
      {
        remove(evId1, kv.first, stationId, phase);
      }
    }
  }

  bool has(unsigned evId1,
           unsigned evId2,
           const std::string &stationId,
           const std::string &phase) const
  {
    try
    {
      _entries.at(evId1).at(stationId).at(phase).at(evId2);
      return true;
    }
    catch (...)
    {
      return false;
    }
  }

  bool has(unsigned evId1,
           const std::string &stationId,
           const std::string &phase) const
  {
    try
    {
      _entries.at(evId1).at(stationId).at(phase);
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
                   const std::string &phase) const
  {
    return _entries.at(evId1).at(stationId).at(phase).at(evId2);
  }

  const std::unordered_map<unsigned, Entry> &get(unsigned evId1,
                                                 const std::string &stationId,
                                                 const std::string &phase) const
  {
    return _entries.at(evId1).at(stationId).at(phase);
  }

  using ForEachCallback = std::function<void(unsigned ev1,
                                             unsigned ev2,
                                             const std::string &stationId,
                                             const std::string &phase,
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
               const std::string &phase,
               const ForEachCallback &c) const
  {
    try
    {
      for (const auto &kv4 : _entries.at(evId1).at(stationId).at(phase))
      {
        c(evId1, kv4.first, stationId, phase, kv4.second);
      }
    }
    catch (...)
    {}
  }

  void writeToFile(const Catalog &cat, const std::string &file) const;

  static XCorrCache readFromFile(const Catalog &cat, const std::string &file);

private:
  std::unordered_map<
      unsigned, // indexed by event id
      std::unordered_map<
          std::string,                    // indexed by station id
          std::unordered_map<std::string, // indexed by phase type
                             std::unordered_map<unsigned, // indexed by event id
                                                Entry>>>>
      _entries;
};

} // namespace HDD

#endif

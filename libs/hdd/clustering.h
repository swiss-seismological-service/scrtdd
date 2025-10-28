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

#ifndef __HDD_CLUSTERING_H__
#define __HDD_CLUSTERING_H__

#include "catalog.h"
#include <fstream>
#include <list>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace HDD {

class Neighbours
{
public:
  Neighbours(unsigned refEvId) : _refEvId(refEvId){};

  Neighbours(const Neighbours &other)            = default;
  Neighbours &operator=(const Neighbours &other) = default;

  Neighbours(Neighbours &&other)            = default;
  Neighbours &operator=(Neighbours &&other) = default;

  ~Neighbours() = default;

  unsigned referenceId() const { return _refEvId; }
  void setReferenceId(unsigned refEvId) { _refEvId = refEvId; }

  const std::unordered_set<unsigned> &ids() const { return _ids; }

  void add(unsigned neighbourId,
           const std::string &stationId,
           const Catalog::Phase::Type &phase)
  {
    _ids.insert(neighbourId);
    _phases[neighbourId][stationId].insert(phase);
  }

  void remove(unsigned neighbourId)
  {
    _ids.erase(neighbourId);
    _phases.erase(neighbourId);
  }

  bool has(unsigned neighbourId) const
  {
    return _ids.find(neighbourId) != _ids.end();
  }

  bool has(unsigned neighbourId, const std::string stationId) const
  {
    const auto &neighPhases = _phases.find(neighbourId);
    if (neighPhases != _phases.end())
      return neighPhases->second.find(stationId) != neighPhases->second.end();
    return false;
  }

  bool has(unsigned neighbourId,
           const std::string stationId,
           Catalog::Phase::Type type) const
  {
    const auto &neighPhases = _phases.find(neighbourId);
    if (neighPhases != _phases.end())
    {
      const auto &neighPhaseTypes = neighPhases->second.find(stationId);
      if (neighPhaseTypes != neighPhases->second.end())
        return neighPhaseTypes->second.find(type) !=
               neighPhaseTypes->second.end();
    }
    return false;
  }

  std::unordered_map<std::string, std::unordered_set<Catalog::Phase::Type>>
  allPhases() const;

  std::unordered_map<std::string, std::unordered_set<Catalog::Phase::Type>>
  allPhases(unsigned neighbourId) const;

  std::unique_ptr<Catalog> toCatalog(const Catalog &catalog,
                                     bool includeRefEv = false) const;

  void writeToFile(const Catalog &cat, const std::string &file) const;
  void writeToFile(const Catalog &cat, std::ostream &os) const;

  static std::unordered_map<unsigned, std::unique_ptr<Neighbours>>
  readFromFile(const Catalog &cat, const std::string &file);

private:
  unsigned _refEvId;
  std::unordered_set<unsigned> _ids; // neighbouring event id
  std::unordered_map<
      unsigned,                       // indexed by event id
      std::unordered_map<std::string, // indexed by station id
                         std::unordered_set<Catalog::Phase::Type>>>
      _phases;
};

std::unique_ptr<Neighbours>
selectNeighbouringEvents(const Catalog &catalog,
                         const Catalog::Event &refEv,
                         const Catalog &refEvCatalog,
                         double minPhaseWeight   = 0,
                         double minESdis         = 0,
                         double maxESdis         = -1, // -1 = no limits
                         double minEStoIEratio   = 0,
                         unsigned minDTperEvt    = 1,
                         unsigned maxDTperEvt    = 0, // 0 = no limits
                         unsigned minNumNeigh    = 1,
                         unsigned maxNumNeigh    = 0, // 0 = no limits
                         unsigned numEllipsoids  = 5,
                         double maxEllipsoidSize = 10,
                         bool keepUnmatched      = false);

std::unordered_map<unsigned, std::unique_ptr<Neighbours>>
selectNeighbouringEventsCatalog(const Catalog &catalog,
                                double minPhaseWeight,
                                double minESdis,
                                double maxESdis,
                                double minEStoIEratio,
                                unsigned minDTperEvt,
                                unsigned maxDTperEvt,
                                unsigned minNumNeigh,
                                unsigned maxNumNeigh,
                                unsigned numEllipsoids,
                                double maxEllipsoidSize,
                                bool keepUnmatched);

std::list<std::unordered_map<unsigned, std::unique_ptr<Neighbours>>>
clusterizeNeighbouringEvents(
    std::unordered_map<unsigned, std::unique_ptr<Neighbours>>
        &neighboursByEvent);

} // namespace HDD

#endif

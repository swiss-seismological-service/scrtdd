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
#include "kdtree.h"
#include <fstream>
#include <list>
#include <set>
#include <tuple>
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

  size_t amount() const;

  std::unordered_set<unsigned> ids() const;

  void add(unsigned neighbourId,
           const std::string &stationId,
           const std::string &phase);

  void remove(unsigned neighbourId);

  bool has(unsigned neighbourId) const;

  bool has(unsigned neighbourId, const std::string &stationId) const;

  bool has(unsigned neighbourId,
           const std::string &stationId,
           const std::string &phase) const;

  std::unordered_set<std::string> stations() const;

  std::vector<std::tuple<std::string, std::string, unsigned>> phases() const;

  std::vector<std::tuple<std::string, std::string>>
  phases(unsigned neighbourId) const;

  Catalog toCatalog(const Catalog &catalog, bool includeRefEv = false) const;

  std::ofstream writeToFile(const Catalog &cat, const std::string &file) const;
  void appendToStream(const Catalog &cat, std::ostream &os) const;

  static void
  writeToFile(const std::unordered_map<unsigned, Neighbours> &neighboursByEvent,
              const Catalog &cat,
              const std::string &file);

  static std::unordered_map<unsigned, Neighbours>
  readFromFile(const Catalog &cat, const std::string &file);

private:
  unsigned _refEvId;
  std::unordered_map<
      unsigned,                       // indexed by event id
      std::unordered_map<std::string, // indexed by phase type
                         std::unordered_set<std::string>>> // station id
      _phases;
};

using EventTree = KDTree<unsigned>; // store event id

EventTree createEventTree(const Catalog &catalog);

Neighbours selectNeighbouringEvents(const EventTree &tree,
                                    const Catalog &catalog,
                                    const Catalog::Event &refEv,
                                    const Catalog &refEvCatalog,
                                    double minESdis = 0,
                                    double maxESdis = -1, // -1 = no limits
                                    double minEStoIEratio  = 0,
                                    unsigned minPhase      = 1,
                                    unsigned maxPhase      = 0, // 0 = no limits
                                    unsigned minNumNeigh   = 1,
                                    unsigned maxNumNeigh   = 0, // 0 = no limits
                                    unsigned numEllipsoids = 5,
                                    double maxNeighbourDist = 5); // km

std::unordered_map<unsigned, Neighbours>
selectNeighbouringEventsCatalog(const EventTree &tree,
                                const Catalog &catalog,
                                double minESdis,
                                double maxESdis,
                                double minEStoIEratio,
                                unsigned minPhase,
                                unsigned maxPhase,
                                unsigned minNumNeigh,
                                unsigned maxNumNeigh,
                                unsigned numEllipsoids,
                                double maxNeighbourDist);

std::list<std::unordered_map<unsigned, Neighbours>>
clusterizeNeighbouringEvents(
    std::unordered_map<unsigned, Neighbours> neighboursByEvent);

} // namespace HDD

#endif

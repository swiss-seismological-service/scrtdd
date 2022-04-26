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

#ifndef __HDD_INDEX_H__
#define __HDD_INDEX_H__

#include <unordered_map>

namespace HDD {

/*
 *  Convert some hashable id of type T (e.g. `std::string`) to a
 *  a sequentially growing integer starting from 0
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

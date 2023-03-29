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

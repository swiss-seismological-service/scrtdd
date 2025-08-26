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

#ifndef __HDD_TTT_KDTREE_H__
#define __HDD_TTT_KDTREE_H__

#include "../utils.h"

#include <memory>
#include <stdexcept>
#include <vector>

namespace HDD {

namespace TTT {

template <typename T> class KDTree
{
public:
  struct Point
  {
    double latitude, longitude, elevation; // km
    T data;
  };

  KDTree() = default;

  KDTree(const std::vector<Point> &points)
  {
    _points = points;

    std::vector<size_t> indices(points.size());
    std::iota(std::begin(indices), std::end(indices), 0);

    _root = buildRecursive(indices.data(), points.size(), 0);
  }

  const Point &search(const Point &query) const
  {
    size_t idx;
    if (!searchRecursive(query, _root.get(), &idx))
    {
      throw std::range_error("There is no such point in the kd-tree");
    }
    return _points[idx];
  }

  const Point &search(double latitude, double longitude, double elevation) const
  {
    Point query;
    query.latitude  = latitude;
    query.longitude = longitude;
    query.elevation = elevation;
    return search(query);
  }

private:
  struct Node
  {
    size_t idx = -1;               // index to the original point
    int axis   = -1;               // dimension's axis
    std::unique_ptr<Node> next[2]; // pointers to the child nodes
  };

  std::unique_ptr<Node>
  buildRecursive(size_t *indices, size_t npoints, int depth)
  {
    if (npoints == 0) return nullptr;

    const int axis   = depth % 3; // lat, lon, elevation
    const size_t mid = (npoints - 1) / 2;

    std::nth_element(
        indices, indices + mid, indices + npoints, [&](int lhs, int rhs) {
          if (axis == 0)
            return _points[lhs].latitude < _points[rhs].latitude;
          else if (axis == 1)
            return _points[lhs].longitude < _points[rhs].longitude;
          else if (axis == 2)
            return _points[lhs].elevation < _points[rhs].elevation;
          else
            throw std::runtime_error("KDTree internal logic error");
        });

    std::unique_ptr<Node> node(new Node());
    node->idx  = indices[mid];
    node->axis = axis;

    node->next[0] = buildRecursive(indices, mid, depth + 1);
    node->next[1] =
        buildRecursive(indices + mid + 1, npoints - mid - 1, depth + 1);

    return node;
  }

  bool searchRecursive(const Point &query, const Node *node, size_t *idx) const
  {
    if (!node)
    {
      return false;
    }

    const Point &curr = _points[node->idx];

    double dist =
        computeDistance(curr.latitude, curr.longitude, curr.elevation,
                        query.latitude, query.longitude, query.elevation);
    if (dist < 0.0001) // 10 cm precision
    {
      *idx = node->idx;
      return true;
    }

    const int axis = node->axis;
    const Node *next_node;
    if (axis == 0)
      next_node = node->next[query.latitude < curr.latitude ? 0 : 1].get();
    else if (axis == 1)
      next_node = node->next[query.longitude < curr.longitude ? 0 : 1].get();
    else if (axis == 2)
      next_node = node->next[query.elevation < curr.elevation ? 0 : 1].get();
    else
      throw std::runtime_error("KDTree internal logic error");

    return searchRecursive(query, next_node, idx);
  }

private:
  std::unique_ptr<Node> _root;
  std::vector<Point> _points;
};

} // namespace TTT

} // namespace HDD

#endif

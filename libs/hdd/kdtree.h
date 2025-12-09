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

#include "utils.h"

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
    double latitude;
    double longitude;
    double depth;
    T data;
  };

  KDTree() = default;

  KDTree(const std::vector<Point> &points)
  {
    _points = points;
    _nodes.resize(points.size());

    std::vector<size_t> indices(_points.size());
    std::iota(std::begin(indices), std::end(indices), 0);

    _root = buildRecursive(indices.data(), _points.size(), 0);
  }

  const Point &search(const Point &query) const
  {
    size_t idx;
    if (!searchRecursive(query, _root, idx))
    {
      throw std::range_error("There is no such point in the kd-tree");
    }
    return _points[idx];
  }

  const Point &search(double latitude, double longitude, double depth) const
  {
    Point query;
    query.latitude  = latitude;
    query.longitude = longitude;
    query.depth     = depth;
    return search(query);
  }

  const Point &nnSearch(const Point &query, double &dist) const
  {
    size_t idx;
    dist = std::numeric_limits<double>::max();
    if (!nnSearchRecursive(query, _root, idx, dist))
    {
      throw std::range_error("kd-tree is empty");
    }
    return _points[idx];
  }

  const Point &
  nnSearch(double latitude, double longitude, double depth, double &dist) const
  {
    Point query;
    query.latitude  = latitude;
    query.longitude = longitude;
    query.depth     = depth;
    return nnSearch(query, dist);
  }

private:
  struct Node
  {
    size_t idx;    // index to the original point
    int axis;      // dimension's axis
    Node *next[2]; // index to the child nodes
  };

  Node *buildRecursive(size_t *indices, size_t npoints, int depth)
  {
    if (npoints == 0) return nullptr;

    const int axis   = depth % 3; // lat, lon, depth
    const size_t mid = (npoints - 1) / 2;

    std::nth_element(
        indices, indices + mid, indices + npoints, [&](int lhs, int rhs) {
          if (axis == 0)
            return _points[lhs].latitude < _points[rhs].latitude;
          else if (axis == 1)
            return _points[lhs].longitude < _points[rhs].longitude;
          else if (axis == 2)
            return _points[lhs].depth < _points[rhs].depth;
          else
            throw std::runtime_error("KDTree internal logic error");
        });

    const size_t idx = indices[mid];

    Node *node    = &_nodes.at(idx);
    node->idx     = idx;
    node->axis    = axis;
    node->next[0] = buildRecursive(indices, mid, depth + 1);
    node->next[1] =
        buildRecursive(indices + mid + 1, npoints - mid - 1, depth + 1);

    return node;
  }

  bool searchRecursive(const Point &query, const Node *node, size_t &idx) const
  {
    if (!node)
    {
      return false;
    }

    const Point &curr = _points[node->idx];

    if (curr.latitude == query.latitude && curr.longitude == query.longitude &&
        curr.depth == query.depth)
    {
      idx = node->idx;
      return true;
    }

    const int axis = node->axis;
    int dir;
    bool equalOnAxis;
    if (axis == 0)
    {
      dir         = query.latitude < curr.latitude ? 0 : 1;
      equalOnAxis = query.latitude == curr.latitude;
    }
    else if (axis == 1)
    {
      dir         = query.longitude < curr.longitude ? 0 : 1;
      equalOnAxis = query.longitude == curr.longitude;
    }
    else if (axis == 2)
    {
      dir         = query.depth < curr.depth ? 0 : 1;
      equalOnAxis = query.depth == curr.depth;
    }
    else
    {
      throw std::runtime_error("KDTree internal logic error");
    }

    bool found = searchRecursive(query, node->next[dir], idx);
    if (!found && equalOnAxis)
    {
      return searchRecursive(query, node->next[dir == 0 ? 1 : 0], idx);
    }
    return found;
  }

  bool nnSearchRecursive(const Point &query,
                         const Node *node,
                         size_t &idx,
                         double &minDist) const
  {
    if (!node)
    {
      return false;
    }

    const Point &curr = _points[node->idx];

    double dist = computeDistance(curr.latitude, curr.longitude, curr.depth,
                                  query.latitude, query.longitude, query.depth);

    if (dist < minDist)
    {
      idx     = node->idx;
      minDist = dist;
    }

    const int axis = node->axis;
    double axisDist;
    int dir;
    if (axis == 0)
    {
      dir      = query.latitude < curr.latitude ? 0 : 1;
      axisDist = computeDistance(curr.latitude, 0, query.latitude, 0);
    }
    else if (axis == 1)
    {
      dir      = query.longitude < curr.longitude ? 0 : 1;
      axisDist = computeDistance(0, query.longitude, 0, query.longitude);
    }
    else if (axis == 2)
    {
      dir      = query.depth < curr.depth ? 0 : 1;
      axisDist = std::abs(curr.depth - query.depth);
    }
    else
    {
      throw std::runtime_error("KDTree internal logic error");
    }

    nnSearchRecursive(query, node->next[dir], idx, minDist);

    if (axisDist < minDist)
    {
      nnSearchRecursive(query, node->next[dir == 0 ? 1 : 0], idx, minDist);
    }
    return true;
  }

private:
  Node *_root;
  std::vector<Node> _nodes;
  std::vector<Point> _points;
};

} // namespace TTT

} // namespace HDD

#endif

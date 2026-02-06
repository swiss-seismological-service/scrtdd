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

#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <vector>

namespace HDD {

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

  KDTree(std::vector<Point> points)
  {
    _points = std::move(points);
    _nodes.resize(_points.size());

    std::vector<size_t> indices(_points.size());
    std::iota(std::begin(indices), std::end(indices), 0);

    _root = buildRecursive(indices.data(), _points.size(), 0);
  }

  KDTree()  = default;
  ~KDTree() = default;

  KDTree(const KDTree &other) : KDTree(other.points()) {}

  // Custom copy assignment operator
  KDTree &operator=(const KDTree &other)
  {
    return *this = KDTree(other.points());
  }

  KDTree(KDTree &&other)            = default;
  KDTree &operator=(KDTree &&other) = default;

  const std::vector<Point> &points() const { return _points; }

  const Point &search(double latitude, double longitude, double depth) const
  {
    Point query;
    query.latitude  = latitude;
    query.longitude = longitude;
    query.depth     = depth;
    size_t idx;
    if (!searchRecursive(query, _root, idx))
    {
      throw std::range_error("There is no such point in the kd-tree");
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
    size_t idx;
    dist = std::numeric_limits<double>::max();
    if (!nnSearchRecursive(query, _root, idx, dist))
    {
      throw std::range_error("kd-tree is empty");
    }
    return _points[idx];
  }

  std::multimap<double, size_t> // distance km, point idx
  knnSearch(double latitude, double longitude, double depth, size_t k) const
  {
    Point query;
    query.latitude  = latitude;
    query.longitude = longitude;
    query.depth     = depth;
    std::multimap<double, size_t> queue;
    knnSearchRecursive(query, _root, queue, k);
    return queue;
  }

  std::multimap<double, size_t> // distance km, point idx
  radiusSearch(double latitude,
               double longitude,
               double depth,
               double radius) const
  {
    Point query;
    query.latitude  = latitude;
    query.longitude = longitude;
    query.depth     = depth;
    std::multimap<double, size_t> queue;
    radiusSearchRecursive(query, _root, queue, radius);
    return queue;
  }

  void
  bestFirstSearch(double latitude,
                  double longitude,
                  double depth,
                  std::function<bool(const Point &, double)> callback) const
  {
    if (!_root) return;

    struct Item
    {
      double minDist;
      const Node *node;   // If non-null, this is a node to expand
      const Point *point; // If non-null, this is a result point

      // Priority queue is a max-heap, we want a min-heap (lowest distance
      // first)
      bool operator>(const Item &other) const
      {
        return minDist > other.minDist;
      }
    };

    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;

    Point query{latitude, longitude, depth, T()};

    // Start with the root node
    pq.push({0.0, _root, nullptr});

    while (!pq.empty())
    {
      Item top = pq.top();
      pq.pop();

      // If we popped a Point, it is guaranteed to be the next closest point.
      if (top.point)
      {
        if (callback(*top.point, top.minDist))
        {
          return;
        }
        continue;
      }

      // If we popped a Node, expand it
      const Node *node  = top.node;
      const Point &curr = _points[node->idx];

      // Add the actual point at this node to the queue as a candidate
      double d_point =
          computeDistance(curr.latitude, curr.longitude, curr.depth,
                          query.latitude, query.longitude, query.depth);
      pq.push({d_point, nullptr, &curr});

      int axis = node->axis;
      int dir;
      double axisDist;

      if (axis == 0)
      {
        dir      = query.latitude < curr.latitude ? 0 : 1;
        axisDist = computeDistance(curr.latitude, 0, query.latitude, 0);
      }
      else if (axis == 1)
      {
        dir      = query.longitude < curr.longitude ? 0 : 1;
        axisDist = computeDistance(0, curr.longitude, 0, query.longitude);
      }
      else
      {
        dir      = query.depth < curr.depth ? 0 : 1;
        axisDist = std::abs(curr.depth - query.depth);
      }

      if (node->next[dir])
      {
        // The near volume has the same min distance as the current volume
        pq.push({top.minDist, node->next[dir], nullptr});
      }

      if (node->next[1 - dir])
      {
        // The far volume's min distance is at least the distance to the split
        // plane
        double farMinDist = std::max(top.minDist, axisDist);
        pq.push({farMinDist, node->next[1 - dir], nullptr});
      }
    }
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
        indices, indices + mid, indices + npoints, [&](size_t lhs, size_t rhs) {
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
      return searchRecursive(query, node->next[1 - dir], idx);
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

    const double dist =
        computeDistance(curr.latitude, curr.longitude, curr.depth,
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
      axisDist = computeDistance(0, curr.longitude, 0, query.longitude);
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
      nnSearchRecursive(query, node->next[1 - dir], idx, minDist);
    }
    return true;
  }

  void knnSearchRecursive(const Point &query,
                          const Node *node,
                          std::multimap<double, size_t> &queue, // dist km, idx
                          size_t k) const
  {
    if (!node)
    {
      return;
    }

    const Point &curr = _points[node->idx];

    const double dist =
        computeDistance(curr.latitude, curr.longitude, curr.depth,
                        query.latitude, query.longitude, query.depth);

    if (queue.size() < k)
    {
      queue.emplace(dist, node->idx);
    }
    else if (dist < queue.rbegin()->first)
    {
      queue.emplace(dist, node->idx);
      queue.erase(std::prev(queue.end())); // remove last entry
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
      axisDist = computeDistance(0, curr.longitude, 0, query.longitude);
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

    knnSearchRecursive(query, node->next[dir], queue, k);

    if (queue.size() < k || axisDist < queue.rbegin()->first)
    {
      knnSearchRecursive(query, node->next[1 - dir], queue, k);
    }
  }

  void
  radiusSearchRecursive(const Point &query,
                        const Node *node,
                        std::multimap<double, size_t> &queue, // dist km, idx
                        double radius) const
  {
    if (!node)
    {
      return;
    }

    const Point &curr = _points[node->idx];

    const double dist =
        computeDistance(curr.latitude, curr.longitude, curr.depth,
                        query.latitude, query.longitude, query.depth);
    if (dist <= radius)
    {
      queue.emplace(dist, node->idx);
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
      axisDist = computeDistance(0, curr.longitude, 0, query.longitude);
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

    radiusSearchRecursive(query, node->next[dir], queue, radius);

    if (axisDist <= radius)
    {
      radiusSearchRecursive(query, node->next[1 - dir], queue, radius);
    }
  }

private:
  Node *_root;
  std::vector<Node> _nodes;
  std::vector<Point> _points;
};

} // namespace HDD

#endif

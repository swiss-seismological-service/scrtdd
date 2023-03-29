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

#ifndef __HDD_RANDOM_H__
#define __HDD_RANDOM_H__

#include <random>

namespace HDD {

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

} // namespace HDD

#endif

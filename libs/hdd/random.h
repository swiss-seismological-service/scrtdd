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

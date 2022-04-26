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

#include "xcorr.h"
#include "log.h"

#include <cfenv>
#include <vector>

using namespace std;

namespace HDD {

void crossCorrelation(const double *dataS,
                      const int sizeS,
                      const double *dataL,
                      const int sizeL,
                      double &delayOut,
                      double &coeffOut)
{
  /*
   * Pearson correlation coefficient for time series X and Y of length n
   *
   *              sum((Xi-meanX) * (Yi-meanY))
   * cc = --------------------------------------------------
   *      sqrt(sum((Xi-meanX)^2)) * sqrt(sum((Yi-meanY)^2))
   *
   * Where sum(X)  is the sum of Xi for i=1 until i=n
   *
   * This can be rearranged in a form suitable for a single-pass algorithm
   * (where the mean of X and Y are not needed)
   *
   *                 n * sum(Xi*Yi) - sum(Xi) * sum(Yi)
   * cc = -----------------------------------------------------------
   *      sqrt(n*sum(Xi^2)-sum(Xi)^2) * sqrt(n*sum(Yi^2)-sum(Yi)^2))
   *
   * For cross-correlation, where we have a short trace S which is correlated
   * against a longer trace L at subsequent offset, we can pre-compute the
   * parts that involves S and re-use them at each step of the
   * cross-correlation:
   *
   *   sumS   = sum(Xi)
   *   denomS = sqrt(n*sum(Xi^2)-sum(Xi)^2)
   *
   * For the parts that involves the longer trace L alone we can compute them
   * in a rolling fashion (removing first sample of previous iteration and
   * adding the last sample of the new iteration):
   *
   *   sumL   = sum(Yi)
   *   sumL2  = sum(Yi^2)
   *   denomL = sqrt(n*sumL2-sumL^2))
   *
   * Finally, this is the equation at each step (offset) of cross-correlation:
   *
   *       n * sum(Xi*Yi) - sumS * sumL
   * cc = ------------------------------
   *             denomS * denomL
   *
   * Unfortunately, we cannot optimize sum(Xi*Yi) and this will be a inner
   * loop inside the main cross-correlation loop.
   */

  std::feclearexcept(FE_ALL_EXCEPT);

  // prepare the data before the main xcorr loop
  const int n = sizeS;
  double sumS = 0, sumS2 = 0;
  double sumL = 0, sumL2 = 0;
  for (int i = 0; i < n; i++)
  {
    sumS += dataS[i];
    sumS2 += dataS[i] * dataS[i];
    if (i >= (n - 1)) continue;
    sumL += dataL[i];
    sumL2 += dataL[i] * dataL[i];
  }
  double denomS = std::sqrt(n * sumS2 - sumS * sumS);

  // cross-correlation loop
  coeffOut           = std::nan("");
  double lastSampleL = 0;
  for (int delay = 0; delay <= (sizeL - sizeS); delay++)
  {
    // sumL/sumL2 update: remove the sample that has just exited the
    // current cross-correlation win and add the sample that has just
    // entered
    const double newSampleL = dataL[delay + n - 1];
    sumL += newSampleL - lastSampleL;
    sumL2 += newSampleL * newSampleL - lastSampleL * lastSampleL;
    lastSampleL = dataL[delay]; // prepare for next loop

    const double denomL = std::sqrt(n * sumL2 - sumL * sumL);

    double sumSL = 0;
    for (int i = 0; i < n; i++) sumSL += dataS[i] * dataL[i + delay];

    const double coeff = (n * sumSL - sumS * sumL) / (denomS * denomL);

    if (!std::isfinite(coeffOut) || std::abs(coeff) > std::abs(coeffOut))
    {
      coeffOut = coeff;
      delayOut = delay;
    }
  }

  int fe = fetestexcept(FE_ALL_EXCEPT);
  if ((fe & ~FE_INEXACT) != 0) // we don't care about FE_INEXACT
  {
    logWarning("Floating point exception during cross-correlation:");
    if (fe & FE_DIVBYZERO) logWarning("FE_DIVBYZERO");
    if (fe & FE_INVALID) logWarning("FE_INVALID");
    if (fe & FE_OVERFLOW) logWarning("FE_OVERFLOW");
    if (fe & FE_UNDERFLOW) logWarning("FE_UNDERFLOW");
  }
}

} // namespace HDD

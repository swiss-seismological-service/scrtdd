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

#ifndef __HDD_SCWAVEFORM_H__
#define __HDD_SCWAVEFORM_H__

#include "hdd/waveform.h"

#include <functional>
#include <string>
#include <unordered_map>

namespace HDD {

namespace SeiscompAdapter {

std::unique_ptr<Trace>
loadTraceFromRecordStream(const std::string &recordStreamURL,
                          const TimeWindow &tw,
                          const std::string &networkCode,
                          const std::string &stationCode,
                          const std::string &locationCode,
                          const std::string &channelCode,
                          double tolerance       = 0.1,
                          double minAvailability = 0.95);

void loadTracesFromRecordStream(
    const std::string &recordStreamURL,
    const std::unordered_multimap<std::string, const TimeWindow> &request,
    const std::function<void(const std::string &,
                             const TimeWindow &,
                             std::unique_ptr<Trace>)> &onTraceLoaded,
    const std::function<void(const std::string &,
                             const TimeWindow &,
                             const std::string &)> &onTraceFailed,
    double tolerance       = 0.1,
    double minAvailability = 0.95);

void getComponentsInfo(const Catalog::Phase &ph,
                       Waveform::ThreeComponents &components);

void filter(Trace &trace, const std::string &filterStr);

void writeTrace(const Trace &trace, const std::string &file);
std::unique_ptr<Trace> readTrace(const std::string &file);

} // namespace SeiscompAdapter
} // namespace HDD

#endif

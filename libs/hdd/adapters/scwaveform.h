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
namespace SCAdapter {

class WaveformProxy : public HDD::Waveform::Proxy
{
public:
  WaveformProxy(const std::string &recordStream = "")
      : _recordStreamURL(recordStream)
  {}
  virtual ~WaveformProxy() = default;

  std::unique_ptr<HDD::Trace>
  loadTrace(const HDD::TimeWindow &tw,
            const std::string &networkCode,
            const std::string &stationCode,
            const std::string &locationCode,
            const std::string &channelCode) override;

  void loadTraces(
      const std::unordered_multimap<std::string, const HDD::TimeWindow>
          &request,
      const std::function<void(const std::string &,
                               const HDD::TimeWindow &,
                               std::unique_ptr<HDD::Trace>)> &onTraceLoaded,
      const std::function<void(const std::string &,
                               const HDD::TimeWindow &,
                               const std::string &)> &onTraceFailed) override;

  void getComponentsInfo(const HDD::Catalog::Phase &ph,
                         HDD::Waveform::ThreeComponents &components) override;

  void filter(HDD::Trace &trace, const std::string &filterStr) override;

  void writeTrace(const HDD::Trace &trace, const std::string &file) override;
  std::unique_ptr<HDD::Trace> readTrace(const std::string &file) override;

private:
  const std::string _recordStreamURL;
  static constexpr double tolerance       = 0.5;
  static constexpr double minAvailability = 0.95;
};

} // namespace SCAdapter
} // namespace HDD
#endif

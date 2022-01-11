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

#include "log.h"

#define SEISCOMP_COMPONENT HDD
#include <seiscomp3/logging/file.h>
#include <seiscomp3/logging/log.h>

using namespace std;

namespace HDD {

void Logger::_log(const Level l, const string &s)
{
  if (l == Level::info)
    SEISCOMP_LOG(Seiscomp::Logging::_SCInfoChannel, "%s", s.c_str());
  else if (l == Level::warning)
    SEISCOMP_LOG(Seiscomp::Logging::_SCWarningChannel, "%s", s.c_str());
  else if (l == Level::error)
    SEISCOMP_LOG(Seiscomp::Logging::_SCErrorChannel, "%s", s.c_str());
  else if (l == Level::debug)
    SEISCOMP_LOG(Seiscomp::Logging::_SCDebugChannel, "%s", s.c_str());
}

void Logger::logToFile(const std::string &logFile,
                       const std::vector<Level> &levels)
{
  Seiscomp::Logging::FileOutput processingInfoOutput(logFile.c_str());
  for (auto l : levels)
  {
    if (l == Level::info)
      processingInfoOutput.subscribe(Seiscomp::Logging::_SCInfoChannel);
    else if (l == Level::warning)
      processingInfoOutput.subscribe(Seiscomp::Logging::_SCWarningChannel);
    else if (l == Level::error)
      processingInfoOutput.subscribe(Seiscomp::Logging::_SCErrorChannel);
    else if (l == Level::debug)
      processingInfoOutput.subscribe(Seiscomp::Logging::_SCDebugChannel);
  }
}

} // namespace HDD

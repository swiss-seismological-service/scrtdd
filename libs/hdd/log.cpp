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

void Logger::_log(Level l, const string &s)
{
  Seiscomp::Logging::Channel *logChannel = Seiscomp::Logging::_SCDebugChannel;
  if (l == Level::info) logChannel = Seiscomp::Logging::_SCInfoChannel;
  if (l == Level::warning) logChannel = Seiscomp::Logging::_SCWarningChannel;
  if (l == Level::error) logChannel = Seiscomp::Logging::_SCErrorChannel;
  SEISCOMP_LOG(logChannel, "%s", s.c_str());
}

void Logger::logToFile(const std::string &logFile,
                       const std::vector<Level> &levels)
{
  Seiscomp::Logging::FileOutput processingInfoOutput(logFile.c_str());
  for (auto l : levels)
  {
    Seiscomp::Logging::Channel *logChannel = Seiscomp::Logging::_SCDebugChannel;
    if (l == Level::info) logChannel = Seiscomp::Logging::_SCInfoChannel;
    if (l == Level::warning) logChannel = Seiscomp::Logging::_SCWarningChannel;
    if (l == Level::error) logChannel = Seiscomp::Logging::_SCErrorChannel;
    processingInfoOutput.subscribe(logChannel);
  }
}

} // namespace HDD

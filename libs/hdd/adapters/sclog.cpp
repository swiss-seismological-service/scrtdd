/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU Affero General Public License as published*
 * by the Free Software Foundation, either version 3 of the License, or    *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This program is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#include "sclog.h"
#include "hdd/log.h"

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp/logging/file.h>
#include <seiscomp/logging/log.h>

using namespace std;
using namespace Seiscomp;

using namespace HDD;

namespace HDD {
namespace SCAdapter {

void initLogger()
{
  auto error   = [](const string &msg) { SEISCOMP_ERROR_S(msg); };
  auto warning = [](const string &msg) { SEISCOMP_WARNING_S(msg); };
  auto info    = [](const string &msg) { SEISCOMP_INFO_S(msg); };
  auto debug   = [](const string &msg) { SEISCOMP_DEBUG_S(msg); };

  auto createFileLogger = [](const std::string &logFile,
                             const std::vector<HDD::Logger::Level> &levels) {
    Logging::FileOutput *processingInfoOutput =
        new Logging::FileOutput(logFile.c_str());
    for (auto l : levels)
    {
      if (l == HDD::Logger::Level::info)
        processingInfoOutput->subscribe(Logging::_SCInfoChannel);
      else if (l == HDD::Logger::Level::warning)
        processingInfoOutput->subscribe(Logging::_SCWarningChannel);
      else if (l == HDD::Logger::Level::error)
        processingInfoOutput->subscribe(Logging::_SCErrorChannel);
      else if (l == HDD::Logger::Level::debug)
        processingInfoOutput->subscribe(Logging::_SCDebugChannel);
    }
    return processingInfoOutput;
  };

  auto destroyFileLogger = [](void *handle) {
    Logging::FileOutput *processingInfoOutput =
        static_cast<Logging::FileOutput *>(handle);
    delete processingInfoOutput;
  };

  HDD::Logger::registerLoggers(error, warning, info, debug, createFileLogger,
                               destroyFileLogger);
}
} // namespace SCAdapter
} // namespace HDD

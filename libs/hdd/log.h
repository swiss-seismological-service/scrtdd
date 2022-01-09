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

#ifndef __HDD_LOG_H__
#define __HDD_LOG_H__

#include "utils.h"

#include <initializer_list>
#include <string>
#include <vector>

namespace HDD {

class Logger
{
public:
  static Logger &getInstance()
  {
    static Logger instance; // Guaranteed to be destroyed.
                            // Instantiated on first use.
    return instance;
  }
  Logger(const Logger &) = delete;
  void operator=(const Logger &) = delete;

  enum class Level
  {
    debug,
    info,
    warning,
    error,
  };

  template <typename... Args> void log(Level l, Args &&... args)
  {
    _log(l, strf(std::forward<Args>(args)...));
  }

  void logToFile(const std::string &logFile, const std::vector<Level> &levels);

private:
  Logger() = default;
  void _log(Level, const std::string &);
};

template <typename... Args> void logDebug(Args &&... args)
{
  Logger::getInstance().log(Logger::Level::debug, std::forward<Args>(args)...);
}

template <typename... Args> void logInfo(Args &&... args)
{
  Logger::getInstance().log(Logger::Level::info, std::forward<Args>(args)...);
}

template <typename... Args> void logWarning(Args &&... args)
{
  Logger::getInstance().log(Logger::Level::warning,
                            std::forward<Args>(args)...);
}

template <typename... Args> void logError(Args &&... args)
{
  Logger::getInstance().log(Logger::Level::error, std::forward<Args>(args)...);
}

} // namespace HDD

#endif

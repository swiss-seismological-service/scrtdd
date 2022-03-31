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

#include <functional>
#include <initializer_list>
#include <string>
#include <vector>

namespace HDD {

class Logger
{
public:
  enum class Level
  {
    debug,
    info,
    warning,
    error,
  };

  Logger()               = delete;
  Logger(const Logger &) = delete;
  void operator=(const Logger &) = delete;

  static void registerLoggers(
      std::function<void(std::string)> error,
      std::function<void(std::string)> warning,
      std::function<void(std::string)> info,
      std::function<void(std::string)> debug,
      std::function<void *(const std::string &, const std::vector<Level> &)>
          createFileLogger,
      std::function<void(void *)> destroyFileLogger)
  {
    _error             = error;
    _warning           = warning;
    _info              = info;
    _debug             = debug;
    _createFileLogger  = createFileLogger;
    _destroyFileLogger = destroyFileLogger;
  }

  static void logError(const std::string &s) { _error(s); }
  static void logWarning(const std::string &s) { _warning(s); }
  static void logInfo(const std::string &s) { _info(s); }
  static void logDebug(const std::string &s) { _debug(s); }

  static void log(const Level l, const std::string &s)
  {
    if (l == Level::error)
      logError(s);
    else if (l == Level::warning)
      logWarning(s);
    else if (l == Level::info)
      logInfo(s);
    else if (l == Level::debug)
      logDebug(s);
  }

  class File
  {
  public:
    File(const std::function<void(void)> &cleanup) : _cleanup(cleanup) {}
    ~File() { _cleanup(); }
    File(const File &) = delete;
    void operator=(const File &) = delete;

  private:
    const std::function<void(void)> _cleanup;
  };

  static std::unique_ptr<Logger::File> toFile(const std::string &filename,
                                              const std::vector<Level> &levels)
  {
    void *handle = _createFileLogger(filename, levels);
    auto cleanup = [handle]() { _destroyFileLogger(handle); };
    return std::unique_ptr<Logger::File>(new File(cleanup));
  }

private:
  static std::function<void(const std::string &)> _error;
  static std::function<void(const std::string &)> _warning;
  static std::function<void(const std::string &)> _info;
  static std::function<void(const std::string &)> _debug;
  static std::function<void *(const std::string &, const std::vector<Level> &)>
      _createFileLogger;
  static std::function<void(void *)> _destroyFileLogger;
};

inline void log(const Logger::Level l, const char *s) { Logger::log(l, s); }
inline void log(const Logger::Level l, const std::string &s)
{
  Logger::log(l, s);
}
inline void log(const Logger::Level l, std::string &&s) { Logger::log(l, s); }

template <typename... Args> void log(const Logger::Level l, Args &&... args)
{
  Logger::log(l, strf(std::forward<Args>(args)...));
}

inline void logDebug(const char *s) { Logger::logDebug(s); }
inline void logDebug(const std::string &s) { Logger::logDebug(s); }
inline void logDebug(std::string &&s) { Logger::logDebug(s); }

template <typename... Args> void logDebug(Args &&... args)
{
  Logger::logDebug(strf(std::forward<Args>(args)...));
}

inline void logInfo(const char *s) { Logger::logInfo(s); }
inline void logInfo(const std::string &s) { Logger::logInfo(s); }
inline void logInfo(std::string &&s) { Logger::logInfo(s); }

template <typename... Args> void logInfo(Args &&... args)
{
  Logger::logInfo(strf(std::forward<Args>(args)...));
}

inline void logWarning(const char *s) { Logger::logWarning(s); }
inline void logWarning(const std::string &s) { Logger::logWarning(s); }
inline void logWarning(std::string &&s) { Logger::logWarning(s); }

template <typename... Args> void logWarning(Args &&... args)
{
  Logger::logWarning(strf(std::forward<Args>(args)...));
}

inline void logError(const char *s) { Logger::logError(s); }
inline void logError(const std::string &s) { Logger::logError(s); }
inline void logError(std::string &&s) { Logger::logError(s); }

template <typename... Args> void logError(Args &&... args)
{
  Logger::logError(strf(std::forward<Args>(args)...));
}

/*
 * Use these macros to let the compiler checks the correctness of
 * the arguments with respect to the format specifiers

#define logLevel(level, fmt, ...)   Logger::log(level, strf(fmt,##__VA_ARGS__))

#define logDebug(fmt, ...)   Logger::logDebug(strf(fmt,##__VA_ARGS__))

#define logInfo(fmt, ...)    Logger::logInfo(strf(fmt,##__VA_ARGS__))

#define logWarning(fmt, ...) Logger::logWarning(strf(fmt,##__VA_ARGS__))

#define logError(fmt, ...)   Logger::logError(strf(fmt,##__VA_ARGS__))

*/

} // namespace HDD

#endif

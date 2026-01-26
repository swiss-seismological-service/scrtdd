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

#include "log.h"
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>

using namespace std;

namespace {

using namespace HDD::Logger;

std::string getCurrentTimestamp()
{
  auto now               = std::chrono::system_clock::now();
  std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);

  std::tm tm;
  localtime_r(&now_time_t, &tm); // Thread-safe on POSIX

  std::stringstream ss;
  ss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
  return ss.str();
}

string formatMessage(const Level level, const std::string &message)
{
  std::string timestamp = getCurrentTimestamp();
  switch (level)
  {
  case Level::debug: return timestamp + " [DEBUG] " + message;
  case Level::info: return timestamp + " [INFO ] " + message;
  case Level::warning: return timestamp + " [WARN ] " + message;
  case Level::error: return timestamp + " [ERROR] " + message;
  default: return timestamp + " [     ] " + message;
  }
}

struct FileLogger
{
  std::ofstream fileStream;
  Level minLevel;
};

std::mutex _logMutex;
std::map<std::string, FileLogger> _fileLoggers;

void logToFile(Level level, const std::string &message)
{
  std::lock_guard<std::mutex> lock(_logMutex);

  for (auto &kv : _fileLoggers)
  {
    FileLogger &logger = kv.second;
    if (level >= logger.minLevel)
    {
      logger.fileStream << formatMessage(level, message) << std::endl;
    }
  }
}

LogCallback _userLogger = nullptr;

Level _minLevel = Level::debug;

} // namespace

namespace HDD {

namespace Logger {

void setLevel(Level l) { _minLevel = l; }

Level getLevel() { return _minLevel; }

void setLogger(LogCallback callback) { _userLogger = std::move(callback); }

void clearLogger() { _userLogger = nullptr; }

void log(const Level level, const std::string &message)
{
  if (level >= _minLevel)
  {
    if (_userLogger)
    {
      _userLogger(level, message);
    }
    else // Fallback to cerr if no user logger is set
    {
      std::cerr << formatMessage(level, message) << std::endl;
    }
  }

  // Always log to active file loggers
  logToFile(level, message);
}

void addFileLogger(const std::string &filename, Level minLevel)
{
  std::lock_guard<std::mutex> lock(_logMutex);

  if (_fileLoggers.count(filename))
  {
    // Logger for this file already exists
    return;
  }

  FileLogger newLogger;
  newLogger.fileStream.open(filename, std::ios_base::app);
  if (!newLogger.fileStream.is_open())
  {
    // Handle error: could not open file
    throw Exception("Could not open log file: " + filename);
  }
  newLogger.minLevel = minLevel;

  _fileLoggers[filename] = std::move(newLogger);
}

void removeFileLogger(const std::string &filename)
{
  std::lock_guard<std::mutex> lock(_logMutex);

  if (_fileLoggers.count(filename))
  {
    _fileLoggers.at(filename).fileStream.close();
    _fileLoggers.erase(filename);
  }
}

} // namespace Logger

} // namespace HDD

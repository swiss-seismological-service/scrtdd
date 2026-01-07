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

#ifndef __HDD_LOG_H__
#define __HDD_LOG_H__

#include "utils.h"

#include <cstdarg>
#include <functional>
#include <initializer_list>
#include <string>
#include <vector>

namespace HDD {

namespace Logger {

enum class Level
{
  debug,
  info,
  warning,
  error
};

using LogCallback =
    std::function<void(Level level, const std::string &message)>;

void setLogger(LogCallback callback);
void log(const Level level, const std::string &message);

void addFileLogger(const std::string &filename, Level minLevel);
void removeFileLogger(const std::string &filename);

inline void logDebug(const std::string &s) { log(Level::debug, s); }
inline void logDebug(const char *s) { log(Level::debug, s); }
inline void logDebug(std::string &&s) { log(Level::debug, s); }
__attribute__((format(printf, 1, 2))) inline void logDebugF(const char *fmt,
                                                            ...)
{
  va_list ap;
  va_start(ap, fmt);
  log(Level::debug, strf(fmt, ap));
  va_end(ap);
}

inline void logInfo(const std::string &s) { log(Level::info, s); }
inline void logInfo(const char *s) { log(Level::info, s); }
inline void logInfo(std::string &&s) { log(Level::info, s); }
__attribute__((format(printf, 1, 2))) inline void logInfoF(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  log(Level::info, strf(fmt, ap));
  va_end(ap);
}

inline void logWarning(const std::string &s) { log(Level::warning, s); }
inline void logWarning(const char *s) { log(Level::warning, s); }
inline void logWarning(std::string &&s) { log(Level::warning, s); }
__attribute__((format(printf, 1, 2))) inline void logWarningF(const char *fmt,
                                                              ...)
{
  va_list ap;
  va_start(ap, fmt);
  log(Level::warning, strf(fmt, ap));
  va_end(ap);
}

inline void logError(const std::string &s) { log(Level::error, s); }
inline void logError(const char *s) { log(Level::error, s); }
inline void logError(std::string &&s) { log(Level::error, s); }
__attribute__((format(printf, 1, 2))) inline void logErrorF(const char *fmt,
                                                            ...)
{
  va_list ap;
  va_start(ap, fmt);
  log(Level::error, strf(fmt, ap));
  va_end(ap);
}

} // namespace Logger

} // namespace HDD

#endif

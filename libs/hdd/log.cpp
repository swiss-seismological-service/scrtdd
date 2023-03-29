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

using namespace std;

namespace HDD {

//
// Register no-op handlers: it is up to the user to register something with
// Logger::registerLoggers(...)
//
std::function<void(const std::string &)> Logger::_error =
    [](const std::string &) {};
std::function<void(const std::string &)> Logger::_warning =
    [](const std::string &) {};
std::function<void(const std::string &)> Logger::_info =
    [](const std::string &) {};
std::function<void(const std::string &)> Logger::_debug =
    [](const std::string &) {};
std::function<void *(const std::string &, const std::vector<Logger::Level> &)>
    Logger::_createFileLogger =
        [](const std::string &, const std::vector<Logger::Level> &) {
          return nullptr;
        };
std::function<void(void *)> Logger::_destroyFileLogger = [](void *) {};

} // namespace HDD

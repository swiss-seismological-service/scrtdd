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

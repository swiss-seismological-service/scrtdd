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

#include "rtddtool.h"

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp/logging/log.h>

int main(int argc, char **argv)
{
  int retCode = EXIT_SUCCESS;

  // Create an own block to make sure the application object
  // is destroyed when printing the overall objectcount.
  {
    Seiscomp::RTDD app(argc, argv);
    retCode = app.exec();
  }

  SEISCOMP_DEBUG("EXIT(%d), remaining objects: %d", retCode,
                 Seiscomp::Core::BaseObject::ObjectCount());

  return retCode;
}

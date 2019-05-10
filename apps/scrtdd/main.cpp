/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 *   You can redistribute and/or modify this program under the             *
 *   terms of the "SED Public License for Seiscomp Contributions"          *
 *                                                                         *
 *   You should have received a copy of the "SED Public License for        *
 *   Seiscomp Contributions" with this. If not, you can find it at         *
 *   http://www.seismo.ethz.ch/static/seiscomp_contrib/license.txt         *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   "SED Public License for Seiscomp Contributions" for more details      *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/



#include "rtdd.h"

int main(int argc, char **argv) {
    int retCode = EXIT_SUCCESS;

    // Create an own block to make sure the application object
    // is destroyed when printing the overall objectcount
    {
        Seiscomp::RTDD app(argc, argv);
        retCode = app.exec();
    }

    SEISCOMP_DEBUG("EXIT(%d), remaining objects: %d",
                   retCode, Seiscomp::Core::BaseObject::ObjectCount());

    return retCode;
}

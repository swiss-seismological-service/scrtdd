<pre>
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
 *   Developed by Luca Scarabello, Tobias Diehl                            *
 *                                                                         *
 ***************************************************************************/
</pre>

[![DOI](https://zenodo.org/badge/246001157.svg)](https://zenodo.org/badge/latestdoi/246001157)

# Description

rtDD is a [SeisComP](<https://github.com/SeisComP>) extension module that implements Double-Difference event relocation both in Real-Time, one event at the time, and classic offline mode, where an earthquake catalog is relocated as a whole.

rtDD contains a C++ library for double-difference inversion which doesn't depend on SeisComP and can be used in other systems.

The actual methods are based on the paper "Near-Real-Time Double-Difference Event Location Using Long-Term Seismic Archives, with Application to Northern California" by Felix Waldhauser and "A Double-Difference Earthquake Location Algorithm: Method and Application to the Northern Hayward Fault, California" by Waldhauser & Ellsworth.

The double-difference equation system solver uses [LSQR by Chris Paige, Michael Saunders](<https://web.stanford.edu/group/SOL/software/lsqr/>) and [LSMR by David Fong, Michael Saunders](<https://web.stanford.edu/group/SOL/software/lsmr/>) algorithms from [this](<https://github.com/tvercaut/LSQR-cpp/>) Apache Licensed beautiful implementation by Tom Vercautereen.

rtDD also supports [NonLinLoc by Anthony Lomax](<http://alomax.free.fr/nlloc/>) grid file format alongside the travel time formats natively supported by SeisComP (LOCSAT and libtau). A feature that enables 3D velocity models within the tool.

rtDD acknowledges support from [Geothermica](http://www.geothermica.eu/)

## Documentation

You can find the online documentation at https://docs.gempa.de/scrtdd/current/

## Installation from binaries

You can find compiled version of this module at https://data.gempa.de/packages/Public/. The installation file is a compressed tar archive containing the binary distribution of this module, which can be extracted in your SeisComP installation folder.

## Installation from source code

In order to use this module the sources have to be merged into the *SeisComP* sources, then *SeisComP* can be compiled and installed as usual. Follow the up-to-date procedure at https://github.com/SeisComP/seiscomp

<pre>
#
# get SeisComP (follow official documentation)
#
git clone https://github.com/SeisComP/seiscomp.git seiscomp
cd seiscomp/src/base
git clone https://github.com/SeisComP/common.git
git clone https://github.com/SeisComP/main.git
[...etc...]

#
# Add rtdd into SeisComP
#
cd seiscomp/src/extras
git clone https://github.com/swiss-seismological-service/scrtdd.git
cd scrtdd

#
# Compile code as per SeisComP official documentation
#
...
</pre>

To update rtDD to a new version:

<pre>
#
# go to rtDD folder
#
cd seiscomp/src/extras/scrtdd

#
# fetch the changes that happened on rtDD repository
#
git fetch origin -p --tags

#
# checkout the new version
#
git checkout master
git rebase origin/master
</pre>



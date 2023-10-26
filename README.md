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

Please cite the code as:

"Luca Scarabello & Tobias Diehl (2021). swiss-seismological-service/scrtdd. Zenodo doi: 10.5281/zenodo.5337361"

# Description

rtDD is a [SeisComP](<https://github.com/SeisComP>) extension module that implements
Double-Difference event relocation both in Real-Time, one event at the time, and classic
offline mode, where an earthquake catalog is relocated as a whole.

rtDD contains a C++ library for double-difference inversion which doesn't depend on
SeisComP and can be used in other systems. See for example [pyrtdd](https://github.com/swiss-seismological-service/pyrtdd).

rtDD acknowledges support from [Geothermica](http://www.geothermica.eu/)

## Documentation

You can find the online documentation at https://docs.gempa.de/scrtdd/current/

## Installation from binaries

You can find compiled version of this module at https://data.gempa.de/packages/Public/.
The installation file is a compressed tar archive containing the binary distribution of
this module, which can be extracted in your SeisComP installation folder.

## Installation from source code

In order to use this module the sources have to be merged into the *SeisComP* sources,
then *SeisComP* can be compiled and installed as usual. Follow the up-to-date procedure
at https://github.com/SeisComP/seiscomp

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

## Tests

If tests were enabled during the compilation of SeiComP then `scrtdd` can be tested
running the following command in the build folder:
.
<pre>
seiscomp/build$ ctest -R test_hdd
</pre>

Alternatively `scrtdd` can be tested together with the whole SeisComP test suite with:

<pre>
seiscomp/build$ make test
</pre>

See SeisComP's [unit testing guide](https://docs.gempa.de/seiscomp/current/base/tests.html).

Please note that due to the size of the NonLinLoc grids, not all the transformation are
tested. To let the tests cover all the implemented transformations the script 
`libs/hdd/test/data/nll/generate.ch` should be run, which will generate the missing grids.

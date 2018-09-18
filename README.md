<pre>
/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 *   You can redistribute and/or modify this program under the             *
 *   terms of the SeisComP Public License.                                 *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   SeisComP Public License for more details.                             *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/
</pre>
**Description**

This module implements a Real Time Double-Difference Event Location method as described in the paper "Near-Real-Time Double-Difference Event Location Using Long-Term Seismic Archives, with Application to Northern California" by Felix Waldhauser.

This module can also be used to perform (non-real-time) full event catalog Double-Difference relocation, which is covered by another paper: "A Double-Difference Earthquake Location Algorithm: Method and Application to the Northern Hayward Fault, California" by Felix Waldhauser et al.

The actual Double-Difference inversion method is performed by the hypoDD software, which has to be installed in the system (version 2.1b).

**Compile**

In order to use this module the sources have to be compiled to an executable. Merge them into the Seiscomp3 sources and compile Seiscomp3 as usual.
<pre>
# merge rtdd-addons and Seiscomp3
git clone https://github.com/SeisComP3/seiscomp3.git sc3-src
cd sc3-src
git submodule add -f https://github.com/luca-s/rtdd-addons.git src/rtdd
</pre>
For compiling Seiscomp3, please refer to https://github.com/SeisComP3/seiscomp3#compiling

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
# Description

The SCRTDD is a [Seiscomp3](<https://github.com/SeisComP3/seiscomp3>) extension module that implements both Real-Time Double-Difference Event relocation and classic offline multi-event (catalog) Double-Difference relocation.

The actual methods are based on the paper "Near-Real-Time Double-Difference Event Location Using Long-Term Seismic Archives, with Application to Northern California" by Felix Waldhauser and "A Double-Difference Earthquake Location Algorithm: Method and Application to the Northern Hayward Fault, California" by Waldhauser & Ellsworth.

The double-difference equation system solver uses [LSQR by Chris Paige, Michael Saunders](<https://web.stanford.edu/group/SOL/software/lsqr/>) and [LSMR by David Fong, Michael Saunders](<https://web.stanford.edu/group/SOL/software/lsmr/>) algorithms from [this](<https://github.com/tvercaut/LSQR-cpp/>) Apache Licencesed beautiful implementation by Tom Vercautereen.

## Compile

In order to use this module the sources have to be compiled to an executable. Merge them into the Seiscomp3 sources and compile Seiscomp3 as usual.
<pre>
# merge rtdd-addons and Seiscomp3
git clone https://github.com/SeisComP3/seiscomp3.git sc3-src
cd sc3-src
git submodule add -f https://github.com/swiss-seismological-service/scrtdd.git src/rtdd-addons
</pre>
For compiling Seiscomp3, please refer to https://github.com/SeisComP3/seiscomp3#compiling


# Getting Started

## 1. Defining a background catalog for real time relocations

New origins will be relocated in real time against a background catalog of high quality locations. Those high quality events that form the background catalog can be already present in seiscompo database or not. In the latter case the background catalog has to be generated. Either way, the catalog has to be specified in scrtdd when configuring it.

### 1.1. Using origins from seiscomp database as background catalog

If we already have a high quality catalog, then we can easily specify it in scrtdd configuration as a path to a file containing a list of origin id. The file format must be a csv file in which there must be at least one column called seiscompId, from which the ids will be fetched by scrtdd.

E.g. *file myCatalog.csv*

```
seiscompId
Origin/20181214107387.056851.253104
Origin/20180053105627.031726.697885
Origin/20190121103332.075405.6234534
Origin/20190223103327.031726.346363
[...]
```

![Catalog selection option](/data/catalog-selection1.png?raw=true "Catalog selection from event/origin ids")

### 1.2. Generating a background catalog

In the more general case in which we don't already have a background catalog, we have to build one using scrtdd. Generating a background catalog involves two steps: first we select the candidate events that might have low quality locations (we still want them to be as accurate as possible though, and with many phases. That is we want origins manully relocated using manual picks) and then we use scrtdd to relocated those events to achieve the quality needed for a background catalog. The higher the quality of the background catalog, the higher will be the quality of the real time origin relocations.

#### 1.2.1 Dumping the candidate events to files

Once the candidate origins have been selected, we write their seiscomp id in a file using the same format as the example above (myCatalog.csv). Once that's done, we run the command:

```
scrtdd --dump-catalog myCatalog.csv [db and log options]
```

The above command will generate three files (event.csv, phase.csv and stations.csv) which contain the information needed by scrtdd to relocate the selected events. 

E.g. *file event.csv*

```
id,isotime,latitude,longitude,depth,magnitude,horizontal_err,vertical_err,rms
1,2014-01-10T04:46:47.689331Z,46.262846,7.400132,8.6855,1.63,0.0000,0.0000,0.1815
2,2014-01-19T05:24:26.754208Z,46.264482,7.404143,8.4316,0.94,0.0000,0.0000,0.1740
3,2014-02-21T04:05:27.03289Z,46.266118,7.402066,7.3145,0.37,0.0000,0.0000,0.1177
4,2014-04-02T17:05:28.141739Z,46.262846,7.408248,7.0098,0.42,0.0000,0.0000,0.1319
```

E.g. *file station.csv*

```
id,latitude,longitude,elevation,networkCode,stationCode,locationCode
4DAG01,46.457412,8.079460,2358.0,4D,AG01,
4DAG02,46.460620,8.078122,2375.0,4D,AG02,
4DAG03,46.458288,8.075408,2369.0,4D,AG03,00
4DBSG1,46.107760,7.732020,3378.0,4D,BSG1,AB
```

E.g. *file phase.csv*

```
eventId,stationId,isotime,lowerUncertainty,upperUncertainty,type,networkCode,stationCode,locationCode,channelCode,evalMode
1,CHSIMPL,2014-01-10T04:47:02.000765Z,0.100,0.100,Sg,CH,SIMPL,,HHT,manual
1,CHNALPS,2014-01-10T04:47:06.78218Z,0.100,0.100,P1,CH,NALPS,,HHR,manual
1,CHBNALP,2014-01-10T04:47:05.918759Z,0.200,0.200,P1,CH,BNALP,,HHZ,automatic
1,CHFUSIO,2014-01-10T04:47:04.812236Z,0.100,0.100,Pg,CH,FUSIO,,HHR,manual
1,FRRSL,2014-01-10T04:47:13.089093Z,0.200,0.200,Sg,FR,RSL,00,HHT,manual
1,FRRSL,2014-01-10T04:47:02.689842Z,0.050,0.050,Pg,FR,RSL,00,HHZ,automatic
1,CHGRIMS,2014-01-10T04:47:01.597023Z,0.100,0.100,Pg,CH,GRIMS,,HHR,manual
1,IVMRGE,2014-01-10T04:46:58.219541Z,0.100,0.100,Pg,IV,MRGE,,HHR,manual
```

Now that we have dumped the events (event.csv, phase.csv, stations.csv) we might perform some editing of those files, if required, then we relocate them (next paragraph).
 
See also the `Troubleshooting` paragraph for details on how to specify a database connection via command line and for printing the logs on the console for easy debugging.

Also it is useful the `--dump-catalog-options` option, which allows to use event id instead of origin id in myCatalog.csv. For example we can list all events between two dates and then ask scrtdd to extract the preferred manual origins within our profile region (our catalog):

```
# prepare myCatalog.csv
echo seiscompId > myCatalog.csv
# save all events ids between 2018-11-27 and 2018-12-14 in myCatalog.csv
scevtls [db and log options] --begin "2018-11-27 00:00:00" --end "2018-12-14 00:00:00" >> myCatalog.csv
# create the catalog using only manually reviewed preferred origins within the region defined in myProfile
scrtdd --dump-catalog myCatalog.csv --dump-catalog-options 'preferred,onlyManual,any,none,myProfile' [db and log options]
```

Note: event.csv, phase.csv and stations.csv files are the extended catalog format, useful when the catalog information is fully contained in those files. Compare to the paragraph 2.1 example format (list of origin ids) the extended format doesn't require the catalog to be present on the seiscomp database. This might come in handy when using events coming from a different seiscomp installation e.g the events can be dumped (--dump-catalog) and used in another machine. Also those file can be easily generated from another source.


#### 1.2.2 Relocating the candidate events (multi-event mode) to generate a high quality background catalog

To relocate a catalog (multi-event relocation) we need to create a new profile in scrtdd configuration and then we have to set the generated files (event.csv, phase.csv, stations.csv) as the catalog of the profile.

![Catalog selection option](/data/catalog-selection2.png?raw=true "Catalog selection from raw file format")

At this point we have to configure the other profile options that control the relocation process: `doubleDifferenceObservations`, which control the creation of catalog absolute travel time entries (dt.ct file in HypoDD terminology) and cross correlated differential travel times for pairs of events (dt.cc file in HypoDD terminology). The clustering options should work just fine with the default values, however some adjustments are worth it. For example we force to have only well connected events (`minNumNeigh` and `minObservationPerEvPair`), which helps in avoiding ill-defined double difference system that can be hard to solve. Then we minimize `maxNumNeigh` to reduce computation time (we can get very good results without using all the possible neighbours for every event). Finallly we like to disable the ellipsoid algorithms (`numEllipsoids=0`) since that is mostly useful in single event relocation.

![Relocation options](/data/multiEventStep2options.png?raw=true "Relocation options")

Then it is time to set the cross-correlation parameters, which require a more careful selection and it is covered in the next paragraph.

In the ```data``` folder of this project there are some hypoDD 1.3 and 2.1b configuration example files that have to be adapted to reflect your specific use case, especially the velocity model. We do not reccomend using those example files blindly, a good understanding of HypoDD software is required to obtain good results.

Finally, when the configuration is ready, we can relocate the catalog with the command:

```
scrtdd --reloc-profile profileName [db and log options] 
```

scrtdd will relocated the catalog and will generate another set of files reloc-event.csv reloc-phase.csv and reloc-stations.csv, which together define a new catalog with relocated origins. 

At this point we should check the relocated events (and logs) and see whether the results make sense and are satisfactory. Usually we want to keep tuning the scrtdd settings and relocate the catalog multiple times until we are sure we reached the best relocations. Having a good background catalog is the base for good real-time relocations.

As an example you can see below two catalogs before and after scrtdd relocation:


![Relocation example picture](/data/multiEventRelocationExample.png?raw=true "Relocation example")

The relocation output files (reloc-event.csv reloc-phase.csv and reloc-stations.csv) will become the background catalog used in real-time relocation and this is the only output we need to keep from the relocation process. The profile configuration can be now deleted or, in the case we want to kept it, it has to be removed from the list of active profiles (`scrtdd.activeProfiles`) to avoid interaction with real-time processing.

Now that we have a high quality background catalog we are ready to perform real-time relocation. For that we will create a new profile whose catalog will be the relocation output files:

![Catalog selection option](/data/catalog-selection3.png?raw=true "Catalog selection from raw file format")


Notes:

- It is possible to transform the catalog in xml (see `--dump-catalog-xml` option) and insert the origins in the seiscomp database and configure the background catalog for real-time as a list of seiscomp origin ids. While it is neat to have the background catalog in the seiscomp database, this approach is inconvenient in the common scenario where the background catalog is re-generated from time to time to include more recent seismicity.

- Be aware that the first time you run the command it will be very slow because the waveforms have to be downloaded from the configured recordstream and saved to disk for the next runs, which will be much faster. 


### 1.2.3 Useful options

There are also some other interesting catalong related options:


```
scrtdd --help

MultiEvents:
  --reloc-profile arg                   Relocate the catalog of profile passed 
                                        as argument 
Catalog:
  --dump-catalog arg                    Dump the seiscomp event/origin id file 
                                        passed as argument into a catalog file 
                                        triplet (station.csv,event.csv,phase.cs
                                        v)
  --dump-catalog-options arg            Allows --dump-catalog to accept event 
                                        ids besides origin ids. For each event 
                                        id an origin will be selected following
                                        the provided options whose format is: 
                                        'type,evalmode,includeCreator,excludeCr
                                        eator,region', where 
                                        type=preferred|last|first  
                                        evalmode=any|onlyManual|onlyAutomatic  
                                        includeCreator=any|author|methodID  
                                        excludeCreator=none|author|methodID  
                                        region=any|profileName e.g. to select 
                                        preferred origins within my profile 
                                        region given the input event ids use 
                                        'preferred,any,any,none,myProfile 
  --dump-catalog-xml arg                Convert the input catalog into XML 
                                        format. The input can be a single file 
                                        (containing seiscomp origin ids) or a 
                                        catalog file triplet 
                                        (station.csv,event.csv,phase.csv)
  --merge-catalogs arg                  Merge in a single catalog all the 
                                        catalog file triplets 
                                        (station1.csv,event1.csv,phase1.csv,sta
                                        tion2.csv,event2.csv,phase2.csv,...) 
                                        passed as arguments
  --merge-catalogs-keepid arg           Similar to --merge-catalogs option but 
                                        events keeps their ids. If multiple 
                                        events share the same id, subsequent 
                                        events will be discarded.

```

## 2. Cross-correlation, waveform filtering and signal to noise ratio options

Good cross-correlation results are needed to achieve high resolution travel time differences between event phases at common stations. Those high resolution differences result in high quality relocations.

Also, good cross-correlation results are very useful since scrtdd uses cross-correlation to detect phases for station that have no associated picks. This is especially useful when relocating automatic origins and only few picks are available.

This set of configuration parameters require some trial and error to be correctly set, but once they are tuned they can be kept identical for both multi event and single event relocation without any major changes, except for the parameter ```maxDelay``` that should be increased in real-time relocation to accomodate the larger uncertainty of automatic picks (unless scrtdd is used to only relocate manually reviewed origins). Here is an example configuration:

![Relocation options](/data/xcorr.png?raw=true "Relocation options")


### 2.1 Logs

To help figuring out the right values for cross-correlation the log file (or console output with ```--console=1 --verbosity=3```) come in handy:

```
[info] Cross-correlation statistics: performed 40361, waveforms with Signal to Noise ratio too low 2435, waveforms not available 98
[info] Total xcorr 40361 (P 59%, S 41%) success 28% (11499/40361). Successful P 22% (5300/23844). Successful S 38% (6199/16517)
[info] xcorr on actual picks 24784/40361 (P 60%, S 40%) success 37% (9186/24784). Successful P 31% (4629/14761). Successful S 45% (4557/10023)
[info] xcorr on theoretical picks 15577/40361 (P 58%, S 42%) success 15% (2313/15577). Successful P 7% (671/9083). Successful S 25% (1642/6494)
```

We can see that the statistics are broken down in actual picks and theoretical picks. This is because scrtdd computes theoretical picks so that they can be used in cross-correlation. If the resulting correlation coefficient is good enough the result is used in the double-difference system inversion. This is especially useful to increase the number of observations in real-time when automatic origins have only few automatic picks/phases. A new pick is also created and added to the relocated origin when the cross-correlation results are good enough.


### 2.2 Eval-xcorr command

A more sophisticated method for evaluating the settings is the `--eval-xcorr` command (here we use `--verbosity=2` because the statistics are printed at this log level, useful to skip other information):

```
scrtdd --eval-xcorr profileName --verbosity=2 --console=1
```

Example output:

```
[...]
[warning] <<<Progressive stats>>>
[...]
[warning] <<<Final stats>>>
Cumulative stats: #pha   8882 pha good CC  73% avg coeff 0.81 avg #matches/ph 13
Stats by inter-event distance in 0.10 km step
Inter-event dist 0.00-0.10 [km]: #CC  83750 good CC  28% avg coeff 0.85 avg #matches/ev 10
Inter-event dist 0.10-0.20 [km]: #CC 110997 good CC  24% avg coeff 0.83 avg #matches/ev  9
Inter-event dist 0.20-0.30 [km]: #CC  47531 good CC  21% avg coeff 0.82 avg #matches/ev  8
Inter-event dist 0.30-0.40 [km]: #CC  29676 good CC  22% avg coeff 0.82 avg #matches/ev  8
Inter-event dist 0.40-0.50 [km]: #CC  20583 good CC  16% avg coeff 0.82 avg #matches/ev  6
Inter-event dist 0.50-0.60 [km]: #CC  26951 good CC  15% avg coeff 0.81 avg #matches/ev  6
Inter-event dist 0.60-0.70 [km]: #CC  32748 good CC  17% avg coeff 0.82 avg #matches/ev  6
Inter-event dist 0.70-0.80 [km]: #CC  10117 good CC  15% avg coeff 0.81 avg #matches/ev  5
Inter-event dist 0.80-0.90 [km]: #CC   2949 good CC  14% avg coeff 0.82 avg #matches/ev  4
Inter-event dist 0.90-1.00 [km]: #CC   1409 good CC  14% avg coeff 0.82 avg #matches/ev  3
Inter-event dist 1.00-1.10 [km]: #CC   7467 good CC  11% avg coeff 0.81 avg #matches/ev  3
[...]
Stats by event to station distance in 3.00 km step
Station dist   3-6   [km]: #pha     40 pha good CC  95% avg coeff 0.81 avg #matches/ph 29
Station dist   6-9   [km]: #pha   1704 pha good CC  84% avg coeff 0.81 avg #matches/ph 24
Station dist   9-12  [km]: #pha   1040 pha good CC  68% avg coeff 0.81 avg #matches/ph  9
Station dist  12-15  [km]: #pha   1278 pha good CC  76% avg coeff 0.80 avg #matches/ph  8
Station dist  15-18  [km]: #pha    412 pha good CC  64% avg coeff 0.81 avg #matches/ph  5
Station dist  18-21  [km]: #pha    438 pha good CC  70% avg coeff 0.81 avg #matches/ph 14
Station dist  21-24  [km]: #pha    612 pha good CC  77% avg coeff 0.81 avg #matches/ph 11
Station dist  24-27  [km]: #pha    552 pha good CC  68% avg coeff 0.80 avg #matches/ph  7
Station dist  27-30  [km]: #pha    713 pha good CC  74% avg coeff 0.81 avg #matches/ph 11
Station dist  30-33  [km]: #pha    185 pha good CC  61% avg coeff 0.80 avg #matches/ph  3
Station dist  33-36  [km]: #pha    858 pha good CC  79% avg coeff 0.83 avg #matches/ph 15
Station dist  36-39  [km]: #pha    217 pha good CC  62% avg coeff 0.81 avg #matches/ph  5
Station dist  39-42  [km]: #pha    126 pha good CC  58% avg coeff 0.82 avg #matches/ph  6
[...]
Stats by station
4D.MH36.A   : #pha     37 pha good CC  81% avg coeff 0.79 avg #matches/ph  4
4D.MH44.A   : #pha     38 pha good CC  92% avg coeff 0.82 avg #matches/ph  5
4D.MH48.A   : #pha     27 pha good CC  85% avg coeff 0.81 avg #matches/ph  4
4D.MH54.A   : #pha     25 pha good CC  92% avg coeff 0.83 avg #matches/ph  2
4D.RA43.    : #pha      5 pha good CC   0% avg coeff 0.00 avg #matches/ph  0
8D.AMIDI.   : #pha     14 pha good CC  93% avg coeff 0.84 avg #matches/ph  3
8D.BTAO.    : #pha      5 pha good CC   0% avg coeff 0.00 avg #matches/ph  0
8D.BTNF.    : #pha      6 pha good CC   0% avg coeff 0.00 avg #matches/ph  0
8D.GSF03.   : #pha      7 pha good CC   0% avg coeff 0.00 avg #matches/ph  0
[...]
```

It is especially interesting to compare the results after the catalog has been relocated and loaded back to the profile. The new statistics should see how number of detected phases, quality of correlation coefficient decrease with increasing inter-event distance (if the relocation is good).

### 2.3 Waveforms inspection

A more in-depth source of information for waveform filtering and signal to noise ratio options comes from this option:


```
scrtdd --help
  --debug-wf                            Enable the saving of waveforms 
                                        (filtered/resampled, SNR rejected, ZRT 
                                        projected and scrtdd detected phase) 
                                        into the profile working directory. 
```

Simply adding ```--debug-wf``` to the command line will make scrtdd dump to disk miniseed files for inspection (e.g. scrttv waveformfile.mseed). Just make sure to delete the folder before using this option to make sure to not look at previous relocation output.  This option can be added to any scrtdd commands (e.g. `--relocate-profile`, `--ev`, `--origin-id` ) but it is mostly useful when relocating a single event mode because in multi-event mode there will be way too many waveforms to be able to check them all manually, although we can still do some random check to get an overall feeling of the filtering and snr.

The generated waveforms will be stored in `workingDirectory/profileName/wfdebug/` (e.g. `~/seiscomp3/var/lib/rtdd/myProfile/wfdebug/`) after filtering and resampling, that is the same wavforms used for the cross-correlation. The files follow the patters:

* evNumber.NET.STATION.phaseTime.manual.mseed       (e.g. ev56.CH.SAYF2.S.manual.mseed)       - event 56 manual S phase on CH.SAYF2 station
* evNumber.NET.STATION.phaseTime.automatic.mseed    (e.g  ev4.CH.SAYF2.P.automatic.mseed)     - event 4 automatic P phase on CH.SAYF2 station
* evNumber.NET.STATION.phaseTime.theoretical.mseed  (e.g. ev267.CH.SFRU.Pt.theoretical.mseed) - event 267 theoretical P phase on CH.SFRU station
* evNumber.NET.STATION.phaseTime.snr-rejected.mseed (e.g. ev9.XY.VET02.S.snr-rejected.mseed)  - event 9 S phase not used in cross-correlation due to the configured signal to noise ratio threshold

Also, scrtdd logs tell us the details of the cross-correlation, so that we can see what waveform was cross-correlated with which others and then inspect the corresponding waveform miniseed files e.g.:

```
12:21:52 [info] xcorr: event   267 sta   CH LKBD2 dist 19.78 [km] - 3 P phases, mean coeff 0.71 lag 0.13 (events: 1 78 265 )
12:21:52 [info] xcorr: event   267 sta   CH VANNI dist 20.05 [km] - 4 P phases, mean coeff 0.75 lag -0.14 (events: 1 48 54 265 )
12:21:52 [info] xcorr: event   267 sta   XY AGA01 dist 20.76 [km] - low corr coeff pairs
12:21:52 [info] xcorr: event   267 sta   8D  RAW4 dist 22.49 [km] - 1 P phases, mean coeff 0.73 lag 0.03 (events: 78 )
12:21:52 [info] xcorr: event   267 sta   XP MONT1 dist 22.63 [km] - low corr coeff pairs
12:21:52 [info] xcorr: event   267 sta   CH GRYON dist 24.60 [km] - 2 S phases, mean coeff 0.76 lag -0.21 (events: 17 21 )
12:21:52 [info] xcorr: event   267 sta   CH SZWD2 dist 25.53 [km] - low corr coeff pairs
12:21:52 [info] xcorr: event   267 sta   CH  SGAK dist 26.03 [km] - low corr coeff pairs
12:21:52 [info] xcorr: event   267 sta   CH  SZIM dist 26.35 [km] - low corr coeff pairs
12:21:52 [info] xcorr: event   267 sta   CH  SCOD dist 27.00 [km] - low corr coeff pairs
12:21:52 [info] xcorr: event   267 sta   CH   DIX dist 27.78 [km] - 1 P phases, mean coeff 0.72 lag 0.04 (events: 265 )
12:21:52 [info] xcorr: event   267 sta   CH FULLY dist 28.33 [km] - 1 P phases, mean coeff 0.71 lag 0.09 (events: 111 )
12:21:52 [info] xcorr: event   267 sta   XY NIE01 dist 29.34 [km] - 1 S phases, mean coeff 0.77 lag -0.34 (events: 38 )
12:21:52 [info] xcorr: event   267 sta   CH LAUCH dist 30.22 [km] - 19 P phases, mean coeff 0.84 lag -0.09 (events: 1 7 9 21 29 48 54 55 58 59 70 77 78 87 111 196 225 234 265 )
12:21:52 [info] xcorr: event   267 sta   CH LAUCH dist 30.22 [km] - 5 S phases, mean coeff 0.86 lag -0.25 (events: 11 29 52 111 234 )
12:21:52 [info] xcorr: event   267 sta   XY RAR01 dist 30.59 [km] - low corr coeff pairs
12:21:52 [info] xcorr: event   267 sta   XY RAR02 dist 32.42 [km] - low corr coeff pairs
12:21:52 [info] xcorr: event   267 sta   CH LAVEY dist 32.83 [km] - 1 P phases, mean coeff 0.71 lag 0.07 (events: 9 )
12:21:52 [info] xcorr: event   267 sta   CH LAVEY dist 32.83 [km] - 1 S phases, mean coeff 0.79 lag 0.07 (events: 225 )

```

For comparison we can always find the raw waveforms (not processed) fetched from the configured recordstream and used as a cache in `workingDirectory/profileName/wfcache/` (e.g. `~/seiscomp3/var/lib/rtdd/myProfile/wfcache/`):
* `NET.ST.LOC.CH.startime-endtime.mseed`




## 3. Real time single origin relocation

We first need to create a new profile or we can modify the profile used to create the background catalog if we don't need that anymore. If we decice to go for a new profile we might want to copy the waveform cache folder of the background catalog profile to the new real-time profile cache folder to avoid re-downloading all the catalog waveforms (see the `Performance` paragraph for more details).

The waveform cache folder is located in `workingDirectory/profileName/wfcache/` (e.g. `~/seiscomp3/var/lib/rtdd/myProfile/wfcache/`).

### 3.1 Testing

Real time relocation uses the same configuration we have seen in full catalog relocation, but real time relocation is done in two steps, each one controlled by a specific configuration.

Step 1: location refinement. In this step scrtdd performs a preliminary relocation of the origin using only catalog absolute travel time entries (dt.ct only).

Step 2: the refined location computerd in the previous step is used as starting location to perform a more precise relocation using both catalog absolute travel times (dt.ct) and differential travel times from cross-correlation (dt.cc). If step1 fails, step2 is attempted anyway.

If step2 completes successfully the relocated origin is sent to the messaging system.

You can see below the most important parameters we want to configure. We enable the ellipsoid clustering algorithms selecting `numEllipsoids != 0` and we define the number of neighbour events to use for relocation. We also increased `maxEllipsoidSize` in `doubleDifferenceObservationsNoXcorr` to accomodate the possible large error in the automatic origin locations.
 

![Relocation options](/data/step1options.png?raw=true "Relocation options")
![Relocation options](/data/step2options.png?raw=true "Relocation options")

You might consider testing the configuration on some existing events to make sure the parameters are suitable for real time relocation. To test the real time relocation there are two command line options which relocate existing origins:

```
scrtdd --help

SingleEvent:
  -O [ --origin-id ] arg                Relocate the origin (or multiple 
                                        comma-separated origins) and send a 
                                        message. Each origin will be processed 
                                        accordingly with the matching profile 
                                        region unless --profile option is used
  --ep arg                              Event parameters XML file for offline 
                                        processing of contained origins (imply 
                                        test option). Each contained origin 
                                        will be processed accordingly with the 
                                        matching profile region unless 
                                        --profile option is used. In 
                                        combination with origin-id option this 
                                        produces an xml output
  --test                                Test mode, no messages are sent
  --profile arg                         Force a specific profile to be used 
                                        when relocating an origin. This 
                                        overrides the selection of profiles 
                                        based on region information and the 
                                        initial origin location
```

If we want to process an origin we can run the following command and then check on scolv the relocated origin (the messaging system must be active). This is mostly useful when we want to relocate an origin on a running system and keep the relocation:

```
scrtdd --origin-id someOriginId [db and log options] 
```

For testing purpose we are more likely interested in not interfere with database and messaging system so we can use the `--ep` option and the relocated origin will be sent to the output as xml file. We can finally open the xml file with scolv for inspection:

```
scrtdd --origin-id someOriginId --ep - [db and log options] >  relocated-origin.xml
```

Alternatively the `--ep` option (without `-O`) can process all origins contained in the input XML file:

```
# dump origin
scxmldump -fPAMF -p -O originId -o origin.xml --verbosity=3  --console=1
# relocate
scrtdd --ep origin.xml [db and log options] > relocated-origin.xml
```

Also, as explained in the cross-correlation settings paragraph, we can use the ```--debug-wf``` option to help debugging.

As an example you can see below the single event relocation of several manually reviewed origins (when relocating automatic origins the quality and number of relocated origins is certainly lower).

![Single event relocation example picture](/data/singleEventRelocationExample.png?raw=true "Single Event Relocation example")

Once we are happy with the configuration we can simply enable and start scrtdd as any other module and it will start relocating origins as soon as they arrive in the messsaging system

### 3.2 Phase update

During real-time relocations, once cross-correlation is completed, scrtdd updates automatic phase pick time and uncertainties of the origin is going to relocate. The pick time is updated accordingly to the average lag detected by all the good (above configured threshold) cross-correlation results. Since the real-time events are cross-correlated against catalog events, which have good manual picks, the new pick time is expected to improve pick accuracy. Also the pick uncertainty is derived from the uncertainties of catalog-events. If there is no cross-correlations result above the configured threshold, the pick time is kept untouched.

The same process is applied to theoretical phases, which are computed by scrtdd for each station that has no pick.Every theoretical phase that has at least one good cross-correlation result is added to the relocated origin, with pick time and uncertainties derived from catalog phases. Theoretical phases that have no good cross-correlation results never end up in the relocated origin.

Picks that have been updated are identifiable by a `x` suffix.

Manual picks are never modified.

### 3.3 RecordStream configuration

SeisComP3 applications access waveform data through the RecordStream interface and it is usually configured in global.cfg, for example:

```
recordstream = combined://slink/localhost:18000;sdsarchive//mnt/miniseed
```

This configuration is a combination of seedlink and sds archive, which is very useful because scrtdd catalog waveforms can be retrieved via sds while real-time event data can be fetched via seedlink (much faster since recent data is probably already in memory).

Please note that seedlink service might delay the relocation of incoming origins due to timeouts. For this reason we suggest to pass to configuration option: timeout and retries

e.g. Below we force a timeout of 1 seconds (default is 5 minutes) and do not try to reconnect (`scrtdd` will deal with what data is available):

```
recordstream = combined://slink/localhost:18000?timeout=1&retries=0;sdsarchive//rz_nas/miniseed
```

Also we use `cron.delayTimes` option to re-perform the relocation some minutes later in case we know more waveforms will become available at a later time.

## 4. Performance and waveforms data

scrtdd spends most of the relocation time downloading waveforms (unless the recordstream points to a local disk storage) and the rest is shared betweeen cross-correlation and double difference system inversion. For this reason the waveforms are cached to disk after being downloaded so that there is no need to download them again. This is obviously true for catalog phases waveforms that are re-used over and over, but that's not true for real-time events waveforms and also for some temporary waveforms (e.g. theoretical phases) that are never saved. The cache folder is `workingDirectory/profileName/wfcache/` (e.g. /installation/path/seiscomp3/var/lib/rtdd/myProfile/wfcache/). 

However, during the parameters tuning phase, the user performs relocations (both single-event and multi-event) from the command line several times to find the best configuration. For those speacial cases even the temporary waveforms are saved to make the process much faster. A different folder is used: `workingDirectory/profileName/tmpcache/` which can be deleted after the parameter tuning phase. The commands that store temporary waveforms are: `--reloc-profile`, `--ep`, `-O`, `--origin-id`.

Due to the time required to download waveforms, when relocating events in real-time, two configuration options have a huge impact on the performance: `doubleDifferenceObservations.clustering` and `crosscorrelation.s-phase.components`. `doubleDifferenceObservations.clustering` is relevant because we can specify how many neighbouring events and how many phases we want to use, which consequently determines the number of cross-correlations to perform, the waveforms to download and the size of the input for the double difference inversion. `crosscorrelation.s-phase.components` defines the components we want to use in the cross-correlation of the `S` phases. Usign the `T` components means scrtdd has to download the additional two components to perform the projection to ZRT, which might be more accurate than 'Z' only.


### 4.1 Catalog waveforms preloading

When scrtdd starts for the first time it loads all the catalog waveforms and store them to disk. In this way the waveforms become quickly available in real-time without the need to access the recordstream. However this takes some time.  If for some reasons (debugging?) we need to force the downloading of all waveforms before starting the module we can use the following option (once the waveforms are downloaded they will not be downloaded again unless the files are deleted from disks):

```
scrtdd --help
  --load-profile-wf arg                 Load catalog waveforms from the 
                                        configured recordstream and save them 
                                        into the profile working directory
``` 


## 5. Locator plugin

A (re)locator plugin is also avaiable in the code, which makes scrtdd available via scolv. To enable this plugin just add `rtddloc` to the list of plugins in the global configuration.

![Locator plugin](/data/locator-plugin.png?raw=true "Locator plugin")

Please note that this plugin is not strictly required since `scrtdd` would relocated any manaul origins anyway (if configured to do so) and the relocated origin will appear on `scolv` as soon as ready.

Also scolv doesn't allow to create new picks when performing a relocation, so `scrtdd` plugin disable the cross-correlation on theoretical picks since those picks will not be reported on scolv.

## 6. Troubleshooting

### 6.1. Logs

Check log file: ~/.seiscomp/log/scrtdd.log 

Alternatively, when running scrtdd from the command line use the following options to see the logs on the console:

```
scrtdd [some options] --verbosity=3 --console=1
```

Verbosity 3 should be preferred to level 4, since the debug level 4 makes the logs hard to read due to the huge amount of information. Any useful information to the user is given at level 3 or above, while level 4 is only for debugging.

### 6.2. Database connection

The seiscomp database connection can be configured either inside `global.cfg`  (and thus every module knows what database to use since they inherit global.cfg), or in `scmaster.cfg`, in which case scmaster module passes the database connection string to every module when they connect to the messagin system. Since several scrtdd command line options don't need the messaging system,  scrtdd doesn't connect to it and in those cases we have to pass the database connection string to scrtdd via command line option (Since the database is still required for the inventory).
Also, specifying the database connection via command line is useful to use a local seiscomp installation to access events stored in another seiscomp installation. We can use the -d option to specify the database connection, e.g.

```
scrtdd [some options] -d  mysql://user:password@host/seiscompDbName
```

or in case of a Postgresql database:

```
scrtdd [some options] --plugins dbpostgresql -d postgresql://user:password@host/seiscompDbName
``` 


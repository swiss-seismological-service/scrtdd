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

The SCRTDD is a [SeisComP](<https://github.com/SeisComP>) extension module that implements both Real-Time Double-Difference Event relocation and classic offline Multi-Event Double-Difference relocation. The module allows multi-event and single-event mode. In multi-event mode it relocates event clusters (any set of events stored in the SeisComP database) without interfering with real-time processing and optionally stores the relocated events back in the database or simply keep the results as plain text files for external use. Also the software allows to perform Real-Time Double-Difference relocation of newly located events, one a time, as they become available through SeisComP. This is the single-event mode and in this mode each event is relocated against a background catalog (historical events for the region) that is the base on which the Double-Difference equation system is built.

The actual methods are based on the paper "Near-Real-Time Double-Difference Event Location Using Long-Term Seismic Archives, with Application to Northern California" by Felix Waldhauser and "A Double-Difference Earthquake Location Algorithm: Method and Application to the Northern Hayward Fault, California" by Waldhauser & Ellsworth.

The double-difference equation system solver uses [LSQR by Chris Paige, Michael Saunders](<https://web.stanford.edu/group/SOL/software/lsqr/>) and [LSMR by David Fong, Michael Saunders](<https://web.stanford.edu/group/SOL/software/lsmr/>) algorithms from [this](<https://github.com/tvercaut/LSQR-cpp/>) Apache Licencesed beautiful implementation by Tom Vercautereen.

## Compile

#### Seiscomp3

In order to use this module the sources have to be compiled to an executable. Merge them into the *Seiscomp3* sources and compile *Seiscomp3* as usual.
<pre>
# merge rtdd-addons and Seiscomp3
git clone https://github.com/SeisComP3/seiscomp3.git sc3-src
cd sc3-src
git submodule add -f https://github.com/swiss-seismological-service/scrtdd.git src/rtdd-addons
</pre>
For compiling Seiscomp3, please refer to https://github.com/SeisComP3/seiscomp3#compiling

#### New AGPL SeisComP

Instructions for the *new AGPL seiscomp*:
<pre>
# merge rtdd-addons and SeisComP
git clone https://github.com/SeisComP/seiscomp.git seiscomp
cd seiscomp/src/extras
git clone https://github.com/swiss-seismological-service/scrtdd.git
cd scrtdd
git checkout new_agpl_seiscomp
</pre>
For compiling SeisComP, please refer to https://github.com/SeisComP/seiscomp#build


# Getting Started

## 1. Multi-event relocation

In this section we will explain how to relocate a set of existing events in offline mode: no interaction with the running system or database will take place. The relocated events will be stored in plain text files and can be used externaly to SeisComP. Optionally it is possible to import the relocated events back into the database.

### 1.1 Dumping the candidate events to files

The multi-event relocation involves two steps: first we select the existing candidate origins and then we use scrtdd to relocated those. The selection of the origins is a user decision but what is important is that those origins exists in the seiscomp database and that we keep notes of their ID since scrtdd will us the ID to fetch all necessary information.

Once the candidate origins have been selected, we write their ids to a file like the following (a csv file in which there must be at least one column called seiscompId):

E.g. *file myCatalog.csv*

```
seiscompId
Origin/20181214107387.056851.253104
Origin/20180053105627.031726.697885
Origin/20190121103332.075405.6234534
Origin/20190223103327.031726.346363
[...]
```

Once that's done, we run the command:

```
scrtdd --dump-catalog myCatalog.csv --verbosity=3 --console=1 [db options]
```

The above command will generate three files (event.csv, phase.csv and stations.csv) which contain the information needed by scrtdd to relocate the selected events (see the `Troubleshooting` paragraph for details on how to specify a database connection):

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

Now that we have dumped the events (event.csv, phase.csv, stations.csv) we might perform some editing of those files, if required.

We now need to create a new profile in scrtdd configuration to perform the relocation. In the new profile we have to set the generated files (event.csv, phase.csv, stations.csv) as catalog.

![Catalog selection option](/data/catalog-selection2.png?raw=true "Catalog selection from raw file format")

It might be useful to know that there exists a `--dump-catalog-options` option too, which allows to use event id instead of origin id in myCatalog.csv. For example we can list all events between two dates and then ask scrtdd to extract the preferred manual origins within our profile region (our catalog):

```
# prepare myCatalog.csv
echo seiscompId > myCatalog.csv
# save all events ids between 2018-11-27 and 2018-12-14 in myCatalog.csv
scevtls --begin "2018-11-27 00:00:00" --end "2018-12-14 00:00:00" --verbosity=3 --console=1 [db options] >> myCatalog.csv
# create the catalog using only manually reviewed preferred origins within the region defined in myProfile
scrtdd --dump-catalog myCatalog.csv --dump-catalog-options 'preferred,onlyManual,any,none,myProfile' --verbosity=3 --console=1 [db options]
```
### 1.2 Avoid events dumping to flat files

event.csv, phase.csv and stations.csv files are the extended catalog format, useful when the catalog information is fully contained in those files and it doesn't require the catalog to be present on the seiscomp database. This might come in handy when using events coming from a different seiscomp installation. However It is possible to totally skip the dumping of events to files and select the catalog as a single file containing the origin ids like the following:

![Catalog selection option](/data/catalog-selection1.png?raw=true "Catalog selection from event/origin ids")


### 1.3 Relocating the candidate events

At this point we have to configure the other profile options that control the relocation process: `doubleDifferenceObservations`, which control the creation of catalog absolute travel time entries (dt.ct file in HypoDD terminology) and cross correlated differential travel times for pairs of events (dt.cc file in HypoDD terminology). The clustering options should work just fine with the default values, however some adjustments are worth it. For example we force to have only well connected events (`minNumNeigh` and `minObservationPerEvPair`), which helps in avoiding ill-defined double difference system that can be hard to solve. Then we minimize `maxNumNeigh` to reduce computation time (we can get very good results without using all the possible neighbours for every event). Finallly we like to disable the ellipsoid algorithms (`numEllipsoids=0`) since that is mostly useful in single event relocation.

![Relocation options](/data/multiEventStep2options.png?raw=true "Relocation options")

Then it is time to set the cross-correlation parameters, which require a more careful selection and it is covered in its specific paragraph.

Finally, when the configuration is ready, we can relocate the catalog with the command:

```
scrtdd --reloc-profile profileName --verbosity=3 --console=1 [db options] 
```

scrtdd will relocated the catalog and will generate another set of files reloc-event.csv reloc-phase.csv and reloc-stations.csv, which together define a new catalog with relocated origins. 

At this point we should check the relocated events (and logs) and see whether the results make sense and are satisfactory. Usually we want to keep tuning the scrtdd settings and relocate the catalog multiple times until we are sure we reached the best relocations. Having a good background catalog is the base for good real-time relocations.

As an example you can see below two catalogs before and after scrtdd relocation:


![Relocation example picture](/data/multiEventRelocationExample.png?raw=true "Relocation example")

The relocation output files (reloc-event.csv reloc-phase.csv and reloc-stations.csv) will become the background catalog used in real-time relocation and this is the only output we need to keep from the relocation process. The profile configuration can be now deleted or, in the case we want to kept it, it has to be removed from the list of active profiles (`scrtdd.activeProfiles`) to avoid interaction with real-time processing.

Now that we have a high quality background catalog we are ready to perform real-time relocation. For that we will create a new profile whose catalog will be the relocation output files:

![Catalog selection option](/data/catalog-selection3.png?raw=true "Catalog selection from raw file format")


Notes:

- Be aware that the first time you run the command it will be very slow because the waveforms have to be downloaded from the configured recordstream and saved to disk for the next runs, which will be much faster. 


### 1.4 Useful options

In addition to the options we have already seen, there are also some other interesting catalong related options.

`--dump-catalog-xml` converts a catalog to XML format, which is useful to insert the XML file in the seiscomp database or to share the catalog in a standard format.
 

E.g.
```
scrtdd --dump-catalog-xml station.csv,event.csv,phase.csv
```
or
```
scrtdd --dump-catalog-xml myCatalog.csv
```

`--merge-catalogs` and `--merge-catalogs-keepid` are useful to merge several catalogs in a single one. The difference between the two is that "keepid" properly handle repeted events. e.g.

```
scrtdd --merge-catalogs station1.csv,event1.csv,phase1.csv,station2.csv,event2.csv,phase2.csv
```

Here is a list of all the options we have seen so far:

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

## 2. Real-time single-event relocation

The first step is to create a new profile or we can modify the profile used to create the background catalog, in case we don't need that anymore. If we decide to go for a new profile we might want to copy the waveform cache folder of the background catalog profile to the new real-time profile cache folder to avoid re-downloading all the catalog waveforms (see the `Performance` paragraph for more details).

The waveform cache folder is located in `workingDirectory/profileName/wfcache/` (e.g. `~/seiscomp3/var/lib/rtdd/myProfile/wfcache/`).

In real-time processing scrtdd relocates new origins against a background catalog of high quality events. If we don't already have those high quality events we can generate the background catalog via multi-event relocation coverred in the previous section.
 

### 2.1. Selecting a background catalog from existing origins

If we already have a high quality catalog in the database, we can specify it in scrtdd configuration as a path to a file containing a list of origin ids. The file format must be a csv file in which there must be at least one column called seiscompId, containing the origins ids we like to use as background catalog. 

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

Having the background catalog in the database is easy: once the multi-event relocation is completed, we can transform the relocated catalog to XML format (see --dump-catalog-xml option) and insert XML in the seiscomp database. While it is neat to have the background catalog in the seiscomp database, this approach is inconvenient in the common scenario where the background catalog is re-generated from time to time and it is not fixed.

### 2.2. Generating a background catalog

In the more general case in which we don't already have a background catalog, we have to build one using scrtdd. To do so simply follow the instuctions in the multi-event relocation paragraph and then configure the profile to use as background catalog the relocated multi-event results (reloc-event.csv, phase.csv, station.csv):

![Catalog selection option](/data/catalog-selection3.png?raw=true "Catalog selection from raw file format")
 
### 2.3 Testing

Real time relocation uses the same configuration we have seen in full catalog relocation, but real time relocation is done in two steps, each one controlled by a specific configuration.

Step 1: location refinement. In this step scrtdd performs a preliminary relocation of the origin using only catalog absolute travel time entries (dt.ct only).

Step 2: the refined location computerd in the previous step is used as starting location to perform a more precise relocation using both catalog absolute travel times (dt.ct) and differential travel times from cross-correlation (dt.cc). If step1 fails, step2 is attempted anyway.

If step2 completes successfully the relocated origin is sent to the messaging system.

You can see below the most important parameters we want to configure. We enable the ellipsoid clustering algorithms selecting `numEllipsoids != 0` and we define the number of neighbour events to use for relocation. We also increased `maxEllipsoidSize` in `doubleDifferenceObservationsNoXcorr` to accomodate the possible large error in the automatic origin locations.
 

![Relocation options](/data/step1options.png?raw=true "Relocation options")
![Relocation options](/data/step2options.png?raw=true "Relocation options")

You might consider testing the configuration relocating some existing events to make sure the parameters are suitable for your use case. To test the real time relocation there are two command line options which relocate existing origins:

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

#### 2.3.1 Relocate origin ID and send the relocation to messaging system for further processing

If we want to process an origin we can run the following command and then check on scolv the relocated origin (the messaging system must be active). This is mostly useful when we want to relocate an origin on a running system and keep the relocation:

```
scrtdd --origin-id someOriginId --verbosity=3 --console=1 [db options] 
```

#### 2.3.2 Relocate origin ID but do not send the relocation (debug)

As above but add `--test`

```
scrtdd --origin-id someOriginId --test --verbosity=3 --console=1 [db options]
``` 

#### 2.3.3 Relocate origin ID and store the result to XML file

For testing purpose we are more likely interested in not interfering with the database and the messaging system so we can use the `--ep` option and the relocated origin will be saved as a xml file. We can finally open the xml file with scolv for inspection:

```
scrtdd --origin-id someOriginId --ep - --verbosity=3 --console=1 [db options] >  relocated-origin.xml
```

#### 2.3.4 Relocate XML file and store the result to XML file

Alternatively the `--ep` option (without `-O`) can process all origins contained in the input XML file:

```
scrtdd --ep origin.xml --verbosity=3 --console=1 [db options] > relocated-origin.xml
```

And we can use the scxmldump to dump and existing origin id to file:

```
# dump origin
scxmldump -fPAMF -p -O originId -o origin.xml --verbosity=3  --console=1
```

#### 2.3.5 Notes

Also, as explained in the cross-correlation settings paragraph, we can add the ```--debug-wf``` option to the above commands to dump and instepct the waveforms used during cross-correlation.

As an example you can see below the single event relocation of several manually reviewed origins (when relocating automatic origins the quality and number of relocated origins is certainly lower).

![Single event relocation example picture](/data/singleEventRelocationExample.png?raw=true "Single Event Relocation example")

Once we are happy with the configuration we can simply enable and start scrtdd as any other module and it will start relocating origins as soon as they arrive in the messsaging system

### 2.4 Phase update

scrtdd uses cross-correlation to detect phases at stations with no associated picks and to fix the pick time and uncertainty of automatic picks. Those features are especially useful in real-time to increase the quality and number of double-difference observations when automatic origins have only few picks/phases.

For automatic picks, the pick time is updated accordingly to the average lag detected by all the good (above configured threshold) cross-correlation results. Since the real-time events are cross-correlated against catalog events, which have good manual picks, the updated pick time is expected to improve. The pick uncertainty is derived from the uncertainties of catalog-events. If no cross-correlations above the configured threshold are found, the pick is kept untouched.

For stations with no associated phases, scrtdd computes theoretical picks. Those are then cross-correlated against the catalog event ones. Every theoretical pick that has at least one good cross-correlation result is added to the relocated origin, with pick time and uncertainties derived from catalog phases (similarly to what is done for automatic picks). Those "good" theoretical picks and thus used in the double-difference system inversion. Theoretical picks that have no good cross-correlation results are simply discarded.

Picks that have been updated or created by scrtdd are identifiable by a `x` suffix (Px, Sx).

Manual picks are never modified.

### 2.5 RecordStream configuration

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


## 3. Cross-correlation, waveform filtering and signal to noise ratio options

Good cross-correlation results are needed to achieve high resolution double-difference observations, which in turn results in high quality relocations.

The configuration of cross-correlation parameters require some trial and error, but once they are properly tuned for the multi-event scenario they can be kept identical in the single-event relocation, except for the parameter ```maxDelay``` that should be increased in real-time relocation to accomodate the larger uncertainty of automatic picks (unless scrtdd is used to only relocate manually reviewed origins). 

Here is an example configuration:

![Relocation options](/data/xcorr.png?raw=true "Relocation options")


### 3.1 Logs

Some cross-correlation statistics are printed in both multi-event and single-event mode. Those can be seen in the log file or in the console output if ```--console=1 --verbosity=3```) is used.

```
[info] Cross-correlation statistics: performed 40361, waveforms with Signal to Noise ratio too low 2435, waveforms not available 98
[info] Total xcorr 40361 (P 59%, S 41%) success 28% (11499/40361). Successful P 22% (5300/23844). Successful S 38% (6199/16517)
[info] xcorr on actual picks 24784/40361 (P 60%, S 40%) success 37% (9186/24784). Successful P 31% (4629/14761). Successful S 45% (4557/10023)
[info] xcorr on theoretical picks 15577/40361 (P 58%, S 42%) success 15% (2313/15577). Successful P 7% (671/9083). Successful S 25% (1642/6494)
```

The statistics are broken down in actual picks and theoretical picks. This is because scrtdd computes theoretical picks that are cross-correlated together with detected picks. This is useful to increase the number of double-difference observations.  See "Phase update" paragraph for details.


### 3.2 Eval-xcorr command

A more sophisticated method for evaluating the settings is the `--eval-xcorr` command (here we use `--verbosity=2` because the statistics are printed at this log level, useful to avoid being overwhelmed by too much information):

```
scrtdd --eval-xcorr profileName --verbosity=2 --console=1
```

Example output:

```
[...]
[warning] <<<Final stats>>>
Cumulative stats: #pha  29622 pha good CC  57% avg coeff 0.83 avg goodCC/ph  4.8
Cumulative stats P ph: #pha  16382 pha good CC  54% avg coeff 0.81 avg goodCC/ph  5.1
Cumulative stats S ph: #pha  13240 pha good CC  62% avg coeff 0.86 avg goodCC/ph  4.4

Cross-correlated Phases by inter-event distance in 0.10 km step
 EvDist [km]  #Phases GoodCC AvgCoeff(MAD) GoodCC/Ph(MAD) 
 0.00-0.10      14159    65%  0.88 (0.06)    2.3 ( 1.3)
 0.10-0.20      16015    55%  0.85 (0.06)    2.0 ( 1.0)
 0.20-0.30      11161    41%  0.84 (0.06)    1.5 ( 0.7)
 0.30-0.40       6744    33%  0.83 (0.06)    1.4 ( 0.6)
 0.40-0.50      10579    27%  0.82 (0.05)    1.4 ( 0.6)
 0.50-0.60      13046    30%  0.81 (0.05)    1.6 ( 0.7)
 0.60-0.70       7672    23%  0.81 (0.05)    1.3 ( 0.5)
 0.70-0.80       4994    19%  0.81 (0.05)    1.2 ( 0.3)
 0.80-0.90       4204    18%  0.81 (0.05)    1.1 ( 0.2)
 0.90-1.00      11315    21%  0.81 (0.05)    1.3 ( 0.5)
 1.00-1.10       9150    17%  0.80 (0.05)    1.3 ( 0.5)
 1.10-1.20       6253    16%  0.81 (0.05)    1.2 ( 0.4)
 1.20-1.30       4900    14%  0.80 (0.05)    1.2 ( 0.3)
 1.30-1.40       3732    13%  0.80 (0.05)    1.1 ( 0.2)
 1.40-1.50       2848    14%  0.80 (0.05)    1.1 ( 0.2)
 1.50-1.60       1961    12%  0.79 (0.05)    1.1 ( 0.1)
 1.60-1.70       1403    13%  0.79 (0.05)    1.1 ( 0.2)
 1.70-1.80       1602    12%  0.80 (0.06)    1.0 ( 0.1)
 1.80-1.90       3218    14%  0.80 (0.05)    1.1 ( 0.2)
 1.90-2.00       7937    12%  0.79 (0.05)    1.3 ( 0.4)
 2.00-2.10       6294    11%  0.80 (0.05)    1.2 ( 0.3)
 2.10-2.20       5010    11%  0.80 (0.05)    1.2 ( 0.3)
 2.20-2.30       3465     8%  0.78 (0.04)    1.0 ( 0.1)
 2.30-2.40       3054     9%  0.79 (0.05)    1.1 ( 0.2)
 2.40-2.50       2811     8%  0.79 (0.05)    1.1 ( 0.1)
 2.50-2.60       2393     7%  0.78 (0.05)    1.1 ( 0.1)
 2.60-2.70       2778     7%  0.79 (0.05)    1.1 ( 0.1)
[...]

Cross-correlated Phases by event to station distance in 3.00 km step
StaDist [km]  #Phases GoodCC AvgCoeff(MAD) GoodCC/Ph(MAD) 
  0-3               9    44%  0.74 (0.03)    2.0 ( 1.0)
  3-6             395    55%  0.81 (0.05)    3.3 ( 2.0)
  6-9            2827    75%  0.83 (0.05)    7.5 ( 5.5)
  9-12           3477    69%  0.84 (0.05)    6.0 ( 4.6)
 12-15           3216    65%  0.83 (0.05)    4.1 ( 2.7)
 15-18           1619    56%  0.84 (0.05)    4.1 ( 2.9)
 18-21           1367    52%  0.83 (0.05)    4.5 ( 3.1)
 21-24           1888    61%  0.83 (0.05)    4.5 ( 2.9)
 24-27           1857    55%  0.84 (0.05)    4.2 ( 2.8)
 27-30           2082    62%  0.84 (0.05)    5.1 ( 3.4)
 30-33           1213    53%  0.84 (0.05)    5.0 ( 3.8)
 33-36           1935    59%  0.84 (0.05)    5.4 ( 3.8)
 36-39           1361    49%  0.83 (0.05)    3.5 ( 2.2)
 39-42            871    44%  0.83 (0.05)    3.6 ( 2.4)
 42-45            774    46%  0.84 (0.06)    3.7 ( 2.8)
 45-48            933    54%  0.84 (0.06)    2.8 ( 1.7)
 48-51            452    43%  0.85 (0.06)    3.1 ( 2.0)
 51-54            431    50%  0.85 (0.06)    2.6 ( 1.5)
 54-57            376    48%  0.85 (0.05)    3.1 ( 2.1)
 57-60            420    41%  0.86 (0.06)    2.3 ( 1.3)
 60-63            311    45%  0.86 (0.06)    2.7 ( 1.9)
 63-66            233    42%  0.86 (0.06)    1.9 ( 1.1)
 66-69            182    41%  0.85 (0.06)    1.9 ( 1.1)
 69-72            141    53%  0.86 (0.07)    1.9 ( 1.2)
 72-75            122    35%  0.85 (0.05)    1.6 ( 0.9)
 75-78            123    46%  0.87 (0.06)    1.5 ( 0.7)
 78-81            124    39%  0.86 (0.05)    2.2 ( 1.6)
 81-84             97    46%  0.84 (0.06)    1.3 ( 0.5)
 84-87             89    25%  0.88 (0.06)    1.5 ( 0.7)
 87-90             85     0%  0.00 (0.00)    0.0 ( 0.0)
[...]

Cross-correlations by station
Station       #Phases GoodCC AvgCoeff(MAD) GoodCC/Ph(MAD) 
4D.MH36.A          145    52%  0.85 (0.07)    2.4 ( 1.2)
4D.MH44.A          179    56%  0.86 (0.05)    2.5 ( 1.5)
4D.MH48.A          126    56%  0.86 (0.07)    2.1 ( 1.1)
4D.MH54.A           56    57%  0.86 (0.07)    1.5 ( 0.6)
4D.RA41.            42    60%  0.90 (0.04)    3.8 ( 1.8)
4D.RA42.             3    67%  0.85 (0.00)    1.0 ( 0.0)
4D.RA43.            55    64%  0.90 (0.06)    3.0 ( 2.1)
8C.CI18.00         102    25%  0.85 (0.05)    1.8 ( 0.8)
8C.CI19.00         128    48%  0.87 (0.05)    2.3 ( 1.2)
8C.CI20.00          69    54%  0.90 (0.06)    3.0 ( 1.2)
8D.AMID2.          117    43%  0.86 (0.07)    1.4 ( 0.5)
8D.AMIDI.           79    78%  0.84 (0.05)    6.0 ( 4.8)
8D.BT050.           39    51%  0.85 (0.05)    2.5 ( 1.4)
8D.BT175.           41    46%  0.84 (0.03)    2.2 ( 0.9)
8D.BT225.           36    50%  0.89 (0.04)    3.8 ( 1.9)
8D.BTAO.            64    55%  0.85 (0.06)    2.8 ( 1.3)
8D.BTNF.           106    49%  0.86 (0.06)    4.9 ( 3.5)
8D.CHBM.            27    59%  0.85 (0.02)    3.8 ( 1.1)
8D.ELM0.            46    83%  0.84 (0.05)    3.6 ( 1.6)
8D.ELM1.            36    86%  0.84 (0.06)    4.9 ( 1.9)
8D.GSF03.          132    52%  0.88 (0.06)    8.0 ( 5.6)
8D.LULY1.            4    50%  0.97 (0.00)    1.0 ( 0.0)
8D.LULY2.            2   100%  0.87 (0.00)    1.0 ( 0.0)
8D.LULY3.            1     0%  0.00 (0.00)    0.0 ( 0.0)
8D.MFERR.          534    68%  0.84 (0.05)    5.0 ( 3.5)
8D.NDBB1.            7     0%  0.00 (0.00)    0.0 ( 0.0)
8D.NDBU.             3    33%  0.87 (0.00)    8.0 ( 0.0)
8D.NDMT1.            5    80%  0.82 (0.11)    1.0 ( 0.0)
8D.NDTR1.           37    49%  0.87 (0.07)    6.8 ( 7.0)
8D.NNFS.            30    57%  0.87 (0.08)    2.4 ( 1.6)
8D.NNFT.            78    22%  0.88 (0.07)    4.4 ( 4.7)
[...]
```

It is especially interesting to compare the results before/after the catalog has been relocated. The new statistics should show better performance for events close to each others and gradually worsen with increasing inter-event distance. That is an indirect measure of the quality of the relocation as explained in Waldhauser & Ellsworth's paper.

### 3.3 Waveforms inspection

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
[info] Computing cross-correlation differential travel times for event 12211
[info] xcorr: event 12211 sta    S ESION dist    5.69 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   XY VET01 dist    5.77 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   XY VET02 dist    6.00 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   CH  SIOV dist    6.12 [km] - 2 P phases, mean coeff 0.73 lag 0.16 (events: 4862 3156 )
[info] xcorr: event 12211 sta   CH  SIOM dist    6.27 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   CH  SARC dist    7.25 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   CH  SIOO dist    7.32 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   8D  RAW2 dist    7.98 [km] - 7 P phases, mean coeff 0.78 lag 0.00 (events: 3337 3264 3310 4862 3156 4861 4036 )
[info] xcorr: event 12211 sta   XY LEO01 dist    9.87 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   CH SAYF2 dist    9.95 [km] - 2 S phases, mean coeff 0.77 lag 0.07 (events: 4051 4036 )
[info] xcorr: event 12211 sta   CH  SHER dist   11.43 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   8D  RAW4 dist   12.07 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   CH SENIN dist   12.70 [km] - 4 P phases, mean coeff 0.74 lag 0.00 (events: 3264 3310 4862 4861 )
[info] xcorr: event 12211 sta   XY SAI01 dist   13.67 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   CH STSW2 dist   14.58 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   XY SIE01 dist   15.12 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   CH GRYON dist   15.49 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   8D  RAW3 dist   16.17 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   CH FULLY dist   17.54 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   CH  SIEB dist   18.43 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   8D  RAW1 dist   21.19 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   CH   DIX dist   21.75 [km] - 4 P phases, mean coeff 0.84 lag -0.00 (events: 3264 4862 3447 4861 )
[info] xcorr: event 12211 sta   CH   DIX dist   21.75 [km] - 1 S phases, mean coeff 0.78 lag -0.04 (events: 3168 )
[info] xcorr: event 12211 sta   XP ILL07 dist   22.72 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   CH LAVEY dist   22.85 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   CH VANNI dist   23.40 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   XP ILL15 dist   24.20 [km] - 1 P phases, mean coeff 0.85 lag -0.05 (events: 3264 )
[info] xcorr: event 12211 sta   XP ILL05 dist   24.20 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   XP ILL18 dist   24.27 [km] - 1 P phases, mean coeff 0.76 lag 0.04 (events: 3264 )
[info] xcorr: event 12211 sta   XP ILL08 dist   24.27 [km] - low corr coeff pairs
[info] xcorr: event 12211 sta   CH  SMAO dist   25.22 [km] - low corr coeff pairs
[...]
[info] Event 12211 total phases 27: created 6 (3 P and 3 S) from theoretical picks

[info] Cross correlation performed 275, phases with Signal to Noise ratio too low 66, phases not available 54 (waveforms downloaded 183, waveforms loaded from disk cache 1020)
[info] Total xcorr 275 (P 67%, S 33%) success 11% (31/275). Successful P 14% (26/185). Successful S 6% (5/90)
[info] xcorr on actual picks 115/275 (P 100%, S 0%) success 18% (21/115). Successful P 18% (21/115). Successful S -nan% (0/0)
[info] xcorr on theoretical picks 160/275 (P 44%, S 56%) success 6% (10/160). Successful P 7% (5/70). Successful S 6% (5/90)
```

For comparison we can always find the raw waveforms (not processed) fetched from the configured recordstream and used as a cache in `workingDirectory/profileName/wfcache/` (e.g. `~/seiscomp3/var/lib/rtdd/myProfile/wfcache/`):
* `NET.ST.LOC.CH.startime-endtime.mseed`


## 4. Performance and waveforms data

scrtdd spends most of the relocation time downloading waveforms (unless the recordstream points to a local disk storage) and the rest is shared betweeen cross-correlation and double difference system inversion. For this reason the waveforms are cached to disk after being downloaded so that there is no need to download them again. This is obviously true for catalog phases waveforms that are re-used over and over, but that's not true for real-time events waveforms and also for some temporary waveforms (e.g. theoretical phases) that are never saved. The cache folder is `workingDirectory/profileName/wfcache/` (e.g. /installation/path/seiscomp3/var/lib/rtdd/myProfile/wfcache/). 

However, during the parameters tuning phase, the user performs relocations (both single-event and multi-event) from the command line several times to find the best configuration. For those speacial cases even the temporary waveforms are saved to make the process much faster. A different folder is used: `workingDirectory/profileName/tmpcache/` which can be deleted after the parameter tuning phase. The commands that store temporary waveforms are: `--reloc-profile`, `--ep`, `-O`, `--origin-id`.

Due to the time required to download waveforms, when relocating events in real-time, two configuration options have a huge impact on the performance: `doubleDifferenceObservations.clustering` and `crosscorrelation.s-phase.components`. `doubleDifferenceObservations.clustering` is relevant because we can specify how many neighbouring events and how many phases we want to use, which consequently determines the number of cross-correlations to perform, the waveforms to download and the size of the input for the double difference inversion. `crosscorrelation.s-phase.components` defines the components we want to use in the cross-correlation of the `S` phases. Usign the `T` components means scrtdd has to download the additional two components to perform the projection to ZRT, which might be more accurate than 'Z' only.


### 4.1 Catalog waveforms preloading

When scrtdd starts for the first time it loads all the catalog waveforms and store them to disk (once the waveforms are downloaded they will not be downloaded again unless the files are deleted from disks). In this way the waveforms become quickly available in real-time without the need to access the recordstream (which takes a long time!). However if the option `performance.profileTimeAlive` is greater than 0, the catalog waveforms will be loaded on demand and not when scrtdd starts. In this case we might decide to pre-download all waveforms anyway, before starting the module, using the following option:

```
scrtdd --help
  --load-profile-wf arg                 Load catalog waveforms from the 
                                        configured recordstream and save them 
                                        into the profile working directory
``` 

e.g.

```
scrtdd --load-profile-wf myprofile
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


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

rtDD is a [SeisComP](<https://github.com/SeisComP>) extension module that implements both Real-Time Double-Difference Event relocation and classic offline Multi-Event Double-Difference relocation.

The actual methods are based on the paper "Near-Real-Time Double-Difference Event Location Using Long-Term Seismic Archives, with Application to Northern California" by Felix Waldhauser and "A Double-Difference Earthquake Location Algorithm: Method and Application to the Northern Hayward Fault, California" by Waldhauser & Ellsworth.

The double-difference equation system solver uses [LSQR by Chris Paige, Michael Saunders](<https://web.stanford.edu/group/SOL/software/lsqr/>) and [LSMR by David Fong, Michael Saunders](<https://web.stanford.edu/group/SOL/software/lsmr/>) algorithms from [this](<https://github.com/tvercaut/LSQR-cpp/>) Apache Licensed beautiful implementation by Tom Vercautereen.

rtDD also supports [NonLinLoc by Anthony Lomax](<http://alomax.free.fr/nlloc/>) grid file format alongside the formats natively supported by SeisComP (LOCSAT and libtau). A feature that enables 3D velocity models support within the tool.

## Compile

In order to use this module the sources have to be compiled to an executable. Merge the scrtdd code into the *SeisComP* sources and compile *SeisComP* as usual.

#### AGPL SeisComP

Instructions for the *new AGPL seiscomp*:
<pre>
# get SeisComP
git clone https://github.com/SeisComP/seiscomp.git seiscomp
# merge rtdd to SeisComP
cd seiscomp/src/extras
git clone https://github.com/swiss-seismological-service/scrtdd.git
cd scrtdd
# now checkout the tag pointing to the latest version: vX.Y.Z (e.g. v1.3.0)
git checkout vX.Y.Z
</pre>
For compiling SeisComP, please refer to https://github.com/SeisComP/seiscomp#build

#### Old versions of Seiscomp3

<pre>
# get Seiscomp3
git clone https://github.com/SeisComP3/seiscomp3.git sc3-src
# merge rtdd to Seiscomp3
cd sc3-src
git submodule add -f https://github.com/swiss-seismological-service/scrtdd.git src/rtdd-addons
cd src/rtdd-addons
# now checkout the tag pointing to the latest version: vX.Y.Z_sc3 (_sc3 for Seiscomp3 compatible versions)
git checkout vX.Y.Z_sc3
</pre>
For compiling Seiscomp3, please refer to https://github.com/SeisComP3/seiscomp3#compiling

# Getting Started

* [1.Multi-event relocation](#1-multi-event-relocationi)
* [2.Real-time single-event relocation](#2-real-time-single-event-relocation)
* [3.Cross-correlation](#3-cross-correlation)
* [4.Waveform data and recordStream configuration](#4-waveform-data-and-recordstream-configuration)
* [5.Database connection](#5-database-connection)
* [6.Troubleshooting](#6-troubleshooting) 
* [7.Custom velocity models](#7-custom-velocity-models)
* [8.Scolv Locator plugin](#8-scolv-locator-plugin) 

## 1. Multi-event relocation

Long story short:

* Use `sclistorg` command to select the events to be relocated. E.g. `sclistorg [options] > myCatalog.csv`

* Create a `scrtdd` profile (e.g. use `scconfig` GUI). The profile will contain the settings for the relocation: the default values provided by `scrtdd` are meant to be a good starting choice, so there is no need to tweak every parameter. However it is a good choice to configure a custom velocity model (`solver.travelTimeTable`)

* Use `scrtdd` to relocate the events. E.g. `scrtdd --reloc-catalog  myCatalog.csv --profile myProfile [db and waveforms options]`

* Use the relocation output (`reloc-event.csv`, `reloc-phase.csv` and `reloc-stations.csv`) as you please

To relocate external (non-SeisComP) data refer to [this paragraph](#134-relocating-external-data).

### The long story

The multi-event relocation consists of two steps: selecting the candidate origins and using `scrtdd` to relocate those. For the former task an utility `sclistorg` might come in handy. For the latter we need to configure a `scrtdd` profile where the relocation options are stored and then run the command line option `--reloc-catalog`. That's it.

The output will be another catalog containing the relocated origins. The configuration is an easy task since the default options should already provide sensible solutions. The input origins might come from different sources: a SeisComP database (local or remote), a xml file, or a `scrtdd` specific format (station.csv,event.csv,phase.csv triplet, explained later).

In multi-event mode there is no interaction neither with the running SeisComP modules nor with database (except for reading the inventory). It is a safe operation and allow for easy experimenting. The relocated events will be stored in plain text files, so that they can be analysed  externally to SeisComP, but they can also be stored back into the database. That is optional.

`scrtdd` supports also the standard options (`--inventory-db inventory.xml --config-db config.xml`) that allows to read the station inventory from plain files, making the database fully optional.


### 1.1 How to get the origin ids?

There is a tool that is installed alongside `scrtdd`, called `sclistorg`, that is useful for listing origin ids satisfying certain criteria, such as time period, geographic area, author, agency and so on. E.g.

```sh
# list the preferred origin ids for all events between 2018-11-27 and 2018-12-14
sclistorg --begin "2018-11-27 00:00:00" --end "2018-12-14 00:00:00" --org-type preferred [db options]

# select also the event type and the accepted agencies
sclistorg --begin "2018-11-27 00:00:00" --end "2018-12-14 00:00:00" --org-type preferred \
          --ev-type "earthquake,quarry blast" --inc-agency Agency1,Agency2 [db options]

# select an area of interest, a rectangle minLat,minLon,maxLat,maxLon
sclistorg --begin "2018-11-27 00:00:00" --end "2018-12-14 00:00:00" --org-type preferred \
          --area 46.0,8.5,46.5,8.7 [db options]

```
See `sclistorg --help` for a full list of options.

```
Events:
  --begin arg                   specify the lower bound of the time interval
  --end arg                     specify the upper bound of the time interval
  --modified-after arg          select events modified after the specified time
  --ev-type arg                 include only events whose type is one of the
                                values provided (comma separated list)
  --simple                      Print only origin ids

Origins:
  --org-type arg                preferred, last or first (default is preferred)
  --manual-only                 Include only manual origins
  --auto-only                   Inlude only automatic origins
  --inc-author arg              include only origins whose author is one of the
                                values provided (comma separated list)
  --excl-author arg             exclude origins whose author is one of the
                                values provided (comma separated list)
  --inc-method arg              include only origins whose methodID is one of
                                the values provided (comma separated list)
  --excl-method arg             exclue origins whose methodID is one of the
                                values provided (comma separated list)
  --inc-agency arg              include only origins whose agencyID is one of
                                the values provided (comma separated list)
  --excl-agency arg             exclude origins whose agencyID is one of the
                                values provided (comma separated list)
  --area arg                    Include only origins in the rectangular area
                                provided: MinLat,MinLon,MaxLat,MaxLon
```

### 1.2 Preparing the candidate events

The origin ids of the events to be relocated must be stored in a format that `scrtdd` understands.

One of the compatible formats is a text file containing the origin IDs (`sclistorg` output is compatible with that). `scrtdd` will use the origin IDs in the file to fetch all necessary information from the SeisComP database.

E.g. *file myCatalog.csv* (a mandatory column named `seiscompId` is required, but other column might be present too).

```
seiscompId
Origin/20181214107387.056851.253104
Origin/20180053105627.031726.697885
Origin/20190121103332.075405.6234534
Origin/20190223103327.031726.346363
[...]
```

There is another format we can use to store a catalog. This format contains the full origins information, not only the origin ids. So, once the files are generated, there is no need to access the database anymore; so this format is quite fast to load. We can instruct `scrtdd` to generate such a format with the following command:

```
scrtdd --dump-catalog myCatalog.csv --verbosity=3 --console=1 [db options]
```

The above command will generate three files (*event.csv*, *phase.csv* and *stations.csv*) which contain all the information needed by `scrtdd`. 

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
With this format it is possible to relocate events that are not stored in any SeisComP database, since all the origins information are contained in those files.

Finally, the events to be relocated can also be stored in SeisComP XML format. Please refer to the official SeisComP  documentation of `scxmldump`, a very convenient tool for dumping events to XML file.


### 1.3 Relocating the candidate events

Before performing the relocation we need to create a new profile in the `scrtdd` configuration where it is possible to select the values for the relocation steps: dubble-differene system creation, cross-correlation and solver.

![Profile options](/data/img/configOverview.png?raw=true "Profile options")

The default values provided by `scrtdd` are meant to be a good starting choice, so there is no need to tweak every parameter. However it is a good choice to configure a custom velocity model (`solver.travelTimeTable`). The cross-correlation parameters are described in a dedicated paragraph. Finally, when the configuration is ready, we can relocate the catalog with the following commands...

#### 1.3.1 Relocating a file containing a list of origin ids

```
scrtdd --reloc-catalog myCatalog.csv --profile myProfile \
       --verbosity=3 --console=1 [db options] 
```

E.g. *file myCatalog.csv*

```
seiscompId
Origin/20181214107387.056851.253104
Origin/20180053105627.031726.697885
[...]
```
 
#### 1.3.2 Relocating the station.csv,event.csv,phase.csv triplet

```
# station.csv,event.csv,phase.csv are generated with `scrtdd --dump-catalog`
scrtdd --reloc-catalog station.csv,event.csv,phase.csv --profile myProfile \
       --verbosity=3 --console=1 [db options] 
```

#### 1.3.3 Relocating an XML/SCML file

Events are stored in a XML files in [SCML format](https://www.seiscomp.de/doc/base/glossary.html#term-SCML). It is possible to convert between different formats with [sccnv command](https://www.seiscomp.de/doc/apps/sccnv.html).

```
# events.xml contais the events data (scxmldump command)
# myCatalog.csv contains the origin ids inside events.xml we want relocate
scrtdd --reloc-catalog myCatalog.csv --ep events.xml --profile myProfile \
       --verbosity=3 --console=1 [db options] 
```

#### 1.3.4 Relocating external data

To relocate external (non SeisComP) data three pieces of information need to be provided: events data, waveform data and inventory information:

* events data has to be provided in [SCML format](https://www.seiscomp.de/doc/base/glossary.html#term-SCML). It is possible to convert between different formats with [sccnv command](https://www.seiscomp.de/doc/apps/sccnv.html). Events data can be passed to `scrtdd` via `--ep events.xml` option together with `--reloc-catalog` option
* alternatively the events data can be converted to a station.csv,event.csv,phase.csv file triplet, explained in the previous paragraphs and passed to `scrtdd` via `--reloc-catalog station.csv,event.csv,phase.csv` option
* Waveform data can to be provided via `-I RecordStream` command line option and the RecordStream cab be any of the [SeisComP supported ones](https://www.seiscomp.de/doc/apps/global_recordstream.html#global-recordstream)
* [Inventory information](https://www.seiscomp.de/doc/base/concepts/inventory.html) has be converted from an external format into SeisComP own station meta-data XML format called inventory ML. This can be passed to `scrtdd` via `--inventory-db inventory.xml` (or stored in the SeisComP database)


Examples:

```
scrtdd --reloc-catalog station.csv,event.csv,phase.csv --profile myProfile \
       -I sdsarchive:///home/sysop/seiscomp/var/lib/archive \
       --inventory-db inventory.xml \
       --verbosity=3 --console=1
```

```
# myCatalog.csv contains the origin ids inside events.xml we want relocate
scrtdd --reloc-catalog myCatalog.csv --ep events.xml --profile myProfile \
       -I fdsnws://service.iris.edu:80/fdsnws/dataselect/1/query  \
       --inventory-db inventory.xml \
       --verbosity=3 --console=1
```

### 1.4 Review of results

Independently on how the input events are provided, `scrtdd` will output a set of files *reloc-event.csv*, *reloc-phase.csv* and *reloc-stations.csv*, these contain the relocated catalog.

At this point we should check the relocated events which contains several statistics and the logs to see whether the results make sense and are satisfactory. Usually we want to keep tuning the `scrtdd` settings and relocate the catalog multiple times until we are sure we reached the best relocations. Having a good background catalog is the base for good real-time relocations.

**Notes**:
Be aware that the first time you try to relocate a catalog it can be very slow because the waveforms have to be downloaded from the configured recordStream and saved to disk for the next runs, which will be much faster. 

To enable `scrtdd` in real-time the output files reloc-event.csv reloc-phase.csv and reloc-stations.csv are required as they will become the background catalog for the real-time/single-event relocation.

![Catalog selection option](/data/img/catalog-selection3.png?raw=true "Catalog selection from raw file format")


### 1.5 Useful options

In addition to the options we have already seen, there are also some other useful ones.

`--xmlout` option can be used in combination with `--reloc-catalog` to generate a XML output, which is useful to later insert the relocated catalog in a SeisComP database (e.g. scdb command).

E.g.
```
scrtdd --reloc-catalog myCatalog.csv --profile myProfile \
       --verbosity=3 --console=1 [db options]
       --xmlout > relocated-catalog.xml
```

`--merge-catalogs` and `--merge-catalogs-keepid` are useful to merge several catalogs into a single one. 

```
scrtdd --merge-catalogs station1.csv,event1.csv,phase1.csv,station2.csv,event2.csv,phase2.csv
```

Here is a list of all the options we have seen so far:

```
scrtdd --help

Mode:

  --reloc-catalog arg                   Relocate the catalog passed as argument
                                        in multi-event mode. The input can be a
                                        single file (containing seiscomp origin
                                        ids) or a file triplet
                                        (station.csv,event.csv,phase.csv). For
                                        events stored in a XML files add the
                                        --ep option. Use in combination with
                                        --profile

  --ep arg                              Event parameters XML file for offline
                                        processing of contained origins
                                        (implies --test option). Each contained
                                        origin will be processed in
                                        signle-event mode unless
                                        --reloc-catalog is provided, which
                                        enable multi-event mode.

ModeOptions:

  --profile arg                         To be used in combination with other 
                                        options: select the profile 
                                        configuration to use

  --xmlout                              Enable XML output when combined with 
                                        --reloc-catalog or --oring-id options

Catalog:

  --dump-catalog arg                    Dump the seiscomp event/origin id file 
                                        passed as argument into a catalog file 
                                        triplet (station.csv,event.csv,phase.cs
                                        v).

  --merge-catalogs arg                  Merge in a single catalog all the 
                                        catalog file triplets 
                                        (station1.csv,event1.csv,phase1.csv,sta
                                        tion2.csv,event2.csv,phase2.csv,...) 
                                        passed as arguments.

  --merge-catalogs-keepid arg           Similar to the --merge-catalogs option 
                                        but events keep their ids. If multiple 
                                        events share the same id, subsequent 
                                        events will be discarded.
```

### 1.6 Periodic update of the multi-event relocation

For real-time monitoring it is useful to periodically update the multi-event relocation, so that the news events are continuously included in the double-difference inversion. This is not only useful for keeping track of what is happening in a region, but it is crucial for real-time relocation, where new origins are relocated using a background catalog: the more update the background catalog is, the better the results.

For this purpose it might come in handy [this script](/data/scripts/generate-catalog.sh), that can be easily adapted to the specific use case.

### Examples

Here are two catalogs before and after `scrtdd` relocation:

![Relocation example picture](/data/img/multiEventRelocationExample.png?raw=true "Relocation example") 


The unit testing folder contains the code to generate some tests with synthetic data:

![Synthetic Data Relocation example picture](/data/img/multiEventRelocationSyntDataExample.png?raw=true "Synthetic Data Relocation example") 


## 2. Real-time single-event relocation

Long story short:

* Use the multi-event relocation feature to prepare a background catalog

* Create a `scrtdd` profile or use the same profile used for generating the background catalog, then set the profile background catalog and add the profile to the list of active real-time profiles (`activeProfiles` parameter). The default profile parameter values are meant to be a good starting choice, so there is no need to tweak them heavily. However it is a good choice to configure a custom velocity model (`solver.travelTimeTable`)

* Make sure to read "Avoiding Relocation Loops" paragraph to avoid a potential issue

* Enable and start `scrtdd` (`seiscomp enable scrtdd`, `seiscomp start scrtdd`)

### The long story

To enable the real-time processing a profile should be created and enabled by including it in `scrtdd.activeProfiles` option.
 
In real-time processing `scrtdd` relocates new origins, one a time as they occur, against a background catalog of high quality events. Those high quality events can be generate via multi-event relocation, which has already been covered in the previous sections.

Real time relocation uses the same configuration we have seen in full catalog relocation, but real time relocation is done in two steps:

**Step 1**: location refinement. In this step `scrtdd` performs a preliminary relocation of the origin where the differential travel times in the double-difference system are derived from the pick times.

**Step 2**: the refined location computed in the previous step is used as starting location to perform a more precise relocation using cross-correlation to refine the differential travel times. If step1 fails, step2 is attempted anyway.

If step2 completes successfully the relocated origin is sent to the messaging system. 

### 2.1. Configuring a background catalog

The easiest choice is to use as background catalog the relocated multi-event results; the triplet *reloc-event.csv*, *phase.csv*, *station.csv*:

![Catalog selection option](/data/img/catalog-selection3.png?raw=true "Catalog selection from raw file format")

However, if the catalog is generated in XML format, it can be imported in the SeisComP database. In this case the background catalog can be a file containing just the origin ids. 

![Catalog selection option](/data/img/catalog-selection1.png?raw=true "Catalog selection from event/origin ids")

While it is neat to have the background catalog in the SeisComP database, this approach has few limitations. First it may take a lot of time for `scrtdd` to load a big catalog from the database comparing to loading it from files. Also, since the background catalog should be periodically updated, old events are continuosly updated with new origins, which can lead to a not optimal database performance-wise.

Once the background catalog is configured `scrtdd` can be enabled and started as any other SeisComP module.  New origins will be relocated as soon as they arrive in the messsaging system.

### 2.2 Testing

You might consider testing the configuration relocating some existing events to make sure the parameters are suitable for your use case. To test the real time relocation there are two command line options which relocate existing origins:

```
scrtdd --help

Mode:

  -O [ --origin-id ] arg                Relocate  the origin (or multiple 
                                        comma-separated origins) in 
                                        signle-event mode and send a message. 
                                        Each origin will be processed 
                                        accordingly to the matching profile 
                                        region unless the --profile option  is 
                                        used.
  --ep arg                              Event parameters XML file for offline 
                                        processing of contained origins 
                                        (implies --test option). Each contained
                                        origin will be processed in 
                                        signle-event mode unless 
                                        --reloc-catalog is provided, which 
                                        enable multi-event mode.

ModeOptions:

   --profile arg                        To be used in combination with other 
                                        options: select the profile 
                                        configuration to use

  --test                                Test mode, no messages are sent when 
                                        relocating a single event

  --xmlout                              Enable XML output when combined with 
                                        --reloc-catalog or --oring-id options

```

#### 2.2.1 Relocate origin ID and send the relocation to the messaging system for further processing

If we want to process an origin we can run the following command and then check on `scolv` the relocated origin (the messaging system must be active). This is mostly useful when we want to relocate an origin on a running system and keep the relocation:

```
scrtdd --origin-id someOriginId \
       --verbosity=3 --console=1 [db options] 
```

#### 2.2.2 Relocate origin ID but do not send the relocation (debug)

As above but add `--test` and the origin will not be sent to the messaging system. Useful for troubleshooting when the `scrtdd.saveProcessingFiles` option is enabled to verify the relocation files in `scrtdd.workingDirectory`.

```
scrtdd --origin-id someOriginId --test \
       --verbosity=3 --console=1 [db options]
``` 

#### 2.2.3 Relocate origin ID and store the result to XML file

Adding the `--xmlout` option allows to save the origin as a XML file. We can finally open the ile with `scolv` for inspection:

```
scrtdd --origin-id someOriginId --xmlout \
       --verbosity=3 --console=1 [db options] \
  >  relocated-origin.xml
```

#### 2.2.4 Relocate XML file and store the result to XML file

Similarly to other SeisComP commands the `--ep` option can be used for full offline processing. All origins contained in the input XML file are relocated.

```
scrtdd --ep origin.xml --verbosity=3 --console=1 [db options] \
  > relocated-origin.xml
```

### 2.3 Phase update

`scrtdd` uses cross-correlation to detect phases at stations with no associated picks in order to fix the pick time and uncertainty of automatic picks. Those features are especially useful in real-time to increase the quality and number of double-difference observations when automatic origins have only few picks/phases.

For automatic picks, the pick time is updated according to the average lag detected by all the good (above configured threshold) cross-correlation results. Since the real-time events are cross-correlated against catalog events, which have good manual picks, the updated pick time is expected to improve. The pick uncertainty is derived from the uncertainties of catalog-events. If no cross-correlation coefficients above the configured threshold are found, the pick is kept untouched.

For stations with no associated phases, `scrtdd` computes theoretical picks. Those are then cross-correlated against the catalog event ones. Every theoretical pick that has at least one good cross-correlation result is added to the relocated origin, with pick time and uncertainties derived from catalog phases (similarly to what is done for automatic picks). Those *good* theoretical picks are thus used in the double-difference system inversion. Theoretical picks that have no good cross-correlation results are simply discarded.

Picks that have been updated or created by `scrtdd` are identifiable by a `x` suffix (Px, Sx).

Manual picks are never modified.

### 2.4 Avoiding Relocation Loops

`scrtdd` listens and sends messages to the LOCATION group. In a default installation where the only locator is `scautoloc` that's not an issue: `scautoloc` will send an origin to LOCATION and `scrtdd` will receive it and send an updated origin to LOCATION.  However, when thare are multiple (re)locators (e.g. scanloc, screloc) that listen to LOCATION and send their own updated origin to LOCATION too, then an infinite loop happens! In this case a new messaging group needs to be created, e.g. RELOCATION, so that the origins flow from LOCATION to RELOCATION without going back.

 E.g. of a properly configured system:

```
                      LISTEN                       SEND 
              (MessagingSubscription)      (PrimaryMessagingGroup)
scautoloc             ...                        LOCATION
scanloc       LOCATION, ...                      LOCATION
screloc       LOCATION, ...                     RELOCATION
scrtdd        LOCATION, ...                     RELOCATION
scevent       LOCATION,RELOCATION, ...             ...
scamp         LOCATION,RELOCATION, ...             ...
scmag         LOCATION,RELOCATION, ...             ...
```

### Examples 

Below the single-event relocation of several manually reviewed origins.

![Single event relocation example picture](/data/img/singleEventRelocationExample.png?raw=true "Single Event Relocation example") 

The unit testing folder contains the code to generate some tests with synthetic data.

![Synthetic Data Relocation example picture](/data/img/singleEventRelocationSyntDataExample.png?raw=true "Synthetic Data Relocation example")

## 3. Cross-correlation

Good cross-correlation results are needed to achieve high quality double-difference observations, which in turn results in high resolution relocations. The purpose of the cross-correlation is to find the exact time difference between two picks of an event pair at a common station.

The default values of the cross-correlation parameters are meant to be a good starting point, however it might be useful to optimize them in certain scenarios.

![Relocation options](/data/img/xcorr.png?raw=true "Relocation options")

### 3.1 Eval-xcorr command

The `--eval-xcorr` command can be used to evaluate the cross-correlation parameter. 

`--eval-xcorr` allows to see:
* `#pha`: how many phases have been cross-correlated
* `pha good CC`: how many of those were successful (correlation coefficient above the configured threshold)
* `coeff`: the average correlation coefficient
* `goodCC/ph`: the average number of good matches per phase (each event phase appears in multiple double-difference observations, so multiple event pais, hence multiple cross-correlations)
* `time-diff`: the average pick time difference detected by the cross-correlation
* `+/-`: whenever sensible, it is also indicated the Mean Absolute Deviation of the value

It is especially interesting to compare the results before/after the catalog has been relocated. The new statistics should show better performance for events close to each other and gradually worsen with increasing inter-event distance. That is an indirect measure of the quality of the relocation as explained in Waldhauser & Ellsworth's paper.

```
scrtdd --eval-xcorr station.csv,event.csv,phase.csv --profile myProfile --verbosity=2 --console=1
```

Example output:

```
[...]
13:13:17 [warning] <FINAL STATS>
Cumulative stats: #pha 196006 pha good CC  72% coeff 0.72 (+/-0.09) goodCC/ph  9.9 (+/-4.2) time-diff [msec]  -0 (+/-52)
Cumulative stats P ph: #pha 118343 pha good CC  68% coeff 0.72 (+/-0.10) goodCC/ph  9.7 (+/-4.5) time-diff [msec]   0 (+/-52)
Cumulative stats S ph: #pha  77663 pha good CC  76% coeff 0.72 (+/-0.08) goodCC/ph 10.3 (+/-4.1) time-diff [msec]  -1 (+/-52)

Cross-correlated Phases by inter-event distance in 0.10 km step
 EvDist [km]  #Phases GoodCC AvgCoeff(+/-) GoodCC/Ph(+/-) time-diff[msec] (+/-)
 0.00-0.10      72667    73%  0.85 (0.09)    3.2 ( 2.2)       0 ( 29)
 0.10-0.20      85191    69%  0.81 (0.09)    2.7 ( 1.6)       0 ( 35)
 0.20-0.30      63659    61%  0.79 (0.09)    1.8 ( 0.9)       0 ( 41)
 0.30-0.40      46852    56%  0.77 (0.10)    1.6 ( 0.7)      -0 ( 44)
 0.40-0.50      54217    53%  0.76 (0.10)    1.5 ( 0.7)       2 ( 47)
 0.50-0.60      67184    55%  0.74 (0.09)    1.9 ( 0.9)       2 ( 46)
 0.60-0.70      51496    49%  0.74 (0.09)    1.5 ( 0.7)      -1 ( 48)
 0.70-0.80      36620    46%  0.73 (0.09)    1.3 ( 0.5)      -0 ( 51)
 0.80-0.90      30600    43%  0.73 (0.09)    1.2 ( 0.4)       0 ( 52)
 0.90-1.00      45866    45%  0.72 (0.09)    1.4 ( 0.6)       1 ( 53)
 1.00-1.10      44881    42%  0.72 (0.09)    1.4 ( 0.6)      -0 ( 53)
 1.10-1.20      34038    40%  0.72 (0.09)    1.3 ( 0.4)       1 ( 55)
 1.20-1.30      29119    38%  0.72 (0.09)    1.2 ( 0.4)       0 ( 57)
[...]
Cross-correlated Phases by event to station distance in 3.00 km step
StaDist [km]  #Phases GoodCC AvgCoeff(+/-) GoodCC/Ph(+/-) time-diff[msec] (+/-)
  0-3             134    84%  0.67 (0.06)    4.7 ( 3.0)      -6 ( 95)
  3-6            4616    87%  0.71 (0.07)   12.5 ( 8.1)      -0 ( 40)
  6-9           13307    84%  0.71 (0.07)   11.9 ( 7.2)       0 ( 35)
  9-12          16138    82%  0.71 (0.07)   12.5 ( 8.2)       1 ( 38)
 12-15          15743    81%  0.71 (0.07)   11.1 ( 6.9)      -1 ( 40)
 15-18          11340    78%  0.72 (0.08)   12.4 ( 8.2)      -0 ( 47)
 18-21           9874    75%  0.71 (0.07)   10.9 ( 7.1)      -0 ( 51)
 21-24          12193    74%  0.71 (0.07)   11.3 ( 7.1)      -0 ( 49)
 24-27          10537    73%  0.72 (0.08)   10.3 ( 6.6)      -1 ( 54)
 27-30          11503    75%  0.71 (0.07)   10.8 ( 6.6)      -3 ( 51)
[...]
Cross-correlations by station
Station       #Phases GoodCC AvgCoeff(+/-) GoodCC/Ph(+/-) time-diff[msec] (+/-)
4D.AG01.             2     0%  0.00 (0.00)    0.0 ( 0.0)       0 (  0)
4D.GDA01.           80    60%  0.80 (0.10)    2.8 ( 1.4)      -7 ( 27)
4D.GDA02.           34    94%  0.87 (0.07)    3.0 ( 0.9)       0 ( 27)
4D.GDA03.           70    80%  0.82 (0.11)    3.2 ( 1.2)      -1 ( 17)
4D.MH36.A          232    63%  0.76 (0.09)    3.3 ( 1.8)       8 ( 55)
4D.MH38.A            6    33%  0.67 (0.00)    1.0 ( 0.0)       0 (162)
4D.MH44.A          264    77%  0.77 (0.08)    3.1 ( 1.8)       9 ( 77)
4D.MH48.A          168    71%  0.78 (0.09)    2.9 ( 1.6)       1 ( 58)
4D.MH52.A           29    59%  0.84 (0.06)    2.4 ( 0.8)       2 ( 31)
4D.MH54.A           81    62%  0.79 (0.11)    1.7 ( 0.7)      -1 ( 44)
4D.RA41.            71    58%  0.81 (0.11)    3.0 ( 1.8)      -3 ( 42)
4D.RA42.             5    20%  0.64 (0.00)    1.0 ( 0.0)       2 (  0)
4D.RA43.            94    47%  0.83 (0.09)    3.2 ( 2.4)      -2 ( 43)
[...]
```

### 3.2 Cross-correlation parameters optimization

The `--eval-xcorr` option should be used to properly configure the cross-correlation parameters. The optimization process involves running `--eval-xcorr` with different configuration and analyzes the results. The goal is to have as many matches as possible (increase `GoodCC`) avoiding bad/false matches (very high values of `time-diff` are probably an indication of false matches): this is a trade-off.

The configuration parameters that are relevant for this analysis are:
* Cross-correlation coefficient threshold for P and S
* Cross-correlation window length for P and S
* Cross-correlation maximum lag/delay for P and S
* Signal to Noise ratio (SNR) configuration
* Waveform filter (and resampling)
* Clustering: number of neighbours

The SNR is particularly important to reject bad picks that lead to bad cross-correlations. The SNR signal/noise windows should be chosen so that they satisfies ALL the following 5 conditions:

* pick time too early (we want low SNR)
* pick time too late (we want low SNR)
* pick time perfect  (we want high SNR)
* pick time is early but acceptable (we want high SNR)
* pick time is late but acceptable (we want high SNR)

There are also few more parameters that are less relevant, but that might become important when relocating automatic origins (the automatic location might be very far from the final solution):
* Cross-correlation maximum station distance
* Cross-correlation maximum inter-event distance
* Clustering: maximum inter-event distance
* Clustering: station to inter-event distance ratio  

### 3.3 Logs

Some cross-correlation statistics are printed in both multi-event and single-event mode. Those can be seen in the log file or in the console output if `--console=1 --verbosity=3`) is used.

```
[info] Cross-correlation statistics: performed 40361, waveforms with Signal to Noise ratio too low 2435, waveforms not available 98
[info] Total xcorr 40361 (P 59%, S 41%) success 28% (11499/40361). Successful P 22% (5300/23844). Successful S 38% (6199/16517)
[info] xcorr on actual picks 24784/40361 (P 60%, S 40%) success 37% (9186/24784). Successful P 31% (4629/14761). Successful S 45% (4557/10023)
[info] xcorr on theoretical picks 15577/40361 (P 58%, S 42%) success 15% (2313/15577). Successful P 7% (671/9083). Successful S 25% (1642/6494)
```

The statistics are broken down in actual picks and theoretical picks. This is because `scrtdd` computes theoretical picks that are cross-correlated together with detected picks. This is useful to increase the number of double-difference observations. See the [Phase update](#24-phase-update) paragraph for further details. 

### 3.4 Waveforms inspection

A more in-depth source of information for waveform filtering and SNR options comes from this option:


```
scrtdd --help
  --debug-wf                            Enable saving of processed waveforms 
                                        into the profile working directory for 
                                        inspection.
```

Simply adding `--debug-wf` to the command line will make `scrtdd` dump to disk miniseed files for inspection (e.g. `scrttv` waveformfile.mseed). Just make sure to delete the folder before using this option to make sure to not look at previous relocation output. This option can be added to any `scrtdd` commands (e.g. `--relocate-profile`, `--ev`, `--origin-id` ) but it is mostly useful when relocating a single event mode because in multi-event mode there will be way too many waveforms to be able to check them all manually, although we can still do some random check to get an overall feeling of the filtering and SNR.

The generated waveforms will be stored in `workingDirectory/profileName/wfdebug/` (e.g. `~/seiscomp3/var/lib/rtdd/myProfile/wfdebug/`) after filtering and resampling, that is the same wavforms used for the cross-correlation. The file names follow the patterns:

* *evNumber.NET.STATION.phaseTime.manual.mseed*       (e.g. ev56.CH.SAYF2.S.manual.mseed)       - event 56 manual S phase on CH.SAYF2 station
* *evNumber.NET.STATION.phaseTime.automatic.mseed*    (e.g  ev4.CH.SAYF2.P.automatic.mseed)     - event 4 automatic P phase on CH.SAYF2 station
* *evNumber.NET.STATION.phaseTime.theoretical.mseed*  (e.g. ev267.CH.SFRU.Pt.theoretical.mseed) - event 267 theoretical P phase on CH.SFRU station
* *evNumber.NET.STATION.phaseTime.snr-rejected.mseed* (e.g. ev9.XY.VET02.S.snr-rejected.mseed)  - event 9 S phase not used in cross-correlation due to the configured signal to noise ratio threshold

Also, `scrtdd` logs tell us the details of the cross-correlation, so that we can see what waveform was cross-correlated with which others and then inspect the corresponding waveform miniseed files e.g.:

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

For comparison we can always find the raw waveforms (not processed) fetched from the configured recordStream and used as a cache in `workingDirectory/profileName/wfcache/` (e.g. `~/seiscomp3/var/lib/rtdd/myProfile/wfcache/`):
* `NET.ST.LOC.CH.startime-endtime.mseed`


## 4. Waveform data and recordStream configuration

SeisComP applications access waveform data through the RecordStream interface (see [official SeisComP documentation](https://www.seiscomp.de/doc/base/concepts/recordstream.html) for more details) and it is usually configured in *global.cfg* or passed via command line with `-I URI`. The RecordStream parameters define the service(s) through which SeisComP can access real-time and historical waveforms data (seedlink, fdsn, sds archive, [etc](https://www.seiscomp.de/doc/apps/global_recordstream.html)). A hypothetical RecordStream configuration might look like this:

```
recordstream = combined://slink/localhost:18000;sdsarchive//path/to/miniseed
```

This configuration is a combination of seedlink and sds archive, which allows `scrtdd` to retrieve catalog waveforms via sds and real-time event data via seedlink.

Please note that depending on the responsiveness of the seedlink server the real-time relocations may incur in delays. A couple of configuration options allow to control those delays: *timeout* and *retries*. The example below forces a timeout of 5 seconds (default is 5 minutes) and to not reconnect. In case of a timeout, `scrtdd` will proceed with the available data, without further delays:

```
recordstream = combined://slink/localhost:18000?timeout=5&retries=0;sdsarchive//path/to/miniseed
```


### 4.1 Waveforms data caching

Unless the recordStream points to a local disk storage, downloading waveforms might require a lot of time. For this reason `scrtdd` stores the waveforms to disk (called waveform cache) after downloading them. This applies only to the catalog event waveforms, which are used over and over again. That's not true for the real-time events, whose waveforms are used just once and never cached. The cache folder is `workingDirectory/profileName/wfcache/`.

However, for certain situations (e.g. debugging) it might be useful to cache all the waveforms, even the ones that are normally not cached. For those special cases the option --cache-wf-all can be used (stored in `workingDirectory/profileName/tmpcache/` which can be deleted afterwards).


### 4.2 Catalog waveforms preloading

When `scrtdd` starts for the first time it loads all the catalog waveforms and stores them to disk. However, if the option `performance.profileTimeAlive` is greater than 0, the catalog waveforms will be loaded only when needed (lazy loading) and not at start time. We can also force `scrtdd` to pre-download all waveforms using the following option:

```
scrtdd --load-profile-wf --profile myprofile
```

```
scrtdd --help
  --load-profile-wf                     Load catalog waveforms from the 
                                        configured recordstream and save them 
                                        into the profile working directory. Use
                                        in combination with --profile

  --profile arg                         To be used in combination with other 
                                        options: select the profile 
                                        configuration to use

```

## 5. Database connection

When SeisComP modules need to access the database for reading or writing data (events, picks, magnitudes, etc) they use the connection string configured in either `global.cfg` (which is inherited by every module) or in `scmaster.cfg`, in which case is scmaster module that passes the database connection string to every module when they connect to the messaging system (usually at module startup).

However, when running `scrtdd` from the command line, it doesn't connect to the messaging system and if the database connection is specified via `scmaster.cfg`, the information never reaches `scrtdd`. In this case the database connection must be passed as a command line option:

```
scrtdd [some options] -d  mysql://user:password@host/seiscompDbName
```

or in case of a Postgresql database:

```
scrtdd [some options] --plugins dbpostgresql -d postgresql://user:password@host/seiscompDbName
```

It is worth noting that this feature allows `scrtdd` to connect to remote seiscomp databases too and relocate events stored in other machines without interfering with the real-time processing happening there.
 

## 6. Troubleshooting

Log files are located in ~/.seiscomp/log/scrtdd.log. Alternatively, when running `scrtdd` from the command line, the following options can be used to see the logs on the console:

```
scrtdd [some options] --verbosity=3 --console=1
```

Verbosity 3 should be preferred to level 4, since the debug level 4 makes the logs hard to read due to the huge amount of information.

Also, enabling the `scrtdd.saveProcessingFiles` option makes `scrtdd` generates multiple information files inside `scrtdd.workingDirectory`. Those files can be useful for carefull inspections of the relocations.

### 6.1 Single-event

A typical *single-event* relocation log looks like the followig;

```
[info/RTDD] Performing step 1: initial location refinement (no cross correlation)
[info/RTDD] Selecting Neighbouring Events for event 1 lat 46.902294 lon 9.109304 depth 0.9287
[...]
Details of the Neighbouring Events found
[...]
[info/RTDD] Building and solving double-difference system...
[...]
Details of the solutions for each iteration of the solver
[...]
[info/RTDD] Step 1 relocation successful, new location: lat 46.899737 lon 9.111036 depth 1.3489 time 2020-10-29T20:08:36.572955Z
[info/RTDD] Relocation report:
            Origin changes: location=0.31[km] depth=0.42[km] time=-0.119[sec] 
            Rms change [sec]: 0.100 (before/after 0.415/0.516)
            Neighbours=70 Used Phases: P=13 S=20
            Stations distance [km]: min=4.5 median=36.6 max=65.7
            Neighbours mean distace to centroid [km]: location=6.97 depth=1.95
            Origin distace to neighbours centroid [km]: location=1.32 depth=-2.46 
            DD observations: 696 (CC P/S 0/0 TT P/S 285/411)
            DD observations residuals [msec]: before=-59+/-38.5 after=8+/-11.8

[info/RTDD] Performing step 2: relocation with cross correlation
[info/RTDD] Selecting Neighbouring Events for event 11371 lat 46.899737 lon 9.111036 depth 1.3489
[...]
Details of the Neighbouring Events found
[...] 
[info/RTDD] Computing cross-correlation differential travel times for event 10260
[...]
Details of cross-correlation
[...]  
[info/RTDD] Cross correlation performed 377, phases with Signal to Noise ratio too low 16, phases not available 0 (waveforms downloaded 0, waveforms loaded from disk cache 89)
[info/RTDD] Total xcorr 377 (P 46%, S 54%) success 71% (267/377). Successful P 59% (103/175). Successful S 81% (164/202)
[info/RTDD] Building and solving double-difference system...
[...]
Details of the solutions for each iteration of the solver
[...]
[info/RTDD] Step 2 relocation successful, new location: lat 46.899428 lon 9.110789 depth 1.5173 time 2020-10-29T20:08:36.558864Z
[info/RTDD] Relocation report:
            Origin changes: location=0.04[km] depth=0.17[km] time=-0.014[sec]
            Rms change [sec]: 0.038 (before/after 0.509/0.546)
            Neighbours=46 Used Phases: P=12 S=19 
            Stations distance [km]: min=8.7 median=36.9 max=65.7
            Neighbours mean distace to centroid [km]: location=5.40 depth=1.86
            Origin distace to neighbours centroid [km]: location=2.10 depth=-1.39 
            DD observations: 532 (CC P/S 103/164 TT P/S 117/148) 
            DD observations residuals [msec]: before=-59+/-38.5 after=8+/-18.3
[info/RTDD] Total Changes: location=0.35[km] depth=0.59[km] time=-0.133[sec] Rms=0.131[sec] (before/after 0.415/0.546)
```

### 6.2 Multi-event

A typical *multi-event* relocation log looks like the following:

```
[info/RTDD] Selecting Catalog Neighbouring Events 
[...]
Details of the Neighbouring Events selection [very very long]
[...]
[info/RTDD] Found 3 event clusters
[info/RTDD] Relocating cluster 1 (134 events)
[info/RTDD] Computing cross-correlation differential travel times for event 134
[...]
Details of cross-correlation [very very long]
[...]   
[info/RTDD] Cross-correlation performed 75917, phases with SNR ratio too low 1040, phases not available 12 (waveforms downloaded 0, waveforms loaded from disk cache 6325)
[info/RTDD] Total xcorr 75917 (P 62%, S 38%) success 87% (66110/75917). Successful P 81% (37825/46816). Successful S 97% (28285/29101)
[info/RTDD] Building and solving double-difference system...
[...]
Details of the solutions for each iteration of the solver
[...] 
[info/RTDD] Successfully relocated 134 events. RMS median 0.3288 [sec] median absolute deviation 0.0170 [sec]
[info/RTDD] Events RMS before relocation: median 0.3431 median absolute deviation 0.0318

[info/RTDD] Relocating cluster 2 (83 events)
[...]

[info/RTDD] Relocating cluster 3 (1583 events)
[...] 
```

### 6.3 Relocation statistics

`scrtdd` adds two comments to each relocated origin: `scrtddSourceOrigin` and `scrtddRelocationReport`. They can be both visualized in `scolv` (see official SeisComP documentation on how to do so), or they can be seen on the logs.

`scrtddSourceOrigin` contains the id of the origin that triggered the relocation. `scrtddRelocationReport` contains a summary of the relocation process. E.g.

```
Origin changes: location=0.23[km] depth=1.40[km] time=-0.147[sec]
Rms change [sec]: -0.153 (before/after 0.502/0.349)
Neighbours=80 Used Phases: P=37 S=16
Stations distance [km]: min=15.9 median=57.0 max=99.8
Neighbours mean distace to centroid [km]: location=5.11 depth=5.06
Origin distace to neighbours centroid [km]: location=1.30 depth=3.01
DD observations: 687 (CC P/S 141/47 TT P/S 375/124)
DD observations residuals [msec]: before=-106+/-21.6 after=9+/-26.2
```

To allow a comparison of the RMS before and after the relocation `scrtdd` computes the RMS before and after the relocation. Without that it would be hard to compare the RMSs. Each locator (scautoloc, scanloc, screloc, nonlinloc, scrtdd, etc) computes the RMS with a certain travel time table, that might not be the same as `scrtdd`. Moreover a locator might apply a specific logic on the RMS computation that prevents a comparison between different locators. For example NonLinLoc locator weighs the residuals by each pick weight and the wighting scheme is decided by NonLinLoc.

The following two lines can be a little cryptic to interpret: 
```
Neighbours mean distace to centroid [km]: location=5.11 depth=5.06
Origin distace to neighbours centroid [km]: location=1.30 depth=3.01
```

Their intent is to highlight how far the relocated event is to the neighbours centroid. The idea is that an event that falls within a cluster has a better chance to be properly relocated than an event that is far away from the neighbouring events.

The above information is also stored in the output files (events.csv,phases.csv,station,csv) of the multi-event relocation and it can be used to compute useful statistics for an entire catalog. Those are the column names containing the information:

```
startRms
locChange
depthChange
timeChange
numNeighbours
neigh_meanDistToCentroid
neigh_centroidToEventDist
neigh_meanDepthDistToCentroid
neigh_centroidToEventDepthDist
ph_usedP
ph_usedS
ph_stationDistMin
ph_stationDistMedian
ph_stationDistMax
ddObs_numTTp
ddObs_numTTs
ddObs_numCCp
ddObs_numCCs
ddObs_startResidualMedian
ddObs_startResidualMAD
ddObs_finalResidualMedian
ddObs_finalResidualMAD 
```

## 7. Custom velocity models

### 7.1 LOCSAT

In the `scrtdd` configuration it is possible to select any travel time table installed in SeisComP; this means the default SeisComP travel time tables and any other tables installed by the user. Although, this is a general SeisComP topic and we suggest to refer to the official SeisComP documentation, here is a quick recipe for generating your own travel time table from a custom velocity model.

SeisComP supports `LOCSAT` and `libtau` travel time table formats (1D velocity model). It is possible to generate a custom travel time table in `LOCSAT` format using the [TauP toolkit](https://www.seis.sc.edu/taup). 

First step is to have a velocity model in one of the formats supported by `TauP`. To do so, just make a copy of the SeisComP iasp91 or ak135 velocity models:

```
cp seiscomp_installation/share/ttt/iasp91.tvel mymodel.tvel
```
or
```
cp seiscomp_installation/share/ttt/ak135.tvel mymodel.tvel
```

Then edit the relevant layers to match your velocity model and leave the others untouched. Now run:

```
./TauP-installation/bin/taup_create -tvel mymodel.tvel --verbose
```

That generates `mymodel.taup` file that will be used by the subsequent TauP commands.

Before generating the travel time tables we need to decide the resolution (range of depths and distances). TauP uses a default range of depths and distances that is not very dense. So create a header file containing the desired depths and distances ranges. E.g. file `mymodel.header`

```
n # P,S travel-time tables resolution for mymodel
30     # number of depth samples
   0.00   0.50   1.00   1.50   2.50   3.00   3.50   4.00   4.50   5.00
   6.00   7.00   8.00   9.00  10.00  15.00  20.00  25.00  30.00  40.00
  50.00  75.00 100.00 150.00 200.00 300.00 400.00 500.00 600.00 800.00
130    # number of distances
   0.00   0.05   0.10   0.15   0.20   0.25   0.30   0.35   0.40   0.45
   0.50   0.60   0.70   0.80   0.90   1.00   1.10   1.20   1.30   1.40
   1.50   1.60   1.70   1.80   1.90   2.00   2.50   3.00   4.00   5.00
   5.50   6.00   6.50   7.00   7.50   8.00   8.50   9.00   9.50  10.00
  11.00  12.00  13.00  14.00  15.00  16.00  17.00  18.00  19.00  20.00
  21.00  22.00  23.00  24.00  25.00  26.00  27.00  28.00  29.00  30.00
  31.00  32.00  33.00  34.00  35.00  36.00  37.00  38.00  39.00  40.00
  41.00  42.00  43.00  44.00  45.00  46.00  47.00  48.00  49.00  50.00
  51.00  52.00  53.00  54.00  55.00  56.00  57.00  58.00  59.00  60.00
  61.00  62.00  63.00  64.00  65.00  66.00  67.00  68.00  69.00  70.00
  71.00  72.00  73.00  74.00  75.00  76.00  77.00  78.00  79.00  80.00
  81.00  82.00  83.00  84.00  85.00  86.00  87.00  88.00  89.00  90.00
  95.00 100.00 105.00 110.00 115.00 120.00 125.00 130.00 135.00 140.00
```

Finally we can generate the travel time tables:

```
./TauP-installation/bin/taup_table -mod mymodel -ph ttp+ -locsat -header mymodel.header -o mymodel.P
./TauP-installation/bin/taup_table -mod mymodel -ph tts+ -locsat -header mymodel.header -o mymodel.S
./TauP-installation/bin/taup_table -mod mymodel -ph PcP  -locsat -header mymodel.header -o mymodel.PcP
./TauP-installation/bin/taup_table -mod mymodel -ph Pg   -locsat -header mymodel.header -o mymodel.Pg
./TauP-installation/bin/taup_table -mod mymodel -ph Pn   -locsat -header mymodel.header -o mymodel.Pn
./TauP-installation/bin/taup_table -mod mymodel -ph pP   -locsat -header mymodel.header -o mymodel.pP
./TauP-installation/bin/taup_table -mod mymodel -ph PP   -locsat -header mymodel.header -o mymodel.PP
./TauP-installation/bin/taup_table -mod mymodel -ph pS   -locsat -header mymodel.header -o mymodel.pS
./TauP-installation/bin/taup_table -mod mymodel -ph ScP  -locsat -header mymodel.header -o mymodel.ScP
./TauP-installation/bin/taup_table -mod mymodel -ph Sg   -locsat -header mymodel.header -o mymodel.Sg
./TauP-installation/bin/taup_table -mod mymodel -ph Sn   -locsat -header mymodel.header -o mymodel.Sn
./TauP-installation/bin/taup_table -mod mymodel -ph sP   -locsat -header mymodel.header -o mymodel.sP
```

Last step is to copy the travel time tables to the SeisComP installation folder so that all modules can see the new model:

```
cp mymodel* seiscomp_installation/share/locsat/tables/
```

![LOCSAT TTT](/data/img/locsat-ttt.png?raw=true "LOCSAT TTT")


### 7.2 NonLinLoc

Please refer to [NonLinLoc by Anthony Lomax](<http://alomax.free.fr/nlloc/>) documentation on how to generate grid files. Once you have them you can configure in `scrtdd` in travel time table options.

The following geographic transformations (TRANS statement) are currently supported: GLOBAL 2D, SIMPLE 2D and 3D, SDS 2D and 3D. Also both float and double values are supported as well as byte swapping.

![NLL TTT](/data/img/nll-ttt.png?raw=true "NLL TTT")


## 8. Scolv Locator plugin

A (re)locator plugin is also available in the code, which makes `scrtdd` available via `scolv`. To enable this plugin just add `rtddloc` to the list of plugins in the global configuration.

![Locator plugin](/data/img/locator-plugin.png?raw=true "Locator plugin")

Please note that this plugin is not strictly required since `scrtdd` would relocated any manaul origins anyway (if configured to do so) and the relocated origin will appear on `scolv` as soon as ready.

Also `scolv` doesn't allow to create new picks when performing a relocation, so `scrtdd` plugin disable the cross-correlation on theoretical picks since those picks will not be reported on `scolv`.


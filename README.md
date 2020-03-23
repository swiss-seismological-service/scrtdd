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
 *                                                                         *
 *   Developed by Luca Scarabello, Tobias Diehl                            *
 *                                                                         *
 ***************************************************************************/
</pre>
# Description

The SCRTDD is a [Seiscomp3](<https://github.com/SeisComP3/seiscomp3>) extension module that implements both Real-Time Double-Difference Event relocation and classic offline multi-event (catalog) Double-Difference relocation.

The actual methods are described in the paper "Near-Real-Time Double-Difference Event Location Using Long-Term Seismic Archives, with Application to Northern California" by Felix Waldhauser and "A Double-Difference Earthquake Location Algorithm: Method and Application to the Northern Hayward Fault, California" by Felix Waldhauser et al.

The actual Double-Difference inversion method is currently performed by the hypoDD software, which has to be installed in the system (currently tested on version 1.3 and 2.1b), while scrtdd code deals with the seiscomp3 internal details, allowing the user to quickly and easily perform Double-Difference relocations on existing events or in real-time.


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

## 1. Install hypoDD software

Go to https://www.ldeo.columbia.edu/~felixw/hypoDD.html and dowload hypoDD. Untar the sources, go to src folder, run make. Done! 

Please remember to set the Hypodd array limits (compilation time options available via include/*inc files) accordingly with the size of your problem.


## 2. Defining a background catalog for real time relocations

New origins will be relocated in real time against a background catalog of high quality locations. Those high quality events that form the background catalog can be already present in seiscompo database or not. In the latter case the background catalog has to be generated. Either way, the catalog has to be specified in scrtdd when configuring it.

### 2.1. Selecting a background catalog from existing origins or events

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

### 2.2. Generating a background catalog

In the more general case in which we don't already have a background catalog, we have to build one using scrtdd. Generating a background catalog involves two steps: first we select the candidate events that might have low quality locations (we still want them to be as accurate as possible though, and with many phases. That is we want origins manully relocated using manual picks) and then we use scrtdd to relocated those events to achieve the quality needed for a background catalog. The higher the quality of the background catalog, the higher will be the quality of the real time origin relocations.

#### 2.2.1 Dumping the candidate events to files

Once the candidate origins have been selected, we write their seiscomp id in a file using the same format as the example above (myCatalog.csv). Once that's done, we run the command:

```
scrtdd --dump-catalog myCatalog.csv [db and log options]
```

See also the `Troubleshooting` paragraph for details on how to specify a database connection via command line and for printing the logs on the console for easy debugging.

The above command will generate three files (event.csv, phase.csv and stations.csv) which contain the information needed by scrtdd to relocate the selected events. 

Note: those file use the extended catalog format, useful when the catalog information is fully contained in those files. Compare to the previous example format (list of event/origin ids) the extended format doesn't require the catalog to be present on the seiscomp database. This might come in handy when using events coming from a different seiscomp installation e.g the events can be dumped (--dump-catalog) and used in another machine. Also those file can be easily generated from another source.

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
id,latitude,longitude,elevation,networkCode,stationCode
4DAG01,46.457412,8.079460,2358.0,4D,AG01
4DAG02,46.460620,8.078122,2375.0,4D,AG02
4DAG03,46.458288,8.075408,2369.0,4D,AG03
4DBSG1,46.107760,7.732020,3378.0,4D,BSG1
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

Now that we have dumped the events (event.csv, phase.csv, stations.csv) we might perform some editing of those files, if required, then we relocate them. 

#### 2.2.2 Relocating the candidate events (multi-event mode) to generate a high quality background catalog

To relocate a catalog (multi-event relocation) we need to create a new profile in scrtdd configuration and then we have to set the generated files (event.csv, phase.csv, stations.csv) as the catalog of the profile.

![Catalog selection option](/data/catalog-selection2.png?raw=true "Catalog selection from raw file format")

At this point we have to configure the other profile options that control the relocation process: `step2options` (`step1options` is used only in real-time), which control the creation of catalog absolute travel time entries (dt.ct file in HypoDD terminology) and cross correlated differential travel times for pairs of events (dt.cc file in HypoDD terminology). The clustering options should work just fine with the default values, however some adjustments are worth it. For example:

![Relocation options](/data/multiEventStep2options.png?raw=true "Relocation options")

Then it is time to set the cross-correlation parameters, which require a more careful selection and it is covered in the next paragraph.

In the ```data``` folder of this project there are some hypoDD 2.1b configuration example files (the velocity model has to be changed to reflect your specific use case).

Finally, when the configuration is ready, we can relocate the catalog with the command:

```
scrtdd --reloc-profile profileName [db and log options] 
```

scrtdd will relocated the catalog and will generate another set of files reloc-event.csv reloc-phase.csv and reloc-stations.csv, which together define a new catalog with relocated origins. At this point we should check the relocated events (and HypoDD logs) and see if we are happy with the results. If not, we change scrtdd settings and relocate the catalog again until we are satisfied with the locations.

As an example you can see below two catalogs before and after scrtdd relocation:


![Relocation example picture](/data/multiEventRelocationExample.png?raw=true "Relocation example")

The relocation output files (reloc-event.csv reloc-phase.csv and reloc-stations.csv) will become the background catalog used in real-time relocation. At this point the profile configuration is not needed anymore, so the profile can be removed from the list of active profiles (`scrtdd.activeProfiles`) or it can be updated with real-time configuration.

![Catalog selection option](/data/catalog-selection3.png?raw=true "Catalog selection from raw file format")

Now that we have a high quality background catalog we are ready to perform real time relocation.

Be aware that the first time you run the command it will be very slow because the waveforms have to be downloaded from the configured recordstream and saved to disk for the next runs, which will be much faster. However some temporary waveforms (e.g. theoretical phases) are never saved to disk since they are not part of the catalog phases and thus not useful when `scrtdd` is enabled for real-time relocations. Nevertheless, during the parameters tuning phase the relocation is performed several times and storing the temporary waveforms too makes the process faster. For this reason we can use an option `--debug-wf-cache` that force `scrtdd` to save all waveforms to disk. Just be aware that using this options increase the size of the cache folder (`workingDirectory/profileName/wfcache/` e.g. `~/seiscomp3/var/lib/rtdd/myProfile/wfcache/`) with phases that will not be useful in real-time. so you might want to delete the cache folder at the end of the relocation attempts.

#### 2.2.3 Cross-correlation, waveform filtering and signal to noise ratio options

Those parameters require some trial and error to be correctly set but after that they can be kept identical for real-time relocation without any major change, except for the parameter ```maxDelay``` that should be increased in real-time relocation to accomodate the larger uncertainty of automatic picks compared to manual ones (unless scrtdd is used to only relocate manually reviewed origins).

![Relocation options](/data/xcorr.png?raw=true "Relocation options")

To help figuring out the right values for cross-correlation the log file (or console output with ```--console=1 --verbosity=3```) come in handy:

```
12:32:22 [info] Cross-correlation statistics: performed 40361, waveforms with Signal to Noise ratio too low 2435, waveforms not available 98
12:32:22 [info] Total xcorr 40361 (P 59%, S 41%) success 28% (11499/40361). Successful P 22% (5300/23844). Successful S 38% (6199/16517)
12:32:22 [info] xcorr on actual picks 24784/40361 (P 60%, S 40%) success 37% (9186/24784). Successful P 31% (4629/14761). Successful S 45% (4557/10023)
12:32:22 [info] xcorr on theoretical picks 15577/40361 (P 58%, S 42%) success 15% (2313/15577). Successful P 7% (671/9083). Successful S 25% (1642/6494)
```

scrtdd computes theoretical picks so that they can be used in cross-correlation. If the resulting correlation coefficient is good enough a new pick is also created and added to the relocated origin.

Another source of information for waveform filtering and signal to noise ratio options comes from a command line option:


```
scrtdd --help
  --debug-wf                            Enable the saving of waveforms 
                                        (filtered/resampled, SNR rejected, ZRT 
                                        projected and scrtdd detected phase) 
                                        into the profile working directory. 
```

Simply adding ```--debug-wf``` to the command line will make scrtdd dump to disk miniseed files for inspection (e.g. scrttv waveformfile.mseed). Just make sure to delete the folder before using this option to make sure to not look at previous relocation output. The option is mostly useful when relocating a single event mode because in multi-event mode there will be way too many waveforms to be able to check them manually, but we can still do some random check to get an overall feeling of the filtering and snr.

The generated waveforms will be stored in `workingDirectory/profileName/wfdebug/` (e.g. `~/seiscomp3/var/lib/rtdd/myProfile/wfdebug/`) after filtering and resampling, that is the same wavforms used for the cross-correlation. The files follow the patters:

* evNumber.NET.STATION.phaseTime.manual.mseed       (e.g. ev56.CH.SAYF2.S.manual.mseed)       - event 56 manual S phase on CH.SAYF2 station
* evNumber.NET.STATION.phaseTime.automatic.mseed    (e.g  ev4.CH.SAYF2.P.automatic.mseed)     - event 4 automatic P phase on CH.SAYF2 station
* evNumber.NET.STATION.phaseTime.theoretical.mseed  (e.g. ev267.CH.SFRU.Pt.theoretical.mseed) - event 267 theoretical P phase on CH.SFRU station
* evNumber.NET.STATION.phaseTime.snr-rejected.mseed (e.g. ev9.XY.VET02.S.snr-rejected.mseed)  - event 9 S phase not used in cross-correlation due to the configured signal to noise ratio threshold

The logs tell us the details of the cross-correlation, so that we can see what waveform was cross-correlated with which others e.g.

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


#### 2.2.4 Using ph2dt

It is possible to use ph2dt utility to perform the clustering. It this case the scrtdd configuration `step2options.clusteing.*` will not be used. Instead, ph2dt will be run to generate dt.ct file, and for each entry in the generated dt.ct file the cross-correlation will be performed and the relative dt.cc file created.

```
scrtdd --reloc-profile profileName --use-ph2dt /some/path/ph2dt.inp [--ph2dt-path /some/path/ph2dt]
```

#### 2.2.5 Useful options

There are also some other interesting catalong related options:


```
scrtdd --help

Catalog:
  --dump-catalog arg                    Dump the seiscomp event/origin id file 
                                        passed as argument into a catalog file 
                                        triplet (station.csv,event.csv,phase.cs
                                        v)
  --dump-catalog-xml arg                Convert the input catalog into XML 
                                        format. The input can be a single file 
                                        (containing seiscomp event/origin ids) 
                                        or a catalog file triplet 
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
  --dump-catalog-options arg            Allows --dump-catalog to accept event 
                                        ids besides origin ids. For each event 
                                        id an origin will be selected following
                                        the provided options whose format is: 
                                        'type,evalmode,includeCreator,excludeCr
                                        eator,profile', where 
                                        type:'preferred','last','first' 
                                        evalmode:'any','onlyManual','onlyAutoma
                                        tic' includeCreator:'any' or 
                                        author/methodID  excludeCreator:'none' 
                                        or author/methodID profile:'any' or 
                                        profileName. e.g. to select preferred 
                                        origins from the provided event ids use
                                        'preferred,any,any,none,any'

MultiEvents:
  --reloc-profile arg                   Relocate the catalog of profile passed 
                                        as argument
  --ph2dt-path arg                      Specify path to ph2dt executable
  --use-ph2dt arg                       When relocating a catalog use ph2dt. 
                                        This option requires a ph2dt control 
                                        file
  --no-overwrite                        When relocating a profile don't 
                                        overwrite existing files in the working
                                        directory (e.g. avoid re-computation of
                                        cross correlation file dt.cc and allow
                                        manual editing of hypoDD.inp configuration)

```

## 3. Real time single origin relocation

We first need to create a new profile or we can modify the profile used to create the background catalog. If we decice to go for a new profile we might want to copy the waveform cache folder of the background catalog profile to the new real-time profile cache folder to avoid re-downloading all the catalog phases.

The waveform cache folder is located in `workingDirectory/profileName/wfcache/` (e.g. `~/seiscomp3/var/lib/rtdd/myProfile/wfcache/`).

## 3.1 Testing

Real time relocation uses the same configuration we have seen in full catalog relocation, but real time relocation is done in two steps, each one controlled by a specific configuration.

Step 1: location refinement. In this step scrtdd performs a preliminary relocation of the origin using only catalog absolute travel time entries (dt.ct only).

Step 2: the refined location computerd in the previous step is used as starting location to perform a more precise relocation using both catalog absolute travel times (dt.ct) and differential travel times from cross-correlation (dt.cc). If step1 fails, step2 is attempted anyway.

If step2 completes successfully the relocated origin is sent to the messaging system.


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
scrtdd -O someOriginId [db and log options] 
```

For testing purpose we are more likely interested in not interfere with database and messaging system so we can use the `--ep` option and the relocated origin will be sent to the output as xml file. We can finally open the xml file with scolv for inspection:

```
scrtdd -O someOriginId --ep - [db and log options] >  relocated-origin.xml
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

### 3.2 Tesing speed-up

Since real-time origin waveforms are never saved to disk, if we are testing over and over the same origins to find the best configuration parameters we might consider using `--debug-wf-cache` option that force `scrtdd` to save all waveforms to disk. Just be aware that using this options increase the size of the cache folder (`workingDirectory/profileName/wfcache/` e.g. `~/seiscomp3/var/lib/rtdd/myProfile/wfcache/`) with phases that will not be useful in real-time (only the catalog phases will be used again), so you might want to delete the cache folder at the end of the relocation testing (or copy it beforehand and restore it).

### 3.3 Catalog waveforms preloading

When scrtdd starts for the first time it loads all the catalog waveforms and store them to disk. In this way the waveforms become quickly available in real-time without the need to access the recordstream. However this takes some time.  If for some reasons (debugging?) we need to force the downloading of all waveforms before starting the module we can use the following option (once the waveforms are downloaded they will not be downloaded again unless the files are deleted from disks):

```
scrtdd --help
  --load-profile-wf arg                 Load catalog waveforms from the 
                                        configured recordstream and save them 
                                        into the profile working directory
```

### 3.4 RecordStream configuration

SeisComP3 applications access waveform data through the RecordStream interface and it is usually configured in global.cfg, for example:

```
recordstream = combined://slink/localhost:18000;sdsarchive//mnt/miniseed
```

This configuration is a combination of seedlink and sds archive, which is very useful because scrtdd catalog waveforms can be retrieved via sds while real-time event data can be fetched via seedlink (much faster since recent data is probably already in memory).

Please note that seedlink service might delay the relocation of incoming origins due to timeouts. For this reason we suggest to pass to configuration option: timeout and retries

e.g. Below we force a timeout of 2 seconds(default is 5 minutes) and do not try to reconnect (`scrtdd` will deal with what data is available. We can also configure the `cron.delayTimes` option to re-perform the relocation some minutes later in case we know more waveforms will become available later):

```
recordstream = combined://slink/localhost:18000?timeout=2&retries=0;sdsarchive//rz_nas/miniseed
```

### 3.5 Performance

scrtdd spends most of the relocation time downloading waveforms (real-time events, the catalog waveforms are cached to disk) and the rest is shared betweeen cross-correlation and hypoDD execution. Two configuration options have a huge impact on the performance: `step2options.clustering` and `crosscorrelation.s-phase.components`. `step2options.clustering` is relevant because we can specify how many neighbouring events and how many phases we want to use, which consequently determines the number of cross-correlations to perform and the size of the input for hypodDD. `crosscorrelation.s-phase.components` defines the components we want to use in the cross-correlation of the `S` phases. Usign the `T` components means scrtdd has to download the additional two components to perform the projection to ZRT.

If we want to analyze the performance impact of downloading waveforms during the relocatione we can make use of `--debug-wf-cache` option: the waveforms will be cached the first time we relocate an origin and the second time we relocate the same origin they will be read from disk.

## 4. Locator plugin

A (re)locator plugin is also avaiable in the code, which makes scrtdd available via scolv. To enable this plugin just add `rtddloc` to the list of plugins in the global configuration.

![Locator plugin](/data/locator-plugin.png?raw=true "Locator plugin")

Please note that this plugin is not strictly required since `scrtdd` would relocated any manaul origins anyway (if configured to do so) and the relocated origin will appear on `scolv` as soon as ready.

Also scolv doesn't allow to create new picks when performing a relocation, so `scrtdd` plugin disable the cross-correlation on theoretical picks since those picks will not be reported on scolv.

## 5. Troubleshooting

### 5.1. Logs

Check log file: ~/.seiscomp/log/scrtdd.log 

Alternatively, when running scrtdd from the command line use the following options to see the logs on the console:

```
scrtdd [some options]--verbosity=3 --console=1
```
### 5.2. Database connection

The seiscomp database connection can be configured either inside `global.cfg`  (and thus every module knows what database to use since they inherit global.cfg), or in `scmaster.cfg`, in which case scmaster module passes the database connection string to every module when they connect to the messagin system. Since several scrtdd command line options don't need the messaging system,  scrtdd doesn't connect to it and in those cases we have to pass the database connection string to scrtdd via command line option (Since the database is still required for the inventory).
Also, specifying the database connection via command line is useful to use a local seiscomp installation to access events stored in another seiscomp installation. We can use the -d option to specify the database connection, e.g.

```
scrtdd [some options] -d  mysql://user:password@host/seiscompDbName
```

or in case of a Postgresql database:

```
scrtdd [some options] --plugins dbpostgresql -d postgresql://user:password@host/seiscompDbName
``` 

### 5.3. Working directory and HypoDD input/output files

A useful option we can find in scrtdd configuration is `keepWorkingFiles`, which prevent the deletion of scrtdd processing files from the working directory (e.g. `~/seiscomp3/var/lib/rtdd/myProfile/`). In this way we can access the working folder and check input, output files used for running hypodd. More importantly we can also run hypodd from the command line using the same files generated by `scrtdd` (and possible edit those) and view the console output since the hypodd log file doesn't always reports all the errors. Look at `scrtdd` logs to understand where the working directory is and how to run hypodd manually:

```
[...]
20:51:33 [info] Creating station file ~/seiscomp3/var/lib/rtdd/myProfile/20170103161341_46266_007400_20200225195133_0950/step1/station.dat
20:51:33 [info] Creating event file ~/seiscomp3/var/lib/rtdd/myProfile/20170103161341_46266_007400_20200225195133_0950/step1/event.dat
20:51:33 [info] Creating differential travel time file ~/seiscomp3/var/lib/rtdd/myProfile/20170103161341_46266_007400_20200225195133_0950/step1/dt.ct
20:51:33 [info] Running hypodd...
20:51:33 [info] Working directory ~/seiscomp3/var/lib/rtdd/myProfile/20170103161341_46266_007400_20200225195133_0950/step1
20:51:33 [info] Executing command: /bin/sh -c ~/dev/HYPODD_2.1b/src/hypoDD/hypoDD hypoDD.inp >hypoDD.out 2>&1 
20:51:33 [info] Loading catalog relocated by hypodd...
[...]
```

### 5.4. HypoDD limits

Finally, remember to set the Hypodd array limits (compilation time options available via *inc files) accordingly with the size of your problem. If the full catalog relocation doesn't seem to relocate at all, you might have probably hit array limits.

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
scrtdd --dump-catalog myCatalog.csv
```

It is usually convenient to see the logs on the console to detect possible errors. To achieve this we can add the options ```--verbosity=3 --console=1``` to the command.

Depending on your configuration, the seiscomp database could be configured inside ```global.cfg```  and thus every module knows what database to use, or this option can be provided to the modules via scmaster (```scmaster.cfg```) once they connect to the messagin system. In the latter case we have to pass the database connection to scrtdd via command line option since ```--dump-catalog``` detach scrtdd from the messaging. Specifying the database connection via command line is also useful when you want to use your local seiscomp installation to access events that reside on a different host (where a seiscomp database is installed). Let's use the -d option to specify the database connection then:

```
scrtdd --dump-catalog myCatalog.csv  -d  mysql://user:password@host/seiscompDbName --verbosity=3 --console=1
```

or in case you have a Postgresql database:

```
scrtdd --dump-catalog myCatalog.csv --plugins dbpostgresql -d postgresql://user:password@host/seiscompDbName --verbosity=3 --console=1
```

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

Then it is time to set the cross correlation parameters, which require a more careful selection and it is covered in the next paragraph.

In the ```data``` folder of this project there are some hypoDD 2.1b configuration example files (the velocity model has to be changed to reflect your specific use case).

Finally, when the configuration is done, we can relocate the catalog with the command:

```
scrtdd --reloc-profile profileName --verbosity=3 --console=1
```

scrtdd will relocated the catalog and will generate another set of files reloc-event.csv reloc-phase.csv and reloc-stations.csv, which together define a new catalog with relocated origins. At this point we should check the relocated events and see if we are happy with the results. If not, we change scrtdd settings and relocate the catalog again until we are satisfied with the locations. Be aware that the first time you run the command it will be very slow because the waveforms have to be downloaded from the configured recordstream (but they will be saved to disk for the next runs which will be much faster).

 As an example you can see below two catalogs before and after scrtdd relocation:


![Relocation example picture](/data/multiEventRelocationExample.png?raw=true "Relocation example")


The command output files (reloc-event.csv reloc-phase.csv and reloc-stations.csv) will become the background catalog used in real-time relocation. At this point the profile configuration is not needed anymore, so the profile can be removed from the list of active profiles ('activeProfiles') or it can be updated with real-time configuration.

![Catalog selection option](/data/catalog-selection3.png?raw=true "Catalog selection from raw file format")

Now that we have a high quality background catalog we are ready to perform real time relocation.

#### 2.2.3 Cross correlation, waveform filtering and signal to noise ratio options

Those parameters require some trial and error to be correctly set but after that they can be kept identical for real-time relocation without any major change, except for the parameter ```maxDelay``` that should be increased in real-time relocation to accomodate the larger uncertainty of automatic picks compared to manual ones (unless scrtdd is used to only relocate manually reviewed origins).

![Relocation options](/data/xcorr.png?raw=true "Relocation options")

To help figuring out the right values for cross correlation the log file (or console output with ```--console=1 --verbosity=3```) come in handy:

```
12:32:22 [info] Cross correlation statistics: performed 40361, waveforms with Signal to Noise ratio too low 2435, waveforms not available 98
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

The logs tell us the details of the cross correlation, so that we can see what waveform was cross-correlated with which others e.g.

```
15:12:02 [info] xcorr: event   267 on   CH BERNI phase Pt - no good cross correlations pairs
15:12:02 [info] xcorr: event   267 on   CH BIBA phase Pt - no good cross correlations pairs
15:12:02 [info] xcorr: event   267 on   CH BIBA phase St - no good cross correlations pairs
15:12:02 [info] xcorr: event   267 on   CH BNALP phase Pt - average correlation coefficient 0.79 over 1 close-by events (206 )
15:12:02 [info] xcorr: event   267 on   CH BNALP phase St - no good cross correlations pairs
15:12:02 [info] xcorr: event   267 on   CH BOURR phase Pt - average correlation coefficient 0.81 over 2 close-by events (206 265 )
15:12:02 [info] xcorr: event   267 on   CH BRANT phase Pt - no good cross correlations pairs
15:12:02 [info] xcorr: event   267 on   CH DAGMA phase Pt - no good cross correlations pairs
15:12:02 [info] xcorr: event   267 on   CH DAVOX phase Pt - no good cross correlations pairs
15:12:02 [info] xcorr: event   267 on   CH  DIX phase St - average correlation coefficient 0.81 over 10 close-by events (1 25 37 41 43 76 141 205 206 223 )
15:12:02 [info] xcorr: event   267 on   CH EMBD phase Pt - no good cross correlations pairs
15:12:02 [info] xcorr: event   267 on   CH EMBD phase St - no good cross correlations pairs
```

For comparison we can always find the raw waveforms (not processed) fetched from the configured recordstream and used as a cache in `workingDirectory/profileName/wfcache/` (e.g. `~/seiscomp3/var/lib/rtdd/myProfile/wfcache/`):
* `NET.ST.LOC.CH.startime-endtime.mseed`


#### 2.2.4 Using ph2dt

It is possible to use ph2dt utility to perform the clustering. It this case the scrtdd configuration `step2options.clusteing.*` will not be used. Instead, ph2dt will be run to generate dt.ct file, and for each entry in the generated dt.ct file the cross correlation will be performed and the relative dt.cc file created.

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

Real time relocation uses the same configuration we have seen in full catalog relocation, but real time relocation is done in two steps, each one controlled by a specific configuration.

Step 1: location refinement. In this step scrtdd performs a preliminary relocation of the origin using only catalog absolute travel time entries (dt.ct only).

Step 2: the refined location computerd in the previous step is used as starting location to perform a more precise relocation using both catalog absolute travel times (dt.ct) and differential travel times from cross correlation (dt.cc). If step1 fails, step2 is attempted anyway.

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
scrtdd -O someOriginId --verbosity=3 --console=1
```

For testing purpose we are more likely interested in not interfere with database and messaging system so we can use the `--ep` option and save the xml to file, finally open the file with scolv for inspection:

```
scrtdd -O someOriginId --ep - --verbosity=3 --console=1 >  relocated-origin.xml
```

Alternatively the `--ep` option can process origins contained in an XML file:

```
scxmldump -fPAMF -p -O originId -o origin.xml --verbosity=3  --console=1

scrtdd --ep origin.xml --verbosity=3 --console=1 [db option] > relocated-origin.xml
```

Also, as explained in the , we can use the ```--debug-wf``` option to help debugging.

As an example you can see below the single event relocation of several manually reviewed origins (when relocating automatic origins the quality and number of relocated origins is certainly lower).

![Single event relocation example picture](/data/singleEventRelocationExample.png?raw=true "Single Event Relocation example")


Once we are happy with the configuration we can simply enable and start scrtdd as any other module. Please note that when scrtdd starts for the first time it will load all the catalog waveforms and store them to disk, to make them available in real time without the need to access the recordstream. However this takes some time. You can also force the loading of all waveforms before starting the module using the following option:

```
scrtdd --help
  --load-profile-wf arg                 Load catalog waveforms from the 
                                        configured recordstream and save them 
                                        into the profile working directory
```

## 4. RecordStream configuration

SeisComP3 applications access waveform data through the RecordStream interface and it is usually configured in global.cfg, for example:

```
recordstream = combined://slink/localhost:18000;sdsarchive//mnt/miniseed
```

This configuration is a combination of seedlink and sds archive, which is very useful because scrtdd catalog waveforms can be retieved via sds while reat time event data can be fetched via seedlink (much faster since recent data is probably already in memory).

The seedlink service can sometime delay the relocation of incoming origin due to timeouts. For this reason we suggest to pass to configuration option: timeout and retries

e.g. Here we force a timeout of 2 seconds(default is 5 minutes) and do not try to reconnect. 

```
recordstream = combined://slink/localhost:18000?timeout=2&retries=0;sdsarchive//rz_nas/miniseed
```

## 5. Locator plugin

A (re)locator plugin is also avaiable in the code, which makes scrtdd available via scolv. To enable this plugin just add `rtddloc` to the list of plugins in the global configuration.

![Locator plugin](/data/locator-plugin.png?raw=true "Locator plugin")

Please note that this plugin is not strictly required since `scrtdd` would relocated any manaul origins anyway (if configured to do so) and the relocated origin will appear on `scolv` as soon as ready.

Also scolv doesn't allow to create new picks when performing a relocation, so `scrtdd` plugin disable the cross correlation on theoretical picks since those picks will not be reported on scolv.

## 6. Troubleshooting

Check log file: ~/.seiscomp/log/scrtdd.log 

Alternatively, when running scrtdd from the command line use the following options:

```
# set log level to info (3), or debug (4) and log to the console (standard output) insted of log file
--verbosity=3 --console=1
```

A useful option we can find in scrtdd configuration is `keepWorkingFiles`, which prevent the deletion of scrtdd processing files from the working directory (e.g. `~/seiscomp3/var/lib/rtdd/myProfile/). In this way we can access the working folder and check input, output files used for running hypodd. More importantly we can also run hypodd from the command line using the same files generated by `scrtdd` (and possible edit those) and view the console output since the hypodd log file doesn't always reports all the errors. Look at `scrtdd` logs to understand where the working directory is and how to run hypodd manually:

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


Finally, remember to set the Hypodd array limits (compilation time options available via *inc files) accordingly with the size of your problem. If the full catalog relocation doesn't seem to relocate at all, you might have probably hit array limits.

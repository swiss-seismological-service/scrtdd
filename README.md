<pre>
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
 *   Developed by Luca Scarabello, Tobias Diehl                            *
 ***************************************************************************/
</pre>
# Description

The SCRTDD module implements a Real Time Double-Difference Event Location method as described in the paper "Near-Real-Time Double-Difference Event Location Using Long-Term Seismic Archives, with Application to Northern California" by Felix Waldhauser.

This module can also be used to perform (non-real-time) full event catalog Double-Difference relocation, which is covered by another paper: "A Double-Difference Earthquake Location Algorithm: Method and Application to the Northern Hayward Fault, California" by Felix Waldhauser et al.

The actual Double-Difference inversion method is currently performed by the hypoDD software, which has to be installed in the system (currently tested on version 1.3 and 2.1b). We might move away from hypoDD eventually and implement the inversion inside SCRTDD itself. Nevertheless we want to make the transition from hypoDD to SCRTDD software easier and the fact that SCRTDD uses hypoDD internally allows the users to make use of their existing hypoDD experience (configuration, expected behaviour and so on).

## Compile

In order to use this module the sources have to be compiled to an executable. Merge them into the Seiscomp3 sources and compile Seiscomp3 as usual.
<pre>
# merge rtdd-addons and Seiscomp3
git clone https://github.com/SeisComP3/seiscomp3.git sc3-src
cd sc3-src
git submodule add -f https://gitlab.seismo.ethz.ch/lucasca/rtdd-addons.git src/rtdd
</pre>
For compiling Seiscomp3, please refer to https://github.com/SeisComP3/seiscomp3#compiling



# Getting Started

## 1. Install hypoDD software

Go to https://www.ldeo.columbia.edu/~felixw/hypoDD.html and dowload hypoDD. Untar the sources, go to src folder, run make. Done! 


## 2. Define a background catalog for real time relocations

New origins will be relocated in real time against a background catalog of high quality locations. Those high quality events that form the background catalog can be already present in seiscompo database or not. In the latter case the background catalog has to be generated. Either way, the catalog has to be specified in scrtdd when configuring it.

![Catalog selection option](/img/catalog-selection.png?raw=true "Catalog selection")

If we already have a high quality catalog, then we can easily specify it in scrtdd configuration as a path to a file containing a list of origin id or event id (in which case the preferred origin will be used). The file format must be a csv file in which there must be at least one column called seiscompId, from which the ids will be fetched by scrtdd.

E.g. *file myCatalog.csv*

```
seiscompId
event2019ducmfd
Origin/20181214107387.056851.253104
event2019sfamfd
Origin/20180053105627.031726.697885
Origin/20190121103332.075405.6234534
event2019dubnfr
Origin/20190223103327.031726.346363
```

In the more general case in which we don't already have a background catalog, we have to build one using scrtdd. Generating a background catalog involves two steps: first we select the candidate events that might have low quality locations and then we use scrtdd to relocated those to achieve the quality needed for a background catalog. 

Once the candidate events have been selected, we write their seiscomp ids in a file using the same format as the example above (myCatalog.csv). Once that's done, we run the command:

```
scrtdd --dump-catalog myCatalog.csv
```

Or, if the events resides on a different machine we can use the -d option:

```
scrtdd --dump-catalog myCatalog.csv  -d  mysql://user:password@host/seiscompDbName
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
eventId,stationId,isotime,weight,type,networkCode,stationCode,locationCode,channelCode
1,CHNALPS,2014-01-10T04:47:06.78218Z,0.95,P,CH,NALPS,,HHR
1,CHBNALP,2014-01-10T04:47:05.918759Z,0.67,P,CH,BNALP,,HHZ
1,CHFUSIO,2014-01-10T04:47:04.812236Z,0.95,P,CH,FUSIO,,HHR
1,FRRSL,2014-01-10T04:47:02.689842Z,1.06,P,FR,RSL,00,HHZ
1,CHGRIMS,2014-01-10T04:47:01.597023Z,0.95,P,CH,GRIMS,,HHR
1,IVMRGE,2014-01-10T04:46:58.219541Z,0.95,P,IV,MRGE,,HHR
```

Now that we have dumped the events (event.csv, phase.csv, stations.csv) we might perform some editing of those files, if required, then we relocate them. To do so we need to create a new profile inside scrtdd configuration. In this profile we set the generated files (event.csv, phase.csv, stations.csv) as the catalog of the profile. Then we can configure the other profile options that control the relocation process.

![Relocation options](/img/difftraveltime.png?raw=true "Relocation options")
![Relocation options](/img/xcorr.png?raw=true "Relocation options")  

Once we are happy witht he options, we can relocate the catalog with the command:

```
scrtdd --reloc-catalog profileName [--force]
```

scrtdd will relocated the catalog and will generate another set of files event.csv phase.csv and stations.csv (save the previous files somewhere before relocating the catalog). At this point we should check the relocated events and see if we are happy with the results. If not, we change scrtdd settings and relocate the catalog again until we are satisfied with the locations, at which point we can finally set the resulting relocated catalog as background catalog in a scrtdd profile (a new one or the previous one) that we will use for real time relocation.

We are now ready to perform real time relocation!

Note: it is possible to use ph2dt utility to perform the catalog relocation. It this case the scrtdd configuraion for generating dt.ct and dt.cc files will not be used. Instead, ph2dt will be run to generate dt.ct file, and for each entry in the generated dt.ct file the cross correlation will be performed and the relative dt.cc file created. 

```
scrtdd --reloc-catalog profileName --use-ph2dt /some/path/ph2dt.inp [--ph2dt-path /some/path/ph2dt] [--force]
```


## 3. Real time single origin relocation

Real time relocation uses the same configuration we have seen in full catalog relocation, but real time relocation is done in two steps. Each one controlled by a specific hypoDD configuration:

![Relocation options](/img/hypoDDcfg.png?raw=true "Relocation options") 

Step 1: location refinement. In this step hypoDD is used to compute a preliminary relocation of the origin using only catalog absolute travel time entries (dt.ct only).

Step 2: the refined location is used to perform a more precise relocation using both catalog absolute travel times (dt.ct) and differential travel times from cross correlation (dt.cc). 

After step2 the relocated origin is sent to the messaging system. If step2 fails, then the relocated origin from step1 is sent to the messaging system.

If both step1 and step2 fail, then a relocation is reattepted at a later time, accordingly to `delayTimes` option.

Note: when performing the catalog relocation ("scrtdd --reloc-catalog") it is done in a single step

To test the real time relocation we can either run playbacks or use two command line options which relocate existing origins:

```
  -O [ --origin-id ] arg                Reprocess the origin (or multiple comma-separated origins)
                                        and send a message
```

E.g. if we want to process an origin or event, we can run the following command and then check on scolv the relocated origin (the messaging system must be active):

```
scrtdd -O event2019dubnfr [--force]
```

Alternatively we can reprocess an XML file:


```
  --ep arg                              Event parameters XML file for offline 
                                        processing of contained origins (imply 
                                        test option). Ech origin will be 
                                        processed accordingly with the matching
                                        profile configuration
```

E.g.

```
scrtdd --ep event.xml [--force]
```


## 4. Troubleshooting

Check log file: ~/.seiscomp/log/scrtdd.log 

Alternatively, when running scrtdd from the command line use the following options:

```
# set log level to debug and log to the console (standard output) insted of log file
--verbosity=4 --console=1
```

A useful option we can find in scrtdd configuration is `keepWorkingFiles`, which prevent the deletion of scrtdd processing files from the working directory. In this way we can access the working folder and check input, output files used for running hypodd. Make sure to check the `*.out` files, which contain the console output of hypodd (sometimes we can find errors only in there, as they do not appear in hypodd.log file).

Another useuful command line option is --dump-wf, which allows to dump the waveforms used for cross correlation after the filtering and resampling have been applied.

Finally, remember to set the Hypodd array limits (compilation time options available via *inc files) accordingly with the size of your problem. If the full catalog relocation doesn't seem to relocate at all, you might have probably hit array limits.

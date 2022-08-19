.. _single-event-label:

Real-time single-event relocation
=================================

Single-event relocation is used to relocate events in real-time and it requires a background catalog to work.

.. figure:: media/singleEventRelocationSyntDataExample.png
   :width: 800
   
   Test with synthetic data from the unit testing folder. Events from 4 clusters have their locations and times altered - using several normal distributions with non-zero mean - to simulate location/time errors. The single-event double-difference inversion is then applied on those altered events, one at a time against the background catalog and their original locations and times are properly recovered. It is interesting to note that the backgroud catalog is not part of the 4 clusters, but the double-differene inversion is still able to perfectly recover the event locations of those close-by clusters.

------
Summay
------

* Use the multi-event relocation feature to prepare a background catalog
* Create a rtDD profile or use the same profile used for generating the background catalog, then set the profile background catalog and add the profile to the list of active real-time profiles (``activeProfiles`` parameter). The default profile parameter values are meant to be a good starting choice, so there is no need to tweak them heavily. However, it is a good choice to configure a custom velocity model (``solver.travelTimeTable``)
* Make sure to read :ref:`avoid-loop-label` paragraph to avoid a potential issue
* Enable and start rtDD (``seiscomp enable scrtdd``, ``seiscomp start scrtdd``)

--------------
The long story
--------------

To enable the real-time processing a profile should be created and enabled by including it in ``scrtdd.activeProfiles`` option.
 
In real-time processing rtDD relocates new origins, one a time as they occur, against a background catalog of high quality events. Those high quality events can be generated via multi-event relocation, which has already been covered in the previous sections.

Real time relocation uses the same configuration we have seen in full catalog relocation, but real time relocation is done in two steps:

**Step 1**: location refinement. In this step rtDD performs a preliminary relocation of the origin where the differential travel times in the double-difference system are derived from the pick times.

**Step 2**: the refined location computed in the previous step is used as starting location to perform a more precise relocation using cross-correlation to refine the differential travel times. If step1 fails, step2 is attempted anyway.

If step2 completes successfully the relocated origin is sent to the messaging system. 

--------------------------------
Configuring a background catalog
--------------------------------

The easiest choice is to use as background catalog the relocated multi-event results; the triplet *reloc-event.csv*, *phase.csv*, *station.csv*:

.. image:: media/catalog-selection3.png
   :width: 800

However, if the catalog is generated in XML format, it can be imported in the SeisComP database. In this case the background catalog can be a file containing just the origin ids. 

.. image:: media/catalog-selection1.png
   :width: 800

While it is neat to have the background catalog in the SeisComP database, this approach has few limitations. First it may take a lot of time for rtDD to load a big catalog from the database comparing to loading it from files. Also, since the background catalog should be periodically updated, old events are continuously updated with new origins, which can lead to a not optimal database performance-wise.

Once the background catalog is configured rtDD can be enabled and started as any other SeisComP module.  New origins will be relocated as soon as they arrive in the messaging system.

-------
Testing
-------

You might consider testing the configuration relocating some existing events to make sure the parameters are suitable for your use case. To test the real time relocation there are two command line options which relocate existing origins::

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


Relocate origin ID and send the relocation to the messaging system for further processing
-----------------------------------------------------------------------------------------

If we want to process an origin we can run the following command and then check on ``scolv`` the relocated origin (the messaging system must be active). This is mostly useful when we want to relocate an origin on a running system and keep the relocation::

    scrtdd --origin-id someOriginId \
           --verbosity=3 --console=1 [db options] 


Relocate origin ID but do not send the relocation (debug)
---------------------------------------------------------

As above but add ``--test`` and the origin will not be sent to the messaging system. Useful for troubleshooting when the ``scrtdd.saveProcessingFiles`` option is enabled to verify the relocation files in ``scrtdd.workingDirectory``.
::

    scrtdd --origin-id someOriginId --test \
           --verbosity=3 --console=1 [db options]

Relocate origin ID and store the result to XML file
---------------------------------------------------

Adding the ``--xmlout`` option allows to save the origin as a XML file. We can finally open the ile with ``scolv`` for inspection::

    scrtdd --origin-id someOriginId --xmlout \
           --verbosity=3 --console=1 [db options] \
      >  relocated-origin.xml

Relocate XML file and store the result to XML file
--------------------------------------------------

Similarly to other SeisComP commands the ``--ep`` option can be used for full offline processing. All origins contained in the input XML file are relocated::

    scrtdd --ep origin.xml --verbosity=3 --console=1 [db options] \
      > relocated-origin.xml

Relocation log
--------------

Here we report an example *single-event* relocation log::

    [info] Starting DD relocator in single event mode: event 1 lat 46.419079 lon 7.942911 depth 8.9902 mag 0.56 time 2020-10-31T19:46:57.703383Z #phases 22
    [info] Performing step 1: initial location refinement (no cross-correlation)
    [info] Found 22 neighbouring events
    [info] Building and solving double-difference system...
    [...]
         ...details of the solutions for each iteration of the solver
    [...]
    [info] Successfully relocated 1 events, RMS median 0.2865 [sec] median absolute deviation 0.0000 [sec]
    [info] Events RMS before relocation: median 0.3309 median absolute deviation 0.0000
    [info] Step 1 relocation successful, new location: lat 46.419460 lon 7.932872 depth 8.9892 time 2020-10-31T19:46:57.770484Z
    [info] Relocation report: 
           Origin changes: location=0.77[km] depth=-0.00[km] time=0.067[sec] 
           Rms change [sec]: -0.044 (before/after 0.331/0.287) 
           Neighbours=22 
           Used Phases: P=9 S=6 
           Stations distance [km]: min=16.6 median=25.6 max=61.9 
           DD observations: 143 (CC P/S 0/0 TT P/S 88/55) 
           DD residuals [msec]: before=40+/-59.4 after=-4+/-4.9
    
    [info] Performing step 2: relocation with cross-correlation
    [info] Found 30 neighbouring events
    [info] Computing differential times via cross-correlation...
    [info] Cross-correlation performed 101 (P phase 50%, S phase 50%), skipped 89 (47%)
    [info] Cross-correlation success (coefficient above threshold) 73% (74/101). Successful P 86% (44/51). Successful S 60% (30/50)
    [info] Building and solving double-difference system...
    [...]
         ...details of the solutions for each iteration of the solver
    [...]    
    [info] Successfully relocated 1 events, RMS median 0.2834 [sec] median absolute deviation 0.0000 [sec]
    [info] Events RMS before relocation: median 0.2642 median absolute deviation 0.0000
    [info] Step 2 relocation successful, new location: lat 46.418945 lon 7.932328 depth 8.6810 time 2020-10-31T19:46:57.808104Z
    [info] Relocation report:
           Origin changes: location=0.07[km] depth=-0.31[km] time=0.038[sec] 
           Rms change [sec]: 0.019 (before/after 0.264/0.283) 
           Neighbours=30 
           Used Phases: P=9 S=6 
           Stations distance [km]: min=16.4 median=25.4 max=61.7 
           DD observations: 190 (CC P/S 44/30 TT P/S 72/44) 
           DD residuals [msec]: before=40+/-59.4 after=-5+/-6.5


rtDD adds two comments to each relocated origin: ``relocation::sourceOrigin`` and ``relocation::report``. 

``relocation::sourceOrigin`` contains the id of the origin that triggered the relocation. ``relocation::report`` contains a summary of the relocation process. E.g.::

    Origin changes: location=0.23[km] depth=1.40[km] time=-0.147[sec]
    Rms change [sec]: -0.153 (before/after 0.502/0.349)
    Neighbours=80 Used Phases: P=37 S=16
    Stations distance [km]: min=15.9 median=57.0 max=99.8
    DD observations: 687 (CC P/S 141/47 TT P/S 375/124)
    DD residuals [msec]: before=-106+/-21.6 after=9+/-26.2

They can be both visualized in ``scolv`` as additional columns adding the following settings to ``scolv.cfg``::

    # SCRTDD: display source origin that generated a scrtdd relocation
    eventlist.customColumn.default = ""
    eventlist.customColumn.originCommentID = relocation::sourceOrigin
    eventlist.customColumn = triggeringOrigin
    
    # SCRTDD: display origin comment containing rtdd relocation report
    eventedit.customColumn.default = ""
    eventedit.customColumn.originCommentID = relocation::report
    eventedit.customColumn.pos = 99
    eventedit.customColumn = scrtd


.. _phase-update-label:

------------
Phase update
------------

rtDD uses cross-correlation to fix the pick time and uncertainty of automatic picks. The pick time is updated according to the average lag detected by all the good (above configured threshold) cross-correlation results. Since the real-time events are cross-correlated against catalog events, which have good manual picks, the updated pick time is expected to improve. The pick uncertainty is derived from the uncertainties of catalog-events. If no cross-correlation result is above the configured threshold, the automatic pick is kept untouched.
 
rtDD can also use cross-correlation to detect phases at stations with no associated picks (see ``crossCorrelation.detectMissingPhasesAutoOrigin`` and ``crossCorrelation.detectMissingPhasesManualOrigin``). It firstly computes the theoretical pick time and then cross-correlates it against the catalog event phases. Every theoretical pick that has at least one good cross-correlation result is added to the relocated origin, with pick time and uncertainties derived from catalog phases (similarly to what is done for automatic picks). Those detected picks are thus used in the double-difference system inversion. Theoretical picks that have no good cross-correlation results are simply discarded.

Picks that have been updated or created by rtDD are identifiable by a ``x`` suffix (Px, Sx).

Manual picks are never modified.

.. _avoid-loop-label:

-------------------------
Avoiding Relocation Loops
-------------------------

rtDD listens and sends messages to the LOCATION group. In a default installation where the only locator is ``scautoloc`` that's not an issue: ``scautoloc`` will send an origin to LOCATION and rtDD will receive it and send an updated origin to LOCATION.  However, when there are multiple (re)locators (e.g. scanloc, screloc) that listen to LOCATION and send their own updated origin to LOCATION too, then an infinite loop happens! In this case a new messaging group needs to be created, e.g. RELOCATION, so that the origins flow from LOCATION to RELOCATION without going back.

 E.g. of a properly configured system::


                          LISTEN                       SEND 
                  (MessagingSubscription)      (PrimaryMessagingGroup)
    scautoloc             ...                        LOCATION
    scanloc       LOCATION, ...                      LOCATION
    screloc       LOCATION, ...                     RELOCATION
    scrtdd        LOCATION, ...                     RELOCATION
    scevent       LOCATION,RELOCATION, ...             ...
    scamp         LOCATION,RELOCATION, ...             ...
    scmag         LOCATION,RELOCATION, ...             ...


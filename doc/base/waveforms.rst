.. _waveform-label:

Waveforms data
==============

SeisComP applications access waveform data through the RecordStream interface (see `official SeisComP documentation <https://www.seiscomp.de/doc/base/concepts/recordstream.html>`_ for more details) and it can be configured in *global.cfg* or passed via command line with ``-I URI``. The RecordStream parameters define the service(s) through which SeisComP can access real-time and historical waveforms data (seedlink, fdsn, sds archive, `etc <https://www.seiscomp.de/doc/apps/global_recordstream.html>`_). 

RecordStream configuration for Single-Event (real-time)
-------------------------------------------------------

A hypothetical RecordStream configuration might look like this::

    recordstream = combined://slink/localhost:18000;sdsarchive//path/to/miniseed

This configuration is a combination of seedlink and sds archive, which allows rtDD to retrieve catalog waveforms via sds and real-time event data via seedlink.

Please note that depending on the responsiveness of the seedlink server the real-time relocations may incur in delays. A couple of configuration options allow to control those delays: *timeout* and *retries*. The example below forces a timeout of 5 seconds (default is 5 minutes) and to not reconnect. In case of a timeout, rtDD will proceed with the available data, without further delays::

    recordstream = combined://slink/localhost:18000?timeout=5&retries=0;sdsarchive//path/to/miniseed
 
RecordStream configuration for Multi-Event
------------------------------------------

To access the catalog waveforms a RecorStream that connects to a historical archive is required. Common formats are fdsn and sds archive. However it is also common for users to have multiple mseed files not arranged in any specific way.

Example of accessing waveforms from a FDSN service::

    scrtdd -I fdsnws://service.iris.edu:80/fdsnws/dataselect/1/query [...options...]

Example of accessing waveforms from a sds service::

    scrtdd -I  sdsarchive:///path/to/archive [...options...]

If the waveforms are store in multiple miniseed files, those could be concatenated in a single file and passed to scrtdd like the following, but it might be slow depending on the file size::

    cat file1.mseed file2.mseed ... fileX.mseed > data.mseed
    
    # sort records and remove duplicates (https://www.seiscomp.de/doc/apps/scmssort.html)
    scmssort -u -E -v data.mseed > sorted.mseed 
    
    scrtdd -I file://sorted.mseed [...options...]

A better approach would be to convert the sorted miniseed file to a sds archive and use that archive as waveform data source::

    # creates the sds archive (https://www.seiscomp.de/doc/apps/scart.html)
    mkdir my-sdsarchive
    scart  -I file://./sorted.mseed ./my-sdsarchive
    
    # use the sds archive
    scrtdd -I  sdsarchive://./my-sdsarchive [...options...]

Performance wise it might be unfeasible to concatenate all files (``cat``) and sort the records (``scmssort``) when the amout of data is very large. However, since the SDS archive splits the data by station in daily files it is possible to take advantage of that and perform the operations in subsets of files. That it we can concatenate files and sort the result on a day by day basis (assuming each file contain data within a day and does not overlap multiple days). If each file contains data for a single station only, we can work on even smaller subset of files by working on a single station files only. For each subset of files we can finally run ``scart`` command safely on the sorted file generated.

Waveforms data caching
----------------------

Unless the recordStream points to a local disk storage, downloading waveforms might require a lot of time. For this reason rtDD stores the waveforms to disk (called waveform cache) after downloading them. This applies only to the catalog event waveforms, which are used over and over again. That's not true for the real-time events, whose waveforms are used just once and never cached. The cache folder is ``workingDirectory/profileName/wfcache/`` and the ``workingDirectory`` is configurable.

However, for certain situations (e.g. debugging) it might be useful to cache all the waveforms, even the ones that are normally not cached. For those special cases the option --cache-wf-all can be used (stored in ``workingDirectory/profileName/tmpcache/`` which can be deleted afterwards).


Catalog waveforms preloading
----------------------------

When rtDD starts for real-time processing (``seiscomp start scrtdd``) it loads all the catalog waveforms and stores them to disk (if they are not already there) if the option ``performance.profileTimeAlive`` is 0. Otherwise the catalog waveforms will be loaded only when needed (lazy loading), that is during cross-correlations of real-time events against the background catalog. We can however force rtDD to pre-download all waveforms before starting the module using the following option::

    scrtdd --load-profile-wf --profile myprofile [-I RecordStream]

::

    scrtdd --help
      --load-profile-wf                     Load catalog waveforms from the 
                                            configured recordstream and save them 
                                            into the profile working directory. Use
                                            in combination with --profile

      --profile arg                         To be used in combination with other 
                                            options: select the profile 
                                            configuration to use



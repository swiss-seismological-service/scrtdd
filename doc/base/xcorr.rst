.. _xcorr-event-label:

Cross-Correlation
=================

Cross-correlation is used to refine the differential travel time between two phase picks at a common station for an event pair. If the correlation coefficient exceeds a configured threshold, the cross-correlation lag is used to correct the absolute travel-time difference as the observed differential time in the double-difference system (see :ref:`multi-event-relocation-process-label`). Phase pairs with correlation coefficients below the threshold retain the differential time derived from their original pick times.

Because cross-correlation is computationally intensive, it is recommended to enable it only after other relocation parameters have been optimized.

Please note that certain cross-correlation settings require inventory information, for example if the horizontal components have to be used, the inventory is required to know which are the horizontal ones.

.. _waveform-label:

-------------
Waveform Data
-------------

SeisComP applications access waveform data through the RecordStream interface (see the `official SeisComP documentation <https://www.seiscomp.de/doc/base/concepts/recordstream.html>`_ for more details). This can be configured in ``global.cfg`` or specified via the command line using the ``-I URI`` option. RecordStream parameters define the services through which SeisComP accesses real-time and historical waveform data (e.g., SeedLink, FDSN, SDS archives, etc.).

RecordStream Configuration for Single-Event (Real-Time)
-------------------------------------------------------

A typical RecordStream configuration might look like this::

    recordstream = combined://slink/localhost:18000;sdsarchive//path/to/miniseed

This example combines SeedLink and an SDS archive, allowing ``scrtdd`` to retrieve catalog waveforms from the archive and real-time event data via SeedLink.

Depending on the responsiveness of the SeedLink server, real-time relocations may incur delays. If data is unavailable, SeisComP may attempt to reconnect indefinitely. To prevent excessive delays, use the ``timeout`` and ``retries`` parameters. The following example sets a timeout of 5 seconds (the default is 5 minutes) and disables reconnections for missing data, allowing ``scrtdd`` to proceed with available data::

    recordstream = combined://slink/localhost:18000?timeout=5&retries=0;sdsarchive//path/to/miniseed


RecordStream Configuration for Multi-Event
------------------------------------------

Accessing catalog waveforms requires a RecordStream connected to a historical archive, such as an FDSN service or an SDS archive.

Example accessing waveforms from an FDSN service::

    scrtdd -I fdsnws://service.iris.edu:80/fdsnws/dataselect/1/query [...options...]

Example accessing waveforms from an SDS archive::

    scrtdd -I sdsarchive:///path/to/archive [...options...]


Waveform Data Caching
---------------------

As downloading waveforms can be time-consuming, ``scrtdd`` can cache reference catalog waveforms to disk once they have been fetched. Real-time event waveforms are used immediately and are not cached. Cache behavior and storage locations can be controlled via configuration options.

Waveforms are loaded "lazily," only when required for cross-correlation (e.g., when a phase pair is processed). You can manually force ``scrtdd`` to pre-download all reference catalog waveforms for a profile using the following command::

    scrtdd --load-profile-wf --profile myprofile  \
           --verbosity=3 --console=1  [-I RecordStream] [db options]


.. _reusing-xcorr-label:

---------------------------------
Reusing Cross-Correlation Results
---------------------------------

To save processing time, cross-correlation results can be reused. This is particularly valuable when experimenting with relocation settings while keeping cross-correlation parameters constant. When the ``--dump-diagnostics`` option is used, the relocation process generates a ``xcorr.csv`` file. This file contains computed cross-correlation lags and can be passed back to ``scrtdd`` using the ``--xcorr-cache`` option.

When this option is active, ``scrtdd`` skips cross-correlation for any event pair already present in the cache file. A common use case is evaluating relocation quality across different correlation coefficient thresholds without recomputing waveform cross-correlations.

-------------------
Waveform Inspection
-------------------

The ``--dump-wf`` option instructs ``scrtdd`` to dump the processed waveforms of the specified catalog to disk. To avoid dumping the entire catalog, you can edit the input ``event.csv`` file to include only a small subset of events. Waveforms are saved in miniSEED format and can be inspected using external tools like ``scrttv`` or ObsPy. These files reflect the waveforms after filtering and resampling have been applied.

Example::

    scrtdd --dump-wf station.csv,event.csv,phase.csv --profile myProfile \
           --verbosity=3 --console=1  [-I RecordStream] [db options]


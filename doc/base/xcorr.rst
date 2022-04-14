.. _xcorr-event-label:

Cross-correlation
=================

Good cross-correlation results are needed to achieve high quality double-difference observations, which in turn results in high resolution relocations. The purpose of the cross-correlation is to find the exact time difference between two picks of an event pair at a common station. The cross-correlation is automatically performed by rtDD before the double-difference inversion when ``RecordStream`` is configured, otherwise it is simply skipped, The cross-correlation step can also be disabled setting the configuration parametters ``crossCorrelation.maxStationDistance`` and/or ``crossCorrelation.maxInterEventDistance`` to 0.

------------------
Eval-xcorr command
------------------

The ``--eval-xcorr`` option should be used to properly configure the cross-correlation parameters. The optimization process involves running ``--eval-xcorr`` with different configuration and analyzes the results. The goal is to have as many matches as possible (increase correlation coefficient) avoiding bad/false matches (``lag`` values higher than the expected pick time uncertainty are probably an indication of false matches): this is a trade-off.

Example of ``--eval-xcorr`` command::

    scrtdd --eval-xcorr station.csv,event.csv,phase.csv --profile myProfile --verbosity=3 --console=1

Example output::

    [...]
    ---FINAL STATS--
    Total xcorr P phases   #CC   #Skip meanCC (MAD) medianCC (MAD) meanLag (MAD) medianLag (MAD)
                         41877   61126  0.77 (0.10)    0.77 (0.10)      -0 ( 36)        -0 ( 15)
    Total xcorr S phases   #CC   #Skip meanCC (MAD) medianCC (MAD) meanLag (MAD) medianLag (MAD)
                         30730   66292  0.85 (0.06)    0.85 (0.06)      -1 ( 42)        -0 ( 28)
    Xcorr P phases by inter-event distance in 0.10 km step
     EvDist [km]           #CC   #Skip meanCC (MAD) medianCC (MAD) meanLag (MAD) medianLag (MAD)
     0.00-0.10            4764    4508  0.82 (0.10)    0.84 (0.09)       0 ( 28)        -0 ( 11)
     0.10-0.20           13906   18615  0.78 (0.10)    0.78 (0.10)       0 ( 32)         0 ( 14)
     0.20-0.30           13633   21365  0.76 (0.10)    0.76 (0.09)      -1 ( 41)         0 ( 17)
     0.30-0.40            5445    9591  0.75 (0.09)    0.75 (0.09)      -5 ( 44)        -1 ( 20)
     0.40-0.50            2360    3807  0.75 (0.10)    0.75 (0.09)       5 ( 40)         1 ( 18)
     0.50-0.60            1010    1727  0.75 (0.10)    0.75 (0.10)       9 ( 35)         3 ( 14)
     0.60-0.70             553    1114  0.75 (0.10)    0.73 (0.10)       3 ( 37)         1 ( 14)
     0.70-0.80             142     247  0.74 (0.11)    0.70 (0.09)      20 ( 37)         4 ( 10)
     0.80-0.90              60     142  0.76 (0.12)    0.77 (0.13)      14 ( 36)         2 ( 10)
     0.90-1.00               4      10  0.79 (0.11)    0.80 (0.10)      25 ( 28)         7 (  1)
    Xcorr S phases by inter-event distance in 0.10 km step
     EvDist [km]           #CC   #Skip meanCC (MAD) medianCC (MAD) meanLag (MAD) medianLag (MAD)
     0.00-0.10            4024    4472  0.89 (0.06)    0.91 (0.05)      -0 ( 37)         0 ( 25)
     0.10-0.20           10521   20221  0.85 (0.06)    0.86 (0.06)       0 ( 40)        -0 ( 28)
     0.20-0.30            9961   23181  0.84 (0.06)    0.84 (0.05)      -0 ( 45)        -0 ( 31)
     0.30-0.40            3802   10307  0.83 (0.06)    0.83 (0.05)      -7 ( 47)        -4 ( 30)
     0.40-0.50            1430    4307  0.83 (0.06)    0.83 (0.05)       0 ( 43)        -0 ( 28)
     0.50-0.60             572    2002  0.83 (0.05)    0.83 (0.05)      16 ( 40)        10 ( 30)
     0.60-0.70             327    1357  0.83 (0.05)    0.83 (0.05)      -7 ( 42)        -7 ( 30)
     0.70-0.80              64     287  0.82 (0.05)    0.84 (0.04)       4 ( 33)         2 ( 24)
     0.80-0.90              26     150  0.83 (0.05)    0.84 (0.05)      13 ( 45)        -3 ( 26)
     0.90-1.00               3       8  0.80 (0.05)    0.83 (0.03)       9 (  8)        14 (  2)
    XCorr P phases by event to station distance in 3.00 km step
    StaDist [km]           #CC   #Skip meanCC (MAD) medianCC (MAD) meanLag (MAD) medianLag (MAD)
      3-6                   22      17  0.64 (0.04)    0.64 (0.04)      19 ( 23)         6 (  8)
      6-9                13487   12847  0.77 (0.11)    0.77 (0.11)       0 ( 26)         0 (  9)
      9-12                7942    9451  0.79 (0.11)    0.81 (0.11)      -1 ( 21)         0 (  9)
     12-15                7238   11380  0.73 (0.09)    0.72 (0.08)      -0 ( 31)         0 ( 18)
     15-18                 994     679  0.76 (0.07)    0.76 (0.07)      -5 ( 64)        -2 ( 36)
     18-21                3067    4218  0.78 (0.08)    0.79 (0.07)      -1 ( 35)        -0 ( 21)
     21-24                1318    5599  0.73 (0.10)    0.73 (0.09)       2 ( 90)         2 ( 66)
    [...]
    XCorr S phases by event to station distance in 3.00 km step
    StaDist [km]           #CC   #Skip meanCC (MAD) medianCC (MAD) meanLag (MAD) medianLag (MAD)
      3-6                   35       5  0.89 (0.04)    0.89 (0.04)      -5 ( 16)        -3 ( 14)
      6-9                 8307   14072  0.87 (0.06)    0.88 (0.05)       0 ( 28)         0 ( 22)
      9-12                7740   13058  0.84 (0.06)    0.84 (0.05)      -0 ( 44)        -0 ( 33)
     12-15                3384   10766  0.83 (0.06)    0.83 (0.05)      -1 ( 36)         0 ( 21)
     15-18                1207    1849  0.83 (0.06)    0.83 (0.05)       1 ( 52)         1 ( 37)
     18-21                 527    2601  0.82 (0.07)    0.81 (0.06)      -0 ( 54)         2 ( 39)
     21-24                2411    7921  0.82 (0.06)    0.81 (0.05)      -1 ( 48)        -2 ( 36)
    [...]
    XCorr P phases by station
    Station                #CC   #Skip meanCC (MAD) medianCC (MAD) meanLag (MAD) medianLag (MAD)
    8D.RAW1.              4254    5109  0.75 (0.09)    0.74 (0.09)      -0 ( 35)         1 ( 26)
    8D.RAW2.              6436    4086  0.75 (0.10)    0.74 (0.09)       0 ( 22)         0 ( 10)
    8D.RAW4.              2902    3045  0.79 (0.08)    0.79 (0.07)      -0 ( 30)         0 ( 20)
    C4.CERNS.                0       0  0.00 (0.00)    0.00 (0.00)       0 (  0)         0 (  0)
    CH.AIGLE.               41     211  0.80 (0.08)    0.84 (0.06)      -6 ( 78)         1 ( 31)
    CH.DIX.                687    5012  0.75 (0.09)    0.75 (0.08)       1 ( 48)         0 ( 21)
    [...]
    XCorr S phases by station
    Station                #CC   #Skip meanCC (MAD) medianCC (MAD) meanLag (MAD) medianLag (MAD)
    8D.RAW1.              2539    7670  0.83 (0.06)    0.83 (0.05)      -0 ( 25)        -0 ( 16)
    8D.RAW2.              6861    3732  0.88 (0.06)    0.89 (0.05)      -0 ( 28)         0 ( 23)
    8D.RAW4.                79     838  0.86 (0.06)    0.87 (0.06)      -5 ( 32)        -5 ( 27)
    CH.AIGLE.              113     288  0.82 (0.04)    0.82 (0.04)      -1 ( 85)        -5 ( 82)
    CH.DIX.               2394    2606  0.84 (0.06)    0.85 (0.05)       1 ( 35)         0 ( 22)
    [...]


* ``#CC``: number of cross-correlations performed
* ``#Skip``: number of cross-correlations whose results do not account for the computation of the statistics
* ``coeff``: correlation coefficient between phase waveforms 
* ``lag``: cross-correlation lag between phase waveforms in milliseconds

There could be several reasons why the cross-correlation between 2 phase waveforms is not considered for computing the statistics: the correlation coefficient is below the configured threshold (see ``crossCorrelation.x-phase.minCCCoef``), the SNR of one or both the waveforms is below the configured threshold (see ``crossCorrelation.snr.minSnr``), the waveform data for one or both the phases is not available and in general when the it is not possible to perform the cross-correlation. It is possible to know the exact reason by looking at the logs at debug level (--verbosity=4).

The SNR is particularly important to reject bad automatic picks or picks detected via cross-correlation (see :ref:`phase-update-label`), but but it is not so relevant when relocating manually reviewed origins since the picks are checked already and bad ones discarded.


.. _reusing-xcorr-label:

---------------------------------
Reusing cross-correlation results
---------------------------------

When cross-correlation settings are not changed, it might be useful to reuse the cross-correlation results to save processing time. Both the ``--eval-xcorr`` and ``--reloc-catalog`` options save a ``xcorr.csv`` file after finishing their execution (thay overwrite it if already present!). That file contains the computed cross-correlation results and can be given back to rtDD via the command line option ``--xcorr-cache``. It is safe to change the value of ``crossCorrelation.x-phase.minCCCoef`` and reuse the cross-correlation results to see how performance change at varying correlation-coefficient threshold.

--------------------
Waveforms inspection
--------------------

The ``--dump-wf`` option will make rtDD dump to disk the waveforms of the catalog passed as argument. Those files are in miniseed format and can be viewed with an external tool (e.g. ``scrttv waveform.mseed``) or obspy). The waveforms are written to disk after the filterting and resampling have been applied::

    scrtdd --help
      --dump-wf arg                         Dump processed waveforms of the catalog
                                            passed as argument in the current 
                                            working directory.The catalog can be a 
                                            single file (containing seiscomp origin
                                            ids) or a file triplet 
                                            (station.csv,event.csv,phase.csv). Use 
                                            in combination with --profile.


e.g.::

    scrtdd --dump-wf station.csv,event.csv,phase.csv --profile myProfile --verbosity=3 --console=1
    
    17:59:28 [info] Writing ev1.8D.RAW2..HHT.Sg.manual.mseed
    17:59:28 [info] Writing ev1.CH.SAYF2..HGT.Sg.manual.mseed
    17:59:28 [info] Writing ev1.CH.SENIN..HHT.Sg.manual.mseed
    17:59:28 [info] Writing ev1.XY.LEO01..HHT.Sg.manual.mseed
    17:59:28 [info] Writing ev1.XY.LEO01..HHZ.Sg.manual.mseed
    17:59:28 [info] Writing ev1.FR.OGSI.00.HHZ.Pg.manual.mseed
    17:59:28 [info] Writing ev1.GU.REMY..HHZ.Pg.manual.mseed
    17:59:28 [info] Writing ev1.CH.FIESA..HHZ.Pg.manual.mseed
    17:59:28 [info] Writing ev1.CH.TORNY..HHZ.Pg.manual.mseed
    17:59:28 [info] Writing ev1.8D.AMIDI..EHZ.Pg.manual.mseed
    17:59:28 [info] Writing ev2.CH.DIX..HHT.Sg.manual.mseed
    17:59:28 [info] Writing ev2.8D.RAW2..HHZ.Pg.manual.mseed
    17:59:28 [info] Writing ev2.CH.SAYF2..HGZ.Pg.manual.mseed
    17:59:28 [info] Writing ev2.CH.STSW2..HGZ.Pg.manual.mseed
    [...]


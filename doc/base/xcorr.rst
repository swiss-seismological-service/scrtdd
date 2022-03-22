.. _xcorr-event-label:

Cross-correlation
=================

Good cross-correlation results are needed to achieve high quality double-difference observations, which in turn results in high resolution relocations. The purpose of the cross-correlation is to find the exact time difference between two picks of an event pair at a common station. The cross-correlation is automatically performed by rtDD before the double-difference inversion when `RecordStream` is configured, otherwise it is simply skipped, The cross-correlation step can also be disabled setting the configuration parametters `crossCorrelation.maxStationDistance` and/or `crossCorrelation.maxInterEventDistance` to 0.

----
Logs
----

Some cross-correlation statistics are printed in both multi-event and single-event mode. Those can be seen in the log file or in the console output (with options `--console=1 --verbosity=3`)::

    [info] Cross-correlation statistics: performed 40361, waveforms with Signal to Noise ratio too low 2435, waveforms not available 98
    [info] Total xcorr 40361 (P 59%, S 41%) success 28% (11499/40361). Successful P 22% (5300/23844). Successful S 38% (6199/16517)
    [info] xcorr on actual picks 24784/40361 (P 60%, S 40%) success 37% (9186/24784). Successful P 31% (4629/14761). Successful S 45% (4557/10023)
    [info] xcorr on theoretical picks 15577/40361 (P 58%, S 42%) success 15% (2313/15577). Successful P 7% (671/9083). Successful S 25% (1642/6494)

There could be several reasons why the cross-correlation between 2 phase waveforms is skipped: the waveform data for one or both the phases is not available, the configured components (`crossCorrelation.x-phase.components`) were not found for the phase,the SNR of one or both the waveforms is below the configured threshold (see `crossCorrelation.snr.minSnr`, the phases were detected on different channel codes (see `crossCorrelation.compatibleChannels` configuration option), the waveforms of the two phases use different frequencies and the option `crossCorrelation.waveformFilteringiresampling` is not used. It is possible to know the reason on why a cross-correlation was skipped for a particular phase pair looking at the logs at debug level (--verbosity=4).

The statistics are broken down in actual picks and theoretical picks. This is because rtDD computes theoretical picks that are cross-correlated together with detected picks. This is useful to increase the number of double-difference observations. See the [Phase update](#23-phase-update) paragraph for further details.

------------------
Eval-xcorr command
------------------

The `--eval-xcorr` command can be used to evaluate the cross-correlation parameter. 

It is especially interesting to compare the results before/after a relocation since the statistics on cross-correlation are an indirect measure of the proximity of events: we should see higher coeffient values for events close to each other and gradually worsen with increasing inter-event distance. That can be used as a verification of the quality of the relocation::

    scrtdd --eval-xcorr station.csv,event.csv,phase.csv --profile myProfile --verbosity=3 --console=1

Example output::

    [...]
    13:13:17 [info] <FINAL STATS>
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


* `#pha`: how many phases have been cross-correlated
* `pha good CC`: how many of those were successful (correlation coefficient above the configured threshold)
* `coeff`: the average correlation coefficient
* `goodCC/ph`: the average number of good matches per phase (each event phase appears in multiple double-difference observations, so multiple event pais, hence multiple cross-correlations)
* `time-diff`: the average pick time difference detected by the cross-correlation
* `+/-`: whenever sensible, it is also indicated the Mean Absolute Deviation of the value

The `--eval-xcorr` option should be used to properly configure the cross-correlation parameters. The optimization process involves running `--eval-xcorr` with different configuration and analyzes the results. The goal is to have as many matches as possible (increase `GoodCC`) avoiding bad/false matches (very high values of `time-diff` are probably an indication of false matches): this is a trade-off.

The SNR is particularly important to reject bad picks (automatic picks or picks detected via cross-correlation by rtDD). The SNR signal/noise windows should be chosen so that they satisfies ALL the following 5 conditions:

* pick time too early -> we want low SNR
* pick time too late -> we want low SNR
* pick time perfect -> we want high SNR
* pick time is early but acceptable -> we want high SNR
* pick time is late but acceptable -> we want high SNR

--------------------
Waveforms inspection
--------------------

The `--dump-wf` option will make rtDD dump to disk the waveforms of the catalog passed as argument. Those files are in miniseed format and can be viewed with an external tool (e.g. `scrttv waveform.mseed`) or obspy). The waveforms are written to disk after the filterting and resampling have been applied::

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


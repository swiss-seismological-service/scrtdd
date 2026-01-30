.. _intro-label:

Introduction
============

rtDD has two modes of operation: :ref:`multi-event <multi-event-label>` and :ref:`single-event <single-event-label>`.

The **Multi-Event mode** relocates catalogs of events using the double-difference method. This is an offline mode that can be used to generate high-resolution catalogs for external analysis or to :ref:`periodically generate <continuous-label>` a catalog snapshot. Such snapshots serve as reference catalogs for the Single-Event mode. 

The **Single-Event mode** relocates events one by one as they occur in real-time, using a reference catalog.

The methods in rtDD are based on these papers:

* "Near-Real-Time Double-Difference Event Location Using Long-Term Seismic Archives, with Application to Northern California" by Felix Waldhauser
* "A Double-Difference Earthquake Location Algorithm: Method and Application to the Northern Hayward Fault, California" by Waldhauser & Ellsworth

The original methodology from the papers has been extended to allow event relocation considering both :ref:`their absolute and relative locations <inclusion-tt-residual-label>`.

For travel time calculations, rtDD supports the NonLinLoc grid file format (by Anthony Lomax) in addition to the native SeisComP travel time formats (e.g., LOCSAT and libtau). NLL grids are available as a general SeisComP Travel Time plugin for use by any module. See :ref:`ttt-label` for more details.

The double-difference equation system is solved using the LSQR (by Chris Paige and Michael Saunders) and LSMR (by David Fong and Michael Saunders) algorithms.


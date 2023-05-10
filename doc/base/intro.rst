.. _intro-label:

Introduction
============

rtDD has two mode of operations: :ref:`multi-event <multi-event-label>` mode and :ref:`single-event <single-event-label>` mode. The Multi-Event mode can relocate catalogs of events using the double-difference method. That is not a real-time mode and it can be used to generate data for external processing, or to :ref:`periodically generate <continuous-label>` a catalog snapshot, which can become the reference catalog of the single event-mode, where the events are relocated one at a time as they happen in real-time. The Multi-Event mode includes an interesting option where the events are relocated considering both :ref:`their absolute and relative locations <inclusion-tt-residual-label>`.

The methods developed in rtDD are based on the paper "Near-Real-Time Double-Difference Event Location Using Long-Term Seismic Archives, with Application to Northern California" by Felix Waldhauser and "A Double-Difference Earthquake Location Algorithm: Method and Application to the Northern Hayward Fault, California" by Waldhauser & Ellsworth.

rtDD also supports NonLinLoc by Anthony Lomax grid file format alongside the travel time formats natively supported by SeisComP (LOCSAT and libtau). See :ref:`ttt-label`.

The double-difference equation system solver uses LSQR by Chris Paige, Michael Saunders and LSMR by David Fong, Michael Saunders algorithms.


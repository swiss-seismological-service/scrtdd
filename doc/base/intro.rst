.. _intro-label:

Introduction
============

rtDD has two mode of operations: :ref:`multi-event <multi-event-label>` mode and :ref:`single-event <single-event-label>` mode. In  Multi-Event mode rtDD can relocated event catalogs using the double-difference method. This mode doesn't work in real-time and the resulting catalog can be used for post processing analysis, or it can become the reference catalog of the single event-mode, or real-time mode, where the events are relocated one at a time as they occur.

The methods developed in rtDD are based on the paper "Near-Real-Time Double-Difference Event Location Using Long-Term Seismic Archives, with Application to Northern California" by Felix Waldhauser and "A Double-Difference Earthquake Location Algorithm: Method and Application to the Northern Hayward Fault, California" by Waldhauser & Ellsworth.

rtDD also supports NonLinLoc by Anthony Lomax grid file format alongside the travel time formats natively supported by SeisComP (LOCSAT and libtau). See :ref:`ttt-label`.

The double-difference equation system solver uses LSQR by Chris Paige, Michael Saunders and LSMR by David Fong, Michael Saunders algorithms.


rtDD documentation
=======================

.. image:: base/media/example.png

This documentation covers the rtDD software, version |version|, which consists of the |appname| module and a few other scripts and resources documented within. Full release information can be found `here. <https://github.com/swiss-seismological-service/scrtdd/releases>`_

|appname| is a `SeisComP <https://seiscomp.de>`_ extension module, developed at the `Swiss Seismological Service <http://www.seismo.ethz.ch>`_. It implements Double-Difference event relocation for both real-time, event-by-event processing and classic offline mode, where an entire earthquake catalog is relocated.

.. image:: base/media/logo.png
   :width: 400
   :align: center
   :target: http://www.seismo.ethz.ch

Installation from binaries
--------------------------

Compiled versions of this module can be found at https://data.gempa.de/packages/Public/. The installation files are compressed tar archives containing the binary distribution of this module. These can be extracted into your SeisComP installation folder and used immediately.


Installation from source code
-----------------------------

The development page for the rtDD project is located at https://github.com/swiss-seismological-service/scrtdd. This page provides instructions on how to compile the software from source.


About
#####

.. toctree::
   :maxdepth: 3

   /base/intro
   /base/multievent
   /base/singleevent
   /base/xcorr
   /base/continuous
   /base/ttt
   /apps/scrtdd

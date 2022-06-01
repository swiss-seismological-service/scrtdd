rtDD documentation
=======================

.. image:: base/media/example.png

This is documentation for rtDD software in version |version|, which consists in the module |appname| and few other scripts and resources documented in here. Full release information can be found `here. <https://github.com/swiss-seismological-service/scrtdd/releases>`_

|appname| is a `SeisComP <https://seiscomp.de>`_ extension module developed at the `Swiss Seismological Service <http://www.seismo.ethz.ch>`_ that implements Double-Difference event relocation both in Real-Time, one event at the time, and classic offline mode, where an earthquake catalog is relocated as a whole.

.. image:: base/media/logo.png
   :width: 400
   :align: center
   :target: http://www.seismo.ethz.ch

Installation from binaries
--------------------------

You can find compiled version of this module at https://data.gempa.de/packages/Public/. The installation file is a compressed tar archive containing the binary distribution of this module, which can be extracted in your SeisComP installation folder and used right after.


Installation from source code
-----------------------------

The development page of the rtDD project can be found at https://github.com/swiss-seismological-service/scrtdd. There you can find instructions for the compilation of the software.


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
   /base/waveforms
   /base/database
   /apps/scrtdd

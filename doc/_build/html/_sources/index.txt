.. Caton documentation master file, created by
   sphinx-quickstart on Sun Sep  6 15:19:35 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

***********************************************
Caton: software for spike sorting
***********************************************

Caton (http://caton.googlecode.com) is a software package that automates the process of spike sorting, starting with raw data. It is specially designed to be used on data from multi-site probes with up to 32 sites or more. The software was named after Richard Caton, the 19th century Englishman who first recorded electrical signals in the brain.

I (John Schulman) wrote this software while I was working in the |Buzsaki lab|_ in the summer of 2009.

This program serves as a complete tool for the spike sorting. The outputs are in formats that can be viewed in Klusters_ and Neuroscope_.

.. _Klusters: http://klusters.sourceforge.net
.. _Neuroscope: http://neuroscope.sourceforge.net
.. |a| unicode:: U+00E1
   :trim:
.. |Buzsaki lab| replace:: Buzs |a| ki lab
.. _`Buzsaki lab`: http://osiris.rutgers.edu/frontmid/indexmid.html

User's guide
=================================

.. toctree::
   :maxdepth: 2
   
   installation.rst
   usage.rst
   probe_files.rst
   
   
Technical documentation
==============================

.. toctree::
   :maxdepth: 2

   detect_extract.rst
   clustering.rst
   benchmarks.rst

References
==========================

.. toctree::

   refs.rst
   
Appendix
===============

.. toctree::
   :maxdepth: 2

   other_scripts.rst
   mog.rst
   why_caton.rst


* :ref:`search`


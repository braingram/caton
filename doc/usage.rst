******
Usage
******

Data file format
=================

The input and output data are the same as the formats used by Klusters and Neuroscope. These are described at http://neuroscope.sourceforge.net/UserManual/data-files.html


Input to ``cluster_from_raw_data.py`` is a binary file, which is a contiguous block of 16-bit integers::
     
     sample 1, channel 1
     sample 1, channel 2
     ...
     sample 1, channel C
     sample 2, channel 1
     sample 2, channel 2
     ...
     sample 2 channel C
     ...

Overview
===============

There are two different processes you can use to sort your data.

      1. Batch

      	 * Do a large batch job to process all of your data at once.

      2. Batch, then generalize:

      	 * Process a subset of your data to generate a clustering model (a bunch of clusters and their parameters),
	 * Possibly modify that model manually (e.g. with Klusters_)
	 * Process the rest of your data, and use this model to classify the spikes.

Batch
-------

	* Navigate to the directory containing your ``.dat`` file.
  	* Make one or more ``.probe`` files in the directory of the ``.dat`` file (see :ref:`probefiles`).
	* Type the following::
	    
		cluster_from_raw_data.py your_dat_file.dat
	    
	* The output will be in the directory ``your_dat_file_tetrode_batch/`` (or something like that, depending on the name of your probe file.) If there are multiple probe files, there will be one output directory for each probe file in the directory of your ``.dat`` file.


Generalize
-----------

**NOT CURRENTLY FUNCTIONAL**

	* Follow the batch instructions above, except restrict the number of spikes extracted::

	  	 cluster_from_raw_data.py your_dat_file.dat -n 100000

	* Then run the generalize script, and tell it what directory to look in for the results of the batch clustering.::
	       	 
		 generalize_from_raw_data.py your_dat_file.dat

	* The output will be in ``your_dat_file_tetrode_generalize/``




.. _usage:

Script options reference
=========================

The following are the outputs of ``[script name].py -h``

cluster_from_raw_data.py
------------------------

.. literalinclude:: generated_data/cluster_from_raw_data.help.txt


generalize_from_raw_data.py
----------------------------

.. literalinclude:: generated_data/generalize_from_raw_data.help.txt


.. _Klusters: http://klusters.sourceforge.net

.. _probefiles:

Probe files
============
To cluster a data file, you must create a ``.probe`` file, which tells the program about the channel mapping and the probe layout.

Each line of the probe file describes one recording site. The lines can have three different formats, whether the probe is zero-dimensional, one-dimensional, or two-dimensional. What follows are sample probe files. The parenthesized expressions give x and y coordinates (scale is irrelevant.)

You can use the script ``plot_probe.py`` to check that the layout of the sites is correct and the adjacency graph (represented by blue lines) is reasonable.

Zero-dimensional probe
-----------------------

Format: ``name	dat-channel``

``caton/probes/tetrode.probe``:

.. literalinclude:: ../probes/tetrode.probe

.. figure:: generated_data/tetrode.png
   :width: 10cm
   :align: center

One-dimensional probe
----------------------

Format: ``name	dat-channel	(depth)``

``caton/probes/mar16.probe``:

.. literalinclude:: ../probes/mar16.probe

.. figure:: generated_data/mar16.png
   :width: 15cm
   :align: center

Two-dimensional probe
----------------------

Format: ``name	dat-channel	(x y)``

``caton/probes/mar32.probe``

.. literalinclude:: ../probes/mar32.probe

.. figure:: generated_data/mar32.png
   :width: 6cm
   :align: center




Advanced probe file usage
--------------------------

The probe file determines two things:
1. The edges in the adjacency graph of the sites (blue lines in figures above)
2. The groups of nearby sites that will be clustered together

You may notice that after you run `cluster_from_raw_data.py` or `plot_probe.py` (any script that uses the probe file), the probe file will have a bunch of extra lines in it. Here is the zero-dimensional probe above--the tetrode--after this modification:

.. literalinclude:: ../probes/tetrode2.probe

You may want to change the edges and the groups manually. Simply add or remove lines. For example, if you want to sort all of your channels separately, you can remove all of the edge lines and put the sites each in their own group.

.. literalinclude:: ../probes/tetrode3.probe


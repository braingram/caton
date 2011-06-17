**********************
Other included scripts
**********************


check_crosstalk.py
-------------------

This script gives an easy way of figuring out if there is crosstalk between your channels. If there is crosstalk between two faraway channels, you should remove both of them (delete those lines from the .probe file)

Command-line usage:

.. literalinclude:: generated_data/check_crosstalk.help.txt

This script plots a colormap of the inverse-covariance matrix of a segment of a few seconds of filtered data. Actually, it takes the above-diagonal portion of the matrix, and zeros the diagonal and the below-diagonal part. (This is because the diagonal entries are by far the largest.)

The following comes from a 32-site probe with two rows of channels:

.. figure:: generated_data/mar32_before_crosstalk.png
   :width: 15cm

One matrix element (between SP25 and SP27) dwarves all of the others. We delete these two lines from the probe file (or comment them out with ``#``)

.. figure:: generated_data/mar32_after_crosstalk.png
   :width: 15cm

The blue dots indicate pairs of sites that are considered spatially adjacent when the channel adjacency graph is generated (see :ref:`adjacency`).

plot_probe.py
---------------

Plots the sites of the probe and the channel adjacency graph. Used to create the probe plots throught this document.

.. literalinclude:: generated_data/plot_probe.help.txt


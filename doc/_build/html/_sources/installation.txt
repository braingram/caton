*************
Installation
*************

Caton should work on any `*nix` operating system (including Linux and Mac OS X). Caton has a number of dependencies. 


1. Python and packages:

   * python 2.5+
   * numpy
   * scipy
   * matplotlib
   * cython
   * pytables

2. KlustaKwik

First follow the instructions on how to install the python dependencies for your particular system. Then follow the KlustaKwik installation instructions. 


If you have no experience installing python packages, or you don't want to use sudo, the easiest option is to use sage_ for the python dependencies. Otherwise you will have to install the packages individually. Last, follow the installation instructions for Caton itself.

.. _downloads: http://code.google.com/p/caton/downloads/
.. _Klustakwik: http://klustakwik.sourceforge.net



Python dependencies on Mac
=======================================
All dependencies are included in the `Enthought Python Distribution <http://www.enthought.com/products/getepd.php>`_

Python dependencies on Debian/Ubuntu Linux
===============================================

First install a bunch of packages::

   sudo apt-get install python-numpy python-scipy python-matplotlib cython python-setuptools
   sudo apt-get build-dep python-tables
   sudo easy_install pyrex
   sudo easy_install tables

If you get an error about ``utilsExtension`` try this::

   easy_install http://www.pytables.org/download/preliminary/pytables-2.2b3/tables-2.2b3.tar.gz

.. _sageinst:

Python dependencies with sage
==============================

sage_ is a self-contained distribution with various math/science python packages. This should work on any `*nix` operating system.

.. _sage: http://www.sagemath.org/

If hdf5 is not installed on your computer, you will need to install it.
1. Download source from http://www.hdfgroup.org/ftp/HDF5/current16/src
2. Unpack it, compile, and install::
   
   tar -xvf hdf5-1.6.10.tar.gz
   ./configure
   make
   make install
   
This will produce directory hdf5-1.6.10/hdf5

3. Set an environment variable to tell pytables where hdf5 libraries are during installation::

    export HDF5_DIR=/path/to/hdf5-1.6.10/hdf5

4. I had to do one more thing, otherwise I got "ImportError: libhdf5.so.0: cannot open shared object file". Put the following line in your `sage` script (it's in the top-level sage directory)::

    export LD_LIBRARY_PATH=/path/to/hdf5-1.6.10/hdf5/lib

**Important**: before you install or use caton, first enter the sage directory and type::

    ./sage -sh

This starts a *subshell* where all environment variables (e.g. paths) are set, and programs will use sage's python instead of your system's default python.

After starting the subshell, type::

    easy_install pyrex
    easy_install http://www.pytables.org/download/preliminary/pytables-2.2b3/tables-2.2b3.tar.gz


Installing KlustaKwik
=======================

If you are on Debian/Ubuntu Linux, install the .deb file on the downloads page.

Otherwise, follow these steps:

1. Download the archive from http://klustakwik.sourceforge.net and unarchive it

2. Compile it by typing `make`

3. Move the compiled binary to a directory on your PATH (e.g. /usr/local/bin, or /home/username/local/bin, or sage/local/bin).

Installing Caton
==================

1. Download the most recent source tarball from the downloads_ page and unarchive it.

2. Go to the Caton directory and type::

      sudo python setup.py install

If you used sage, omit the sudo since Caton is being installed in your sage directory.


Overview
========

ChaosMagPy is a simple Python package for evaluating the
`CHAOS <http://www.spacecenter.dk/files/magnetic-models/CHAOS-6/>`_ geomagnetic
field model. To get quickly started, download a complete working example
including the latest model `here <http://www.spacecenter.dk/files/magnetic-models/CHAOS-6/chaosmagpy_package.zip>`_.

References
----------

*Finlay, C.C., Olsen, N., Kotsiaros, S., Gillet, N. and Toeffner-Clausen, L.
(2016)*, Recent geomagnetic secular variation from Swarm and ground observatories
as estimated in the CHAOS-6 geomagnetic field model Earth Planets Space,
Vol 68, 112. doi: 10.1186/s40623-016-0486-1

*Olsen, N., Luehr, H., Finlay, C.C., Sabaka, T. J., Michaelis, I., Rauberg, J.
and Toeffner-Clausen, L. (2014)*, The CHAOS-4 geomagnetic field model,
Geophys. J. Int., Vol 197, 815-827, doi: 10.1093/gji/ggu033.

*Olsen, N.,  Luehr, H.,  Sabaka, T.J.,  Mandea, M. ,Rother, M., Toeffner-Clausen, L.
and Choi, S. (2006)*, CHAOS — a model of Earth's magnetic field derived from CHAMP,
Ørsted, and SAC-C magnetic satellite data, Geophys. J. Int., vol. 166 67-75

Documentation
-------------

The latest documentation of the package is available on
`Read the Docs <https://chaosmagpy.readthedocs.io/en/latest/>`_.

License
=======

MIT License

Copyright (c) 2019 Clemens Kloss

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Installation
============

ChaosMagPy relies on the following:

* python>=3.6
* numpy
* scipy
* pandas
* cython
* cartopy>=0.17
* matplotlib>=3
* cdflib
* hdf5storage

Specific installation steps using the conda/pip package managers are as follows:

1. Install packages with conda:

   >>> conda install python numpy scipy pandas cython cartopy matplotlib

2. Install packages with pip:

   >>> pip install cdflib hdf5storage

3. Finally install ChaosMagPy either with pip from PyPI:

   >>> pip install chaosmagpy

   or, if you have downloaded the `package files <https://pypi.org/project/chaosmagpy/#files>`_
   to the current working directory, with:

   >>> pip install chaosmagpy-x.x-py3-none-any.whl

   or, alternatively

   >>> pip install chaosmagpy-x.x.tar.gz

   replacing ``x.x`` with the correct version.

Contents
========

The directory contains the files/directories:

1. "chaosmagpy-x.x*.tar.gz": pip installable archive of the chaosmagpy package
   (version x.x*)

2. "chaos_examples.py": executable Python script containing several examples
   that can be run by changing the examples in line 16, save and run in the
   command line:

   >>> python chaos_examples.py

   example 1: Calculate CHAOS model field predictions from input coordinates
              and time and output simple data file
   example 2: Calculate and plot residuals between CHAOS model and
              Swarm A data (from L1b MAG cdf data file, example from May 2014).
   example 3: Calculate core field and its time derivatives for specified times
              and radii and plot maps
   example 4: Calculate static (i.e. small-scale crustal) magnetic field and
              plot maps (may take a few seconds)
   example 5: Calculate timeseries of the magnetic field at a ground
              observatory and plot
   example 6: Calculate external and associated induced fields described in SM
              and GSM reference systems and plot maps

3. "data/CHAOS-6-x9.mat": mat-file containing CHAOS-6 model (extension 9)

4. "SW_OPER_MAGA_LR_1B_20180801T000000_20180801T235959_PT15S.cdf":
   cdf-file containing Swarm A magnetic field data from August 1, 2018.

5. directory called "html" containing the built documentation as
   html-files. Open "index.html" in your browser to access the main site.


Clemens Kloss (ancklo@space.dtu.dk)


Changelog
=========

Version 0.1
-------------
| **Date:** May 10, 2019
| **Release:** v0.1

Features
^^^^^^^^
* New CHAOS class method ``synth_euler_angles`` to compute euler angles for
  the satellites from the CHAOS model (used to rotate vectors from
  magnetometer frame to the satellite frame).
* Added CHAOS class methods ``synth_values_tdep``, ``synth_values_static``,
  ``synth_values_gsm`` and ``synth_values_sm`` for field value computation.
* RC index file now stored in HDF5 format.
* Filepaths and other parameters are now handled by ``basicConfig``.
* Added extrapolation keyword to CHAOS class method ``synth_coeffs``, linear by
  default.
* ``mjd2000`` now also accepts datetime class instances.
* ``load_RC_datfile`` downloads latest RC-index file from the website if no
  file is given.

Bugfixes
^^^^^^^^
* Resolved issue in ``model_utils.degree_correlation``.
* Changed the date conversion to include hours and seconds not just the day
  when plotting the timeseries.

Version 0.1a3
-------------
| **Date:** February 19, 2019
| **Release:** v0.1a3

Features
^^^^^^^^
* New CHAOS class method ``save_matfile`` to output MATLAB compatible
  files of the CHAOS model (using the ``hdf5storage`` package).
* Added ``epoch`` keyword to basevector input arguments of GSM, SM and MAG
  coordinate systems.

Bugfixes
^^^^^^^^
* Fixed problem of the setup configuration for ``pip`` which caused importing
  the package to fail although installation was indicated as successful.

Version 0.1a2
-------------
| **Date:** January 26, 2019
| **Release:** v0.1a2

Features
^^^^^^^^
* ``mjd_to_dyear`` and ``dyear_to_mjd`` convert time with microseconds
  precision to prevent round-off errors in seconds.
* Time conversion now uses built-in ``calendar`` module to identify leap year.

Bugfixes
^^^^^^^^
* Fixed wrong package requirement that caused the installation of
  ChaosMagPy v0.1a1 to fail with ``pip``. If installation of v0.1a1 is needed,
  use ``pip install --no-deps chaosmagpy==0.1a1`` to ignore faulty
  requirements.


Version 0.1a1
-------------
| **Date:** January 5, 2019
| **Release:** v0.1a1

Features
^^^^^^^^
* Package now supports Matplotlib v3 and Cartopy v0.17.
* Loading shc-file now converts decimal year to ``mjd2000`` taking leap years
  into account by default.
* Moved ``mjd2000`` from ``coordinate_utils`` to ``data_utils``.
* Added function to compute degree correlation.
* Added functions to compute and plot the power spectrum.
* Added flexibility to the function synth_values: now supports NumPy
  broadcasting rules.
* Fixed CHAOS class method synth_coeffs_sm default source parameter: now
  defaults to ``'external'``.

Deprecations
^^^^^^^^^^^^
* Optional argument ``source`` when saving shc-file has been renamed to
  ``model``.
* ``plot_external_map`` has been renamed to ``plot_maps_external``
* ``synth_sm_field`` has been renamed to ``synth_coeffs_sm``
* ``synth_gsm_field`` has been renamed to ``synth_coeffs_gsm``
* ``plot_static_map`` has been renamed to ``plot_maps_static``
* ``synth_static_field`` has been renamed to ``synth_coeffs_static``
* ``plot_tdep_maps`` has been renamed to ``plot_maps_tdep``
* ``synth_tdep_field`` has been renamed to ``synth_coeffs_tdep``


Version 0.1a0
-------------
| **Date:** October 13, 2018
| **Release:** v0.1a0

Initial release to the users for testing.

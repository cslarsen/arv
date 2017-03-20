arv — a fast 23andMe parser for Python
======================================
|travis-status| |versions| |license| |pypi|

Arv (Norwegian; "heritage" or "inheritance") is a Python module for parsing raw
23andMe genome files. It lets you lookup SNPs from RSIDs.

.. code:: python

  from arv import load, unphased_match as match

  genome = load("genome.txt")

  print("You are a {gender} with {color} eyes and {complexion} skin.".format(
    gender     = "man" if genome.y_chromosome else "woman",
    complexion = "light" if genome["rs1426654"] == "AA" else "dark",
    color      = match(genome["rs12913832"], {"AA": "brown",
                                              "AG": "brown or green",
                                              "GG": "blue"})))

For my genome, this little program produces::

    You are a man with blue eyes and light skin.

The parser is insanely fast, having been written in finely tuned C++, exposed
via Cython. A 2013 Xeon machine I've tested on parses a 24 Mb file into a hash
table in 70 ms!

Works with Python 2.7+ and 3+. Installable with pip!

Installation
============

The recommended way is to install from PyPi.

.. code:: bash

    $ pip install arv

This will most likely build Arv from source. The package requires Cython, but
it doesn't check if you have a C++ compiler. Currently, it expects that you
have clang++ or g++.

If you have problems running ``pip install arv``, please open an issue on
GitHub with as much detail as possible (``g++/clang++ --version``, ``uname
-a``, ``python --version`` and so on).

If you set the environment variable ``ARV_DEBUG``, it will build with full
warnings and debug symbols.

License
=======

Copyright 2017 Christian Stigen Larsen

Distributed under the GNU GPL v3 or later. See the file COPYING for the full
license text. This software makes use of open source software; see LICENSES for
details.

.. |travis-status| image:: https://travis-ci.org/cslarsen/arv.svg?branch=master
    :alt: Travis build status
    :scale: 100%
    :target: https://travis-ci.org/cslarsen/arv

.. |license| image:: https://img.shields.io/badge/license-GPL%20v3%2B-blue.svg
    :target: http://www.gnu.org/licenses/old-licenses/gpl-3.en.html
    :alt: Project License

.. |versions| image:: https://img.shields.io/badge/python-2%2B%2C%203%2B-blue.svg
    :target: https://pypi.python.org/pypi/arv/
    :alt: Supported Python versions

.. |pypi| image:: https://badge.fury.io/py/arv.svg
    :target: https://badge.fury.io/py/arv

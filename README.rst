arv — a fast 23andMe parser for Python
======================================
|travis-status| |versions| |license|

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

Note that the pip install builds from source. You'll need not only Cython, but
also a C++11 capable compiler. I might distribute binary wheels in time. If the
installation doesn't work for you, please file a GitHub issue with as much
detail as you can.

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

arv — a fast 23andMe parser for Python
======================================

Arv (Norwegian; "heritage" or "inheritance") is a Python module for parsing raw
23andMe genome files. It lets you lookup SNPs from RSIDs.

.. code:: python

  from arv import load, unphased_match as match

  genome = load("genome.txt")

  print("You are a {gender} with {color} eyes and {complexion} skin.".format(
    gender     = "man" if genome.y_chromosome else "woman",
    complexion = "light" if genome["rs1426654"] == "AA" else "dark",
    eyecolor   = match(genome["rs12913832"], {"AA": "brown",
                                              "AG": "brown or green",
                                              "GG": "blue"})))

For my genome, this little program produces::

    You are a man with blue eyes and light skin.

The parser is insanely fast, having been written in finely tuned C++. A Xeon
machine I've tested on parses a 24 Mb file into a hash table in 70 ms.

Works with Python 2 and 3.

Installation
============

The recommended way is to install from PyPi.

.. code:: bash

    $ pip install arv

Note that the pip install builds from source. You'll need not only Cython, but
also a C++11 capable compiler. I might distributed binary wheels in time. If
the installation doesn't work for you, please file a GitHub issue with as much
detail as you can.

License
=======

Copyright 2017 Christian Stigen Larsen

Distributed under the GNU GPL v3 or later.

See the file COPYING for the full license text. This software makes use of open
source software; see LICENSES for details.

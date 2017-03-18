arv — a fast 23andMe parser for Python
======================================

Arv (Norwegian; "inheritance" or "heritage") is a Python module for parsing raw
23andMe genome files. It lets you lookup SNPs from RSIDs::

    from arv import load, unphased_match as match

    genome = load("genome.txt")

    print("You are a {gender} with {color} eyes and {complexion} skin.".format(
      gender     = "man" if genome.y_chromosome else "woman",
      complexion = "light" if genome["rs1426654"] == "AA" else "dark",
      eyecolor   = match(genome["rs12913832"], {"AA": "brown",
                                                "AG": "brown or green",
                                                "GG": "blue"})))

In my case, this little program produces::

    You are a man with blue eyes and light skin.

It's insanely fast: On a 2013 Xeon machine, a 24 Mb file is fully parsed and
put into a hash table in less than 70 ms. Its guts are written in finely tuned
C++ and is exposed to Python via Cython.

Status
======

> This project is currently just a work in progress! I intend to wrap dna-traits
under a new name, using Cython to interface with dna-traits much more easily.
>
> For a working (but old) version, see https://github.com/cslarsen/dna-traits*

Installation
============

The recommended way is to install from PyPi::

    $ pip install arv

**NOTE: This is not yet available!**

You can also build from source. In that case, you need a C++ compiler with
C++11 support, CMake and Cython::

    $ mkdir build; cd build
    $ cmake .. -GNinja \
      -DARV_PYTHON_VERSION=2 -DCMAKE_INSTALL_PREFIX=../install
    $ ccmake . # optional, for changing settings
    $ ninja

You can change which Python version to compile for with the CMake cache
variable `ARV_PYTHON_VERSION`, for example `cmake -DARV_PYTHON_VERSION=2
[...]`. If you have several versions of Python installed on your machine, e.g.
3.3 and 3.4 and so on, you are advised to specify both major and minor version
here, for example `ARV_PYTHON_VERSION=3.3`.

Requirements
------------

  * A C++ compiler with C++11 support
  * A UNIX system
  * Optional: CMake 3.5 or newer
  * Cython: `pip install cython`

License
=======

Copyright 2014, 2016, 2017 Christian Stigen Larsen  
Distributed under the GNU GPL v3 or later.

See the file COPYING for the full license text. This software makes use of open
source software; see LICENSES for details.

This code is largely based on the project dna-traits, by the same author.

arv — a fast 23andMe parser for Python
======================================

Arv (Norwegian; "inheritance") is a Python module for parsing raw 23andMe
genome files. It lets you lookup SNPs from RSIDs:

    import arv

    genome = arv.load("genome.txt")

    print("You are a {gender} with {eyecolor} and {complexion} skin.".format(
      gender     = "man" if genome.y_chromosome else "woman",
      complexion = "light" if genome["rs1426654"] == "AA" else "dark",
      eyecolor   = arv.unphased_match(genome["rs12913832"], {
                      "AA": "brown eyes",
                      "AG": "brown or green eyes",
                      "GG": "blue eyes"})))

In my case, this little program produces

    You are a man with blue eyes and light skin.

It's insanely fast: On a 2013 Xeon machine, a 24 Mb file is fully
parsed and put into a hash table in less than 70 ms. Its guts are written in
finely tuned C++ and is exposed to Python via Cython.

Status
======

*This project is currently just a work in progress! I intend to wrap dna-traits
under a new name, using Cython to interface with dna-traits much more easily.*

*For a working (but old) version, see https://github.com/cslarsen/dna-traits*

Installation
============

The recommended way is to install from PyPi:

    $ pip install arv

*NOTE: This is not yet available!*

You can also build from source. In that case, you need a C++ compiler with
C++11 support, CMake and Cython:

    $ mkdir build; cd build
    $ cmake .. -GNinja \
      -DARV_PYTHON_VERSION=2 -DCMAKE_INSTALL_PREFIX=../install
    $ ccmake . # optional, for changing settings
    $ ninja

You can change which Python version to compile for with the CMake cache
variable `ARV_PYTHON_VERSION`, for example `cmake -DARV_PYTHON_VERSION=2
[...]`.

Requirements
------------

  * A C++ compiler
  * A UNIX system
  * CMake 3.5 or newer
  * Cython: `pip install cython`

License
=======

Copyright 2014, 2016, 2017 Christian Stigen Larsen  
Distributed under the GNU GPL v3 or later.

See the file COPYING for the full license text. This software makes use of open
source software; see LICENSES for details.

This code is largely based on the project dna-traits, by the same author.

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

.. code:: bash

    $ pip install --upgrade arv

See below for software requirements.

Important disclaimer
====================

It's very important to tell you that I, the author of arv, am merely a
*hobbyist*! I'm a professional software developer, not a geneticist, medical
doctor or anything like that.

Because of that, this software may not only look weird to people in the field,
it may also contain serious errors. If you find any problem whatsoever, please
submit a GitHub issue!

This a slightly modified version of what I wrote for the original software
called "dna-traits", and the same goes this software:

In addition to the GPL v3 license terms, and given that this code deals with
health-related issues, I want to stress that the provided code most likely
contains errors, or invalid genome reports. Results from this code must be
interpreted as HIGHLY SPECULATIVE and may even be downright INCORRECT. Always
consult an expert (medical doctor, geneticist, etc.) for guidance. I take NO
RESPONSIBILITY whatsoever for any consequences of using this code, including
but not limited to loss of life, money, spouses, self-esteem and so on. Use at
YOUR OWN RISK.

The indended use is for casual, educational purposes. If this code is used for
research purposes, please cross-check key results with other software: The
parser code may contain serious errors, for example.

An interesting story about the research part: I once released a pretty good
Mersenne Twister PRNG that ended up being used in research. Turned out the
engine had bugs, and by the time I had fixed them, the poor researcher already
had results (hopefully not published; I don't know). The guy had to go back and
fix hi stuff, and I felt terribly bad about it.

So beware!

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

Usage
=====

First you need to dump the raw genome file from 23andMe. You'll find it under
the raw genome browser, and download the file. You may have to unzip it first:
The parser works on the pure text files.

Then you load the genome in Python with

... code:: python

    >>> genome = arv.load("filename.txt")
    >>> genome
    <Genome: SNPs=960614, name='filename.txt'>

To see if there are any Y-chromosomes present in the genome,

.. code:: python

    >>> genome.y_chromosome
    True

The genome provides a ``dict``-like interface. To get the genotype of a given SNP, just enter the RSID. It will return it as a string.

.. code:: python

    >>> genome["rs123"]
    'AA'

You can also access the SNP as an object:

.. code:: python

    >>> genome.get_snp("rs123")
    >>> snp
    <SNP: chromosome=7 position=24966446 genotype='AA'>
    >>> snp.chromosome
    7
    >>> snp.position
    24966446
    >>> snp.genotype
    'AA'

The last line actually returns a ``PyGenotype`` object, but its ``repr``
returns something that looks like a string. This lets you perform a few
operations on the nucleotides. For example, you can get its complement with the
``~``-operator.

.. code:: python

    >>> type(snp.genotype)
    <type '_arv.PyGenotype'>
    >>> snp.genotype
    'AA'
    >>> ~snp.genotype
    'TT'

The complement is important due to eah SNPs orientation. All of 23andMe SNPs
are oriented towards the positive ("plus") strand, based on the GRCh37
reference human genome assembly build. But some SNPs on SNPedia are given with
the `minus orientation <http://snpedia.com/index.php/Orientation>`.

For example, to determine if the human in question is likely lactose tolerant
or not, we can look at `rs4988235 <http://snpedia.com/index.php/Rs4988235>`.
SNPedia reports its _Stabilized_ orientation to be minus, so we need to use the
complement:

.. code:: python

    >>> genome.get_snp("rs4988235").genotype
    'AA'
    >>> ~genome.get_snp("rs4988235").genotype
    'TT'

By reading a few `GWAS <>` research papers, we can build a rule to determine a
human's likelihood for lactose tolerance:

.. code:: python

    >>> arv.unphased_match(~genome.get_snp("rs4988235").genotype, {
        "TT": "Likely lactose tolerant",
        "TC": "Likely lactose tolerant",
        "CC": "Likely lactose intolerant",
        None: "Unable to determine (genotype not present)"})

Note that for non-professionals, reading GWAS papers can be a bit tricky. . To
create a
usually requi
Note that reading GWAS papers can be a bit tricky, and it is very tempting to
jump to conclusions or not understand all of it. But for casual users, it's
very educational and most of all very fun to try your hand at inferring various
results.

Command line interface
======================

You can also invoke ``arv`` from the command line:

.. code:: bash

		$ python -m arv --help

For example, you can drop into a Python REPL like so:

.. code:: bash

		$ python -m arv --repl genome.txt
		genome.txt ... 960614 SNPs, male
		Type `genome` to see the parsed 23andMe raw genome file
		>>> genome
		<Genome: SNPs=960614, name='genome.txt'>
		>>> genome["rs123"]
		'AA'

If you specify several files, you can access them through the variable
``genomes``.

The example at the top of this document can be run with ``--example``:

.. code:: bash

		$ python -m arv --example genome.txt
		genome.txt ... 960614 SNPs, male

		genome.txt ... A man with blue eyes and light skin

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

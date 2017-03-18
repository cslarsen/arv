"""
A fast 23andMe raw genome file parser.

To load a genome file,

    >>> import arv
    >>> genome = arv.load("genome.txt")

You can then look up genotypes from RSIDs

    >>> genome["rs123"]
    'AA'

You can also access SNPs

    >>> genome.get_snp("rs123")
    <SNP: chromosome=1 position=112233 genotype='AA'>

By using ``unphased_match``, you can match genotypes while disregarding the
ordering of the two nucleotides. For example, ``AT`` and ``TA`` would be
considered equal. Here is an example usage:

        genotype = genome["rs12913832"]
        color = unphased_match(igenotype, {
                    "AA": "brown",
                    "AG": "brown or green",
                    "GG": "blue",
                    None: "unknown"})

A full example would be

    import arv

    genome = arv.load("genome.txt")
    print("You are a {gender} with {color} eyes and {complexion} skin.".format(
        gender     = "man" if genome.y_chromosome else "woman",
        complexion = "light" if genome["rs1426654"] == "AA" else "dark",
        eyecolor   = unphased_match(genome["rs12913832"], {
                        "AA": "brown",
                        "AG": "brown or green",
                        "GG": "blue"})))

Copyright 2014, 2016, 2017 Christian Stigen Larsen
Distributed under the GPL v3 or later.
"""

from _arv import (
    _sizes,
    load,
    PyGenome as Genome,
    PyGenotype as Genotype,
    PySNP as SNP,
)

from .match import unphased_match
from . import traits
from . import util

__author__ = "Christian Stigen Larsen"
__copyright__ = "Copyright 2014, 2016, 2017 Christian Stigen Larsen"
__email__ = "csl@csl.name"
__license__ = "GNU General Public License v3"
__version__ = "0.1"

__all__ = [
    "_sizes",
    "Genome",
    "Genotype",
    "load",
    "SNP",
    "unphased_match",
]

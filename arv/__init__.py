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
        eyecolor = unphased_match(genotype, {
                        "AA": "brown",
                        "AG": "brown or green",
                        "GG": "blue",
                        None: "unknown"})

A full example would be:

    import arv

    genome = arv.load("genome.txt")

    print("You are a {gender} with {color} eyes and {complexion} skin.".format(
        gender     = "man" if genome.y_chromosome else "woman",
        complexion = "light" if genome["rs1426654"] == "AA" else "dark",
        color      = unphased_match(genome["rs12913832"], {
                        "AA": "brown",
                        "AG": "brown or green",
                        "GG": "blue"})))

For a given genome, this might print

    You are a man with blue eyes and light skin.

Copyright 2014, 2016, 2017 Christian Stigen Larsen
Distributed under the GNU GPL v3 or later.
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
__copyright__ = "Copyright 2017 Christian Stigen Larsen"
__credits__ = ["Christian Stigen Larsen", "Google"]
__email__ = "csl@csl.name"
__license__ = "GNU General Public License v3"
__maintainer__ = "Christian Stigen Larsen"
__status__ = "Prototype"
__version__ = "0.4"

__all__ = [
    "_sizes",
    "Genome",
    "Genotype",
    "load",
    "SNP",
    "unphased_match",
]

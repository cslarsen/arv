"""
A fast 23andMe raw genome file parser.

Copyright 2014, 2016, 2017 Christian Stigen Larsen
Distributed under the GPL v3 or later.
"""

from _arv import (
    PyGenome as Genome,
    load,
)

from match import unphased_match

__author__ = "Christian Stigen Larsen"
__copyright__ = "Copyright 2014, 2016, 2017 Christian Stigen Larsen"
__email__ = "csl@csl.name"
__license__ = "GNU General Public License v3"
__version__ = "1.0"

__all__ = [
    "Genome",
    "load",
    "unphased_match",
]

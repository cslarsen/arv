# Copyright 2017 Christian Stigen Larsen
# Distributed under the GNU GPL v3 or later; see COPYING.

from libcpp cimport bool
from libcpp.string cimport string

cdef extern from "dnatraits.hpp":
    cdef cppclass Genome:
        Genome() except +
        Genome(const size_t) except +
        Genome(const Genome&) except +
        Genome& operator=(const Genome&)
        double load_factor() const
        size_t size() const
        bool y_chromosome

    cdef void parse_file(const string&, Genome&) except +

cdef class PyGenome:
    cdef Genome _genome

    def __cinit__(self, size_t size=1000003):
        self._genome = Genome(size)

    def load_factor(self):
        return self._genome.load_factor()

    def __len__(self):
        return self._genome.size()

    def ychromo(self):
        return self._genome.y_chromosome

def load(string filename):
    """Loads given 23andMe raw genome file.

    Arguments:
        filename: Name of file to load.

    Raises:
        RuntimeError

    Returns:
        A PyGenome.
    """
    genome = PyGenome(1000003)
    parse_file(filename, genome._genome)
    return genome

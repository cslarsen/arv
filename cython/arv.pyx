# arv
# Copyright 2017 Christian Stigen Larsen
# Distributed under the GNU GPL v3 or later; see COPYING.

from libc.stdint cimport uint32_t
from libcpp cimport bool
from libcpp.string cimport string

cdef extern from "dnatraits.hpp":
    ctypedef uint32_t Position
    ctypedef uint32_t RSID

    cdef cppclass Genome:
        Genome() except +
        Genome(const size_t) except +
        Genome(const Genome&) except +
        Genome& operator=(const Genome&)
        double load_factor() const
        size_t size() const
        bool y_chromosome
        string genotype(const RSID& id) const

    cdef void parse_file(const string&, Genome&) except +

cdef class PyGenome:
    cdef Genome _genome

    def __cinit__(PyGenome self, size_t size=1000003):
        self._genome = Genome(size)

    cpdef double load_factor(PyGenome self):
        return self._genome.load_factor()

    def __len__(PyGenome self):
        return self._genome.size()

    def __getitem__(PyGenome self, key):
        cdef RSID rsid
        if isinstance(key, str):
            rsid = int(key[2:])
        elif isinstance(key, int):
            rsid = key
        else:
            raise TypeError("Expected str or int but got %s" %
                    type(key).__name__)
        cdef string s = self._genome.genotype(rsid)
        if s.empty():
            raise KeyError(key)
        else:
            return s

    @property
    def ychromo(PyGenome self):
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

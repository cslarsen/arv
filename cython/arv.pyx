# Copyright 2017 Christian Stigen Larsen
# Distributed under the GNU GPL v3 or later; see COPYING.

cdef extern from "dnatraits.hpp":
    cdef cppclass Genome:
        Genome()
        Genome(const size_t)
        Genome(const Genome&)
        double load_factor() const
    cdef Genome parse(const char* filename)

cdef class PyGenome:
    cdef Genome _genome

    def __cinit__(self, size_t size=1000000):
        self._genome = Genome(size)

    def load_factor(self):
        return self._genome.load_factor()

    @staticmethod
    def load(filename):
        cdef Genome genome = parse(filename)
        r = PyGenome()
        r._genome = genome
        return r

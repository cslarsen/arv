# arv
# Copyright 2017 Christian Stigen Larsen
# Distributed under the GNU GPL v3 or later; see COPYING.

from libc.stdint cimport uint32_t
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "dnatraits.hpp":
    ctypedef uint32_t Position
    ctypedef uint32_t RSID

    cdef enum Nucleotide:
        NONE, A, G, C, T, D, I

    cdef enum Chromosome:
        NO_CHR=0,
        CHR1, CHR2, CHR3, CHR4, CHR5,
        CHR6, CHR7, CHR8, CHR9, CHR10,
        CHR11, CHR12, CHR13, CHR14, CHR15,
        CHR16, CHR17, CHR18, CHR19, CHR20,
        CHR21, CHR22, CHR_MT, CHR_X, CHR_Y

    cdef cppclass Genotype:
        Nucleotide first
        Nucleotide second

        Genotype(const Nucleotide&, const Nucleotide&)
        bool operator==(const Genotype&) const
        bool operator<(const Genotype&) const

    cdef cppclass SNP:
        Chromosome chromosome
        Position position
        Genotype genotype

        SNP(const Chromosome& = NO_CHR,
            const Position& = 0,
            const Genotype& = Genotype(NONE, NONE))
        SNP(const SNP&)

        SNP& operator=(const SNP&);
        bool operator!=(const SNP&) const;
        bool operator<(const SNP&) const;
        bool operator<=(const SNP&) const;
        bool operator==(const Genotype&) const;
        bool operator==(const SNP&) const;
        bool operator>(const SNP&) const;
        bool operator>=(const SNP&) const;

    cdef cppclass RsidSNP:
        RSID rsid
        SNP snp
        bool operator==(const RsidSNP&) const

    cdef cppclass GenomeIterator:
        GenomeIterator(const GenomeIterator&);
        GenomeIterator& operator=(const GenomeIterator&);
        GenomeIterator& operator++();
        bool operator==(const GenomeIterator&);
        bool operator!=(const GenomeIterator&);
        const RsidSNP operator*();

    cdef const SNP NONE_SNP

    cdef cppclass Genome:
        Genome() except +
        Genome(const size_t) except +
        Genome(const Genome&) except +
        Genome& operator=(const Genome&)

        double load_factor() const
        size_t size() const

        string genotype(const RSID& id) const
        const SNP& operator[](const RSID&) const
        bool has(const RSID&) const
        void insert(const RSID&, const SNP&)
        vector[RSID] rsids() const
        vector[SNP] snps() const

        bool operator==(const Genome&) const
        bool operator!=(const Genome&) const

        GenomeIterator begin() const
        GenomeIterator end() const

        bool y_chromosome
        RSID first
        RSID last

    cdef void parse_file(const string&, Genome&) except +

cdef class PyGenome:
    cdef Genome _genome
    cdef int orientation
    cdef string name

    def __cinit__(PyGenome self, size_t size=1000003):
        self._genome = Genome(size)
        self.orientation = 0
        self.name = ""

    cpdef double load_factor(PyGenome self):
        return self._genome.load_factor()

    def __len__(PyGenome self):
        return self._genome.size()

    def __repr__(PyGenome self):
        return "<Genome: SNPs=%d, name=%r>" % (self.__len__(), self.name)

    def __getitem__(PyGenome self, key):
        """Retrieves genotype keyed by its RSID.

        Arguments:
            key: An RSID as a string or integer.

        Raises:
            KeyError - SNP not found.

        Returns:
            A one or two-character string containing the genotype for this
            RSID. The order of the two characters is the same as given in the
            source file.

            The characters that are used are `A`, `T`, `C`, `G` as well as `D`,
            `I` and `-`. The first four are nucleotides, while `D` indicates
            deletion, `I` insertion. SNPs that are present in the genome, but
            with no result, are always indicated by two dashes, `--`. Some SNPs
            only have one nucleotide, for example from the Y chromosome or
            mitochondrial DNA.

        Usage:
            >>> genome["rs2534636"]
            'C'
            >>> genome["rs123"]
            'AA'
            >>> genome["rs28504042"]
            '--'
            >>> genome["rs3135027"]
            'G'
        """
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
    def y_chromosome(PyGenome self):
        """Flag indicating presence of a Y chromosome."""
        return self._genome.y_chromosome

def load(string filename, name=None, size_t initial_size=1000003):
    """Loads given 23andMe raw genome file.

    Arguments:
        filename: Name of file to load.
        name (optional): Name to give the genome
        initial_size (optional): Number of initial empty slots to reserve in
            the underlying hash table.

    Raises:
        RuntimeError

    Returns:
        A ``Genome``.
    """
    genome = PyGenome(initial_size)
    parse_file(filename, genome._genome)
    genome.name = name if name is not None else filename
    return genome

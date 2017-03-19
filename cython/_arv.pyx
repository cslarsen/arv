# cython: c_string_type=unicode, c_string_encoding=utf8

# arv
# Copyright 2017 Christian Stigen Larsen
# Distributed under the GNU GPL v3 or later; see COPYING.

from libc.stdint cimport uint32_t, int32_t
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "arv.hpp" namespace "arv":
    ctypedef uint32_t Position
    ctypedef int32_t RSID

    cdef enum Nucleotide:
        NONE, A, G, C, T, D, I

    cdef enum Chromosome:
        CHR_NO =  0,
        CHR_01 =  1,
        CHR_02 =  2,
        CHR_03 =  3,
        CHR_04 =  4,
        CHR_05 =  5,
        CHR_06 =  6,
        CHR_07 =  7,
        CHR_08 =  8,
        CHR_09 =  9,
        CHR_10 = 10,
        CHR_11 = 11,
        CHR_12 = 12,
        CHR_13 = 13,
        CHR_14 = 14,
        CHR_15 = 15,
        CHR_16 = 16,
        CHR_17 = 17,
        CHR_18 = 18,
        CHR_19 = 19,
        CHR_20 = 20,
        CHR_21 = 21,
        CHR_22 = 22,
        CHR_X  = 23,
        CHR_Y  = 24,
        CHR_MT = 25

    cdef cppclass Genotype:
        Nucleotide first
        Nucleotide second

        Genotype()
        Genotype(const Nucleotide&, const Nucleotide&)
        bool operator==(const Genotype&) const
        bool operator<(const Genotype&) const

        string to_string() const

    cdef cppclass SNP:
        Chromosome chromosome
        Position position
        Genotype genotype

        SNP(const Chromosome& = CHR_NO,
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
    cdef Genotype complement(const Genotype&)

cdef class PyGenotype:
    """A pair of nucleotides."""
    cdef Genotype _genotype

    def __repr__(self):
        return "'%s'" % self._genotype.to_string()

    def __str__(self):
        return str(self._genotype.to_string())

    def __richcmp__(self, obj, int op):
        # The PyGenotype may look like a string, so allow comparisons with
        # strings
        this = str(self)
        that = str(obj)

        if op == 0:
            return this < that
        elif op == 1:
            return this <= that
        elif op == 2:
            return this == that
        elif op == 3:
            return this != that
        elif op == 4:
            return this > that
        elif op == 5:
            return this >= that
        else:
            raise NotImplementedError()

    def __invert__(self):
        gt = PyGenotype()
        gt._genotype = complement(self._genotype)
        return gt

cdef class PySNP:
    """A single nucleotide polymorphism."""
    cdef SNP _snp

    def __cinit__(self):
        self._snp = NONE_SNP

    @property
    def position(self):
        return self._snp.position;

    @property
    def chromosome(self):
        """Returns the chromosome.

        Possible values are the integers 1 through 22, and the strings "MT",
        "X" and "Y".
        """
        cdef Chromosome c = self._snp.chromosome

        if CHR_01 <= c <= CHR_22:
            return c
        else:
            return {CHR_NO: None, CHR_MT: "MT", CHR_X: "X", CHR_Y: "Y"}[c]

    @property
    def genotype(self):
        """Returns the Genotype."""
        gt = PyGenotype()
        gt._genotype = self._snp.genotype
        return gt

    def __repr__(self):
        return "<SNP: chromosome=%r position=%r genotype=%r>" % (
                self.chromosome, self.position, self.genotype)

    def __str__(self):
        return str(self.genotype)

cdef class PyGenome:
    """A collection of SNPs for a human being.

    Provides a dictionary interface, including keys, values, __get_item__ and
    so on.
    """
    cdef Genome _genome
    cdef int orientation
    cdef str name
    cdef str ethnicity

    def __cinit__(PyGenome self, size_t size=1000003):
        self._genome = Genome(size)
        self.orientation = 0
        self.name = ""
        self.ethnicity = ""

    cpdef double load_factor(PyGenome self):
        """The underlying hash table's load factor."""
        return self._genome.load_factor()

    def keys(self):
        return self._genome.rsids()

    cdef vector[SNP] values(self):
        return self._genome.snps()

    def get_snp(self, key):
        """Retrieves given SNP.

        Arguments:
            key: An RSID as a string ("rs123") or an integer (123, for
                 example). Internal IDs ("i123") are represented as negative
                 numbers (-123).

        Returns:
            A ``SNP``.

        Raises:
            KeyError - if RSID was not found.
        """

        cdef RSID rsid

        if isinstance(key, str):
            if key.startswith("rs"):
                rsid = int(key[2:])
            elif key.startswith("i"):
                rsid = -int(key[1:])
            else:
                raise KeyError(key)
        elif isinstance(key, int):
            rsid = key
        else:
            raise TypeError("Expected str or int but got %s" %
                    type(key).__name__)
        cdef SNP s = self._genome[rsid]
        if s == NONE_SNP:
            raise KeyError(key)
        else:
            snp = PySNP()
            snp._snp = s
            return snp

    def __len__(self):
        return self._genome.size()

    def __repr__(self):
        return "<Genome: SNPs=%d, name=%r>" % (self.__len__(), self.name)

    def __getitem__(self, key):
        """Retrieves genotype keyed by its RSID.

        Arguments:
            key: An RSID as a string (e.g. "rs123" or "i123") or a plain
            integer. Positive integers are used for normal RSIDs (e.g. 123 for
            "rs123") and negative integers for internal IDs (e.g. -123 for
            "i123"). The reason we support integers is because ``keys()`` will
            return integers, because it is so much faster to work with and
            takes up much less space.

        Raises:
            KeyError - SNP not found.

        Returns:
            A one or two-character string containing the genotype for this
            RSID. The order of the two characters is the same as given in the
            source file.

            The characters that are used are ``A``, ``T``, ``C``, ``G`` as well
            as ``D``, ``I`` and ``-``. The first four are nucleotides, while
            ``D`` indicates deletion, ``I`` insertion. SNPs that are present in
            the genome, but with no result, are always indicated by two dashes,
            ``--``. Some SNPs only have one nucleotide, for example from the Y
            chromosome or mitochondrial DNA.

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
            if key.startswith("rs"):
                rsid = int(key[2:])
            elif key.startswith("i"):
                rsid = -int(key[1:])
            else:
                raise KeyError(key)
        elif isinstance(key, int):
            rsid = key
        else:
            raise KeyError(key)

        cdef SNP s = self._genome[rsid]
        if s == NONE_SNP:
            raise KeyError(key)
        else:
            return str(s.genotype.to_string())

    @property
    def y_chromosome(PyGenome self):
        """Flag indicating presence of a Y chromosome."""
        return self._genome.y_chromosome

def load(filename, name=None, ethnicity=None, size_t initial_size=1000003):
    """Loads given 23andMe raw genome file.

    Arguments:
        filename:        Name of file to load.
        name (optional): Name to give the genome. Will use filename by default.
        ethnicity:       Provide genome's ethnicity to unlock more reports.
                         Accepted values are typically "european", "african",
                         "asian". None by default.
        initial_size (optional): Number of initial empty slots to reserve in
                                 the underlying hash table.

    Raises:
        RuntimeError - instead of FileNotFoundError etc.

    Returns:
        A ``Genome``.
    """
    genome = PyGenome(initial_size)
    parse_file(filename.encode("utf-8"), genome._genome)
    genome.name = name if name is not None else filename
    genome.ethnicity = ethnicity if ethnicity is not None else ""
    return genome

def _sizes():
    """Returns C++ sizeof() for internal structures."""
    return {
        "Chromosome":     sizeof(Chromosome),
        "Genome":         sizeof(Genome),
        "GenomeIterator": sizeof(GenomeIterator),
        "Genotype":       sizeof(Genotype),
        "Nucleotide":     sizeof(Nucleotide),
        "Position":       sizeof(Position),
        "RSID":           sizeof(RSID),
        "RsidSNP":        sizeof(RsidSNP),
        "SNP":            sizeof(SNP),
    }

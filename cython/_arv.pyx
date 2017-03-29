# cython: c_string_type=unicode, c_string_encoding=utf8

# arv
# Copyright 2017 Christian Stigen Larsen
# Distributed under the GNU GPL v3 or later; see COPYING.

from libc.stdint cimport uint32_t, int32_t
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.utility cimport pair
from libcpp.vector cimport vector

cdef extern from "arv.hpp" namespace "arv":
    ctypedef uint32_t Position
    ctypedef int32_t RSID

    cdef enum Nucleotide:
        NONE, A, G, C, T, D, I

    cdef enum Chromosome:
        CHR_NO =  0, CHR_01 =  1,
        CHR_02 =  2, CHR_03 =  3,
        CHR_04 =  4, CHR_05 =  5,
        CHR_06 =  6, CHR_07 =  7,
        CHR_08 =  8, CHR_09 =  9,
        CHR_10 = 10, CHR_11 = 11,
        CHR_12 = 12, CHR_13 = 13,
        CHR_14 = 14, CHR_15 = 15,
        CHR_16 = 16, CHR_17 = 17,
        CHR_18 = 18, CHR_19 = 19,
        CHR_20 = 20, CHR_21 = 21,
        CHR_22 = 22, CHR_X  = 23,
        CHR_Y  = 24, CHR_MT = 25

    cdef cppclass CGenotype "arv::Genotype":
        Nucleotide first
        Nucleotide second

        CGenotype()
        CGenotype(const Nucleotide&, const Nucleotide&)
        bool operator==(const CGenotype&) const
        bool operator<(const CGenotype&) const
        string to_string() const

    cdef cppclass CSNP "arv::SNP":
        Chromosome chromosome
        Position position
        CGenotype genotype

        CSNP(const Chromosome& = CHR_NO,
             const Position& = 0,
             const CGenotype& = CGenotype(NONE, NONE))
        CSNP(const CSNP&)
        CSNP& operator=(const CSNP&)
        bool operator!=(const CSNP&) const
        bool operator<(const CSNP&) const
        bool operator<=(const CSNP&) const
        bool operator==(const CGenotype&) const
        bool operator==(const CSNP&) const
        bool operator>(const CSNP&) const
        bool operator>=(const CSNP&) const

    ctypedef pair[RSID, CSNP] RsidSNP

    cdef cppclass CGenomeIterator "arv::GenomeIterator":
        CGenomeIterator()
        CGenomeIterator(const CGenomeIterator&)
        CGenomeIterator& operator=(const CGenomeIterator&)
        bool operator==(const CGenomeIterator&) const
        bool operator!=(const CGenomeIterator&) const
        void next()
        RsidSNP value() const

    cdef const CSNP NONE_SNP

    cdef cppclass CGenome "arv::Genome":
        CGenome() except +
        CGenome(const size_t) except +
        CGenome(const CGenome&) except +
        CGenome& operator=(const CGenome&)

        bool y_chromosome

        double load_factor() const
        size_t size() const

        string genotype(const RSID& id) const
        const CSNP& operator[](const RSID&) const
        bool has(const RSID&) const
        void insert(const RSID&, const CSNP&)

        bool operator==(const CGenome&) const
        bool operator!=(const CGenome&) const

        CGenomeIterator begin() const
        CGenomeIterator end() const

    cdef void parse_file(const string&, CGenome&) except +
    cdef CGenotype complement(const CGenotype&)

cdef basestring __rsid2str(const RSID& rsid):
    """Converts RSID to integer."""
    if rsid >= 0:
        return "rs%d" % rsid
    else:
        return "i%d" % -rsid

cdef RSID __rsid2int(key) except? 0:
    """Converts RSID to string."""
    if isinstance(key, int):
        return key
    elif isinstance(key, str):
        if key.startswith("rs"):
            return int(key[2:])
        elif key.startswith("i"):
            return -int(key[1:])
    else:
        return 0

cdef class Genotype(object):
    """Contains a pair of nucleotides.

    The characters that are used are ``A``, ``T``, ``C``, ``G``, in addition to
    ``D`` (deletion), ``I`` (insertion) and ``-`` (no call). The last dash
    always occurs in pairs, as in ``--``, indicating that 23andMe wasn't able
    to call it.

    Some chromosomes naturally have only one nucleotide, for example the
    Y-chromosome and mitochondrial DNA (``MT``).

    You can convert a ``Genotype`` to a string with ``str()``. Since this is
    such a common operation, you can compare directly with strings. For
    example, this is valid:

        >>> gt
        <Genotype 'AT'>
        >>> gt == "AT"
        True
        >>> gt == "AA"
        False
    """
    cdef CGenotype _genotype

    def __cinit__(self):
        self._genotype = CGenotype(NONE, NONE)

    @staticmethod
    cdef _init(CGenotype genotype):
        r = Genotype()
        r._genotype = genotype
        return r

    def __repr__(self):
        return "<Genotype %r>" % str(self)

    def __str__(self):
        return str(self._genotype.to_string())

    cpdef _rich_cmp(self, Genotype obj, int op):
        cdef CGenotype this = self._genotype
        cdef CGenotype that = obj._genotype

        if op == 0:
            return this < that
        elif op == 1:
            return (this == that) or (this < that)
        elif op == 2:
            return this == that
        elif op == 3:
            return not (this == that)
        elif op == 4:
            return that < this
        elif op == 5:
            return (this == that) or (that < this)
        raise NotImplementedError()

    def __richcmp__(self, obj, int op):
        if isinstance(obj, str):
            # String comparison
            this = str(self)
            that = str(obj)
        elif isinstance(obj, Genotype):
            # Genotype comparison
            return self._rich_cmp(obj, op)
        else:
            raise NotImplementedError()

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

        raise NotImplementedError()

    def __invert__(self):
        gt = Genotype()
        gt._genotype = complement(self._genotype)
        return gt

cdef class SNP(object):
    """A single nucleotide polymorphism.

    Contains a ``genotype``, which ``chromosome`` it belongs to, and its
    ``position`` on the reference human genome.

    Comparisons
    -----------

    ``SNP``\ s support rich comparisons operator. Comparing two SNPs is
    equivalent to comparing two tuples of (position, chromosome, genotype).

    You may also compare a ``SNP`` with a string. In that case, it will perform
    a string comparison with its genotype. For example

        >>> snp.genotype
        'AT'
        >>> snp == "AT"
        True
    """
    cdef CSNP _snp

    def __cinit__(self):
        self._snp = NONE_SNP

    @staticmethod
    cdef _init(CSNP snp):
        r = SNP()
        r._snp = snp
        return r

    @property
    def position(self):
        """Location on the reference human genome."""
        return self._snp.position

    @property
    def chromosome(self):
        """Which ``Chromosome`` the SNP belongs to.

        Possible values are the integers 1 through 22 and the strings "MT"
        (mitochondrial DNA) "X" and "Y".
        """
        cdef Chromosome c = self._snp.chromosome

        if CHR_01 <= c <= CHR_22:
            return c
        else:
            return {CHR_NO: None, CHR_MT: "MT", CHR_X: "X", CHR_Y: "Y"}[c]

    @property
    def genotype(self):
        """The SNPs called ``Genotype``."""
        return Genotype._init(self._snp.genotype)

    def __repr__(self):
        return "<SNP: chromosome=%r position=%r genotype=%r>" % (
                self.chromosome, self.position, self.genotype)

    def __str__(self):
        return str(self.genotype)

    cpdef _rich_cmp(self, SNP obj, int op):
        cdef CSNP this = self._snp
        cdef CSNP that = obj._snp

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
        raise NotImplementedError()

    def __richcmp__(self, obj, int op):
        if isinstance(obj, str):
            # String comparison
            this = str(self.genotype)
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
        elif isinstance(obj, SNP):
            # Genotype comparison
            return self._rich_cmp(obj, op)
        raise NotImplementedError()


cdef class GenomeIterator(object):
    """Iterates through the ``SNP`` objects of a ``Genome``."""
    cdef CGenomeIterator _cur
    cdef CGenomeIterator _end
    cdef int _type

    def __cinit__(GenomeIterator self, int iterator_type):
        self._type = iterator_type

    @staticmethod
    cdef GenomeIterator _iterate(const CGenome& genome, const int
            iterator_type):
        it = GenomeIterator(iterator_type)
        it._cur = genome.begin()
        it._end = genome.end()
        return it

    def __iter__(self):
        return self

    def next(self):
        # Python 2 compatibility
        return self.__next__()

    def __next__(self):
        return self._dispatch()

    cdef _dispatch(GenomeIterator self):
        if self._type == 0:
            return self._next_key()
        elif self._type == 1:
            return self._next_value()
        else:
            return self._next_item()

    cdef SNP _next_value(GenomeIterator self):
        if self._cur == self._end:
            raise StopIteration()

        snp = SNP._init(self._cur.value().second)
        self._cur.next()
        return snp

    cdef basestring _next_key(GenomeIterator self):
        if self._cur == self._end:
            raise StopIteration()

        cdef RSID rsid = self._cur.value().first
        self._cur.next()
        return __rsid2str(rsid)

    cdef _next_item(GenomeIterator self):
        if self._cur == self._end:
            raise StopIteration()

        cdef basestring rsid = __rsid2str(self._cur.value().first)
        cdef SNP snp = SNP._init(self._cur.value().second)
        self._cur.next()
        return (rsid, snp)


cdef class Genome(object):
    """A collection of SNPs for a human being.

    Implements the dictionary protocol.
    """
    cdef CGenome _genome
    cdef int _orientation
    cdef basestring _name
    cdef basestring _ethnicity

    def __cinit__(Genome self, size_t size=1000003):
        self._genome = CGenome(size)
        self._orientation = 0
        self._name = ""
        self._ethnicity = ""

    cpdef double load_factor(Genome self):
        """The underlying hash table's load factor."""
        return self._genome.load_factor()

    @property
    def orientation(self):
        """The orientation of the genotype call, represented as an integer
        being either +1 (plus strand) or -1 (minus strand).

        For 23andMe genomes, this will always be +1.
        """
        return self._orientation

    @orientation.setter
    def orientation(self, int value):
        if value not in [-1, 1]:
            raise ValueError("Orientation must be either -1 or +1")
        self._orientation = value

    @property
    def name(self):
        """A user-defined name that goes along with the genome.

        If not explicitly set, it will be set to the filename.
        """
        return self._name

    @name.setter
    def name(self, basestring value):
        self._name = value

    @property
    def ethnicity(self):
        """A user-defined ethnicity that goes along with the genome.

        Certain reports are only valid for given ethnic groups, and therefore
        this serves as a way to detect which reports are applicable to this
        genome.

        Usual values are ``european``, ``asian``, ``african``.
        """
        return self._ethnicity

    @ethnicity.setter
    def ethnicity(self, basestring value):
        self._ethnicity = value

    @property
    def y_chromosome(Genome self):
        """A boolean indicating the presence of a Y-chromosome within this
        genome."""
        return self._genome.y_chromosome

    def keys(self):
        return GenomeIterator._iterate(self._genome, 0)

    def values(self):
        return GenomeIterator._iterate(self._genome, 1)

    def items(self):
        return GenomeIterator._iterate(self._genome, 2)

    def __len__(self):
        return self._genome.size()

    def __repr__(self):
        return "<Genome: SNPs=%d, orientation=%d, name=%r, ethnicity=%r>" % (
                self.__len__(), self.orientation, self.name, self.ethnicity)

    def __contains__(self, key):
        return self._genome.has(__rsid2int(key))

    def __getitem__(self, key):
        """Retrieves SNP from its RSID.

        Arguments:
            key: An RSID as a string (e.g. "rs123" or "i123") or a plain
                integer. Positive integers are used for normal RSIDs (e.g. 123
                for "rs123") and negative integers for internal IDs (e.g. -123 for
                "i123").

        Raises:
            KeyError - RSID not found in genome.

        Returns:
            An ``SNP``.

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
        cdef RSID rsid = __rsid2int(key)
        cdef CSNP snp = self._genome[rsid]

        if snp == NONE_SNP:
            raise KeyError(key)

        return SNP._init(snp)


def load(filename, name=None, ethnicity=None, size_t initial_size=1000003,
        orientation=1):
    """Loads given 23andMe raw genome file.

    Arguments:
        filename:        Name of file to load.

        name (optional): Name to give the genome. Will use filename by default.

        ethnicity (optional): Provide genome's ethnicity to unlock more reports.
                              Accepted values are typically "european", "african",
                              "asian". None by default.

        initial_size (optional): Number of initial empty slots to reserve in
                                 the underlying hash table.

        orientation (optional): +1 or -1, corresponding to the plus or minus
                                strand. 23andMe files have always plus
                                orientation.

    Raises:
        RuntimeError - File not found, parser errors, etc.

    Returns:
        A ``Genome``.
    """
    cdef Genome genome = Genome(0)
    parse_file(filename.encode("utf-8"), genome._genome)
    genome.name = name if name is not None else filename
    genome.ethnicity = ethnicity if ethnicity is not None else ""
    genome.orientation = orientation
    return genome

def _sizes():
    """Returns C++ sizeof() for internal structures."""
    return {
        "CGenome":          sizeof(CGenome),
        "CGenomeIterator":  sizeof(CGenomeIterator),
        "CGenotype":        sizeof(CGenotype),
        "Chromosome":       sizeof(Chromosome),
        "CSNP":             sizeof(CSNP),
        "Nucleotide":       sizeof(Nucleotide),
        "Position":         sizeof(Position),
        "RSID":             sizeof(RSID),
        "RsidSNP":          sizeof(RsidSNP),
    }

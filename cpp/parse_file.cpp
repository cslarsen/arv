/*
 * dna-traits
 * Copyright 2014, 2016 Christian Stigen Larsen
 * Distributed under the GNU GPL v3 or later. See COPYING.
 *
 * arv
 * Copyright 2017 Christian Stigen Larsen
 * Distributed under the GNU GPL V3 or later. See COPYING.
 */

#include "dnatraits.hpp"
#include "file.hpp"
#include "filesize.hpp"
#include "mmap.hpp"

static Nucleotide CharToNucleotide[256] = {NONE};

static inline void skip_comments(const char*& s)
{
  while ( *s == '#' )
    while ( *s++ != '\n' )
      ; // loop
}

static inline bool iswhite(const char c)
{
  return c=='\t' || c=='\n' || c=='\r';
}

static inline const char*& skipwhite(const char*& s)
{
  while ( iswhite(*s) ) ++s;
  return s;
}

static inline uint32_t parse_uint32(const char*& s)
{
  uint32_t n = 0;

  while ( isdigit(*s) )
    n = n*10 - '0' + *s++;

  return n;
}

static inline int32_t parse_int32(const char*& s)
{
  int32_t n = 0;

  while ( isdigit(*s) )
    n = n*10 - '0' + *s++;

  return n;
}

static inline Nucleotide parse_nucleotide(const char*& s)
{
  return CharToNucleotide[static_cast<const std::size_t>(*s++)];
}

static inline Chromosome parse_chromo(const char*& s)
{
  if ( isdigit(*s) )
      return static_cast<Chromosome>(parse_uint32(s));

  switch ( *s++ ) {
    case 'M': ++s; // skip T in "MT"
              return CHR_MT;
    case 'X': return CHR_X;
    case 'Y': return CHR_Y;
    default:  return NO_CHR;
  }
}

static inline Genotype parse_genotype(const char*& s)
{
  Nucleotide first = parse_nucleotide(s);
  Nucleotide second = parse_nucleotide(s);
  return Genotype(first, second);
}

static inline void skipline(const char*& s)
{
  while ( *s != '\n' ) ++s;
}

/**
 * Reads a 23andMe-formatted genome file.  It currently uses reference human
 * assembly build 37 (annotation release 104).
 */
void parse_file(const std::string& name, Genome& genome)
{
  CharToNucleotide[static_cast<unsigned>('A')] = A;
  CharToNucleotide[static_cast<unsigned>('G')] = G;
  CharToNucleotide[static_cast<unsigned>('C')] = C;
  CharToNucleotide[static_cast<unsigned>('T')] = T;
  CharToNucleotide[static_cast<unsigned>('D')] = D;
  CharToNucleotide[static_cast<unsigned>('I')] = I;

  File fd(name.c_str(), O_RDONLY);
  MMap fmap(0, filesize(fd), PROT_READ, MAP_PRIVATE, fd, 0);
  auto s = static_cast<const char*>(fmap.ptr());

  skip_comments(s);
  bool ychromo = false;

  // Local cache of SNPs and RSIDs, for more locality and hence more speed. Its
  // size is somewhat arbitrary, but shouldn't be too big.
  #define SIZE 200
  SNP snps[SIZE];
  RSID rsids[SIZE];
  int i=0;
  bool internal=false;

  for ( ; *s; ++s ) {
    if (*s == 'i')
      internal = true;
    else if ( *s == 'r')
      internal = false;
    else {
      skipline(s);
      continue;
    }

    RSID& rsid = rsids[i];
    SNP& snp = snps[i];

    rsid = parse_int32(s+=internal? 1 : 2); // skip "i"/"rs"-prefix
    if (internal)
      rsid = -rsid;

    if ( rsid < genome.first ) genome.first = rsid;
    if ( rsid > genome.last ) genome.last = rsid;

    snp.chromosome = parse_chromo(skipwhite(s));
    snp.position = parse_uint32(skipwhite(s));
    snp.genotype = parse_genotype(skipwhite(s));

    ychromo |= (snp.chromosome==CHR_Y && snp.genotype.first!=NONE);

    // Ordinarly, we would just call `genome.insert(rsid, snp)` here, but it's
    // a tad faster to stage them in an array first, and then flush it to the
    // hash map when it's full.

    if ( ++i == SIZE ) {
      i = 0;
      for ( int n = 0; n < SIZE; ++n )
        genome.insert(rsids[n], snps[n]);
    }
  }

  // flush the rest
  for ( int n=0; n < i; ++n )
    genome.insert(rsids[n], snps[n]);

  genome.y_chromosome = ychromo;
}

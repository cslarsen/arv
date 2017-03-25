/*
 * dna-traits
 * Copyright 2014, 2016 Christian Stigen Larsen
 * Distributed under the GNU GPL v3 or later. See COPYING.
 *
 * arv
 * Copyright 2017 Christian Stigen Larsen
 * Distributed under the GNU GPL V3 or later. See COPYING.
 */

#include "arv.hpp"
#include "file.hpp"
#include "filesize.hpp"
#include "mmap.hpp"

#ifdef __GNUC__
#define likely(x)   __builtin_expect((x),1)
#define unlikely(x) __builtin_expect((x),0)
#else
#define likely(x)   (x)
#define unlikely(x) (x)
#endif

namespace arv {

static Nucleotide CharToNucleotide[256] = {NONE};

static void skip_comments(const char*& s)
{
  while ( *s == '#' )
    while ( *s++ != '\n' )
      ; // loop
}

static bool iswhite(const char c)
{
  return c=='\t' || c=='\n' || c=='\r';
}

static const char*& skipwhite(const char*& s)
{
  while ( iswhite(*s) ) ++s;
  return s;
}

static uint32_t parse_uint32(const char*& s)
{
  uint32_t n = 0;

  while ( isdigit(*s) )
    n = n*10 - '0' + *s++;

  return n;
}

static int32_t parse_int32(const char*& s)
{
  int32_t n = 0;

  while ( isdigit(*s) )
    n = n*10 - '0' + *s++;

  return n;
}

static Nucleotide parse_nucleotide(const char*& s)
{
  return CharToNucleotide[static_cast<short>(*s++)];
}

static Chromosome parse_chromo(const char*& s)
{
  if ( likely(isdigit(*s)) )
      return static_cast<Chromosome>(parse_uint32(s));

  const char c = *s++;

  if ( c == 'X' )
    return CHR_X;

  if ( c == 'M' ) {
    ++s; // skip T in "MT"
    return CHR_MT;
  }

  if ( c == 'Y' )
    return CHR_Y;

  return CHR_NO;
}

static Genotype parse_genotype(const char*& s)
{
  const Nucleotide first = parse_nucleotide(s);
  const Nucleotide second = parse_nucleotide(s);
  return Genotype(first, second);
}

static void skipline(const char*& s)
{
  while ( *s != '\n' ) ++s;
}

/**
 * Reads a 23andMe-formatted genome file.  It currently uses reference human
 * assembly build 37 (annotation release 104).
 */
void parse_file(const std::string& name, Genome& genome)
{
  using namespace arv;

  CharToNucleotide[static_cast<short>('-')] = NONE;
  CharToNucleotide[static_cast<short>('A')] = A;
  CharToNucleotide[static_cast<short>('C')] = C;
  CharToNucleotide[static_cast<short>('D')] = D;
  CharToNucleotide[static_cast<short>('G')] = G;
  CharToNucleotide[static_cast<short>('I')] = I;
  CharToNucleotide[static_cast<short>('T')] = T;

  File fd(name.c_str(), O_RDONLY);
  MMap fmap(0, filesize(fd), PROT_READ, MAP_PRIVATE, fd, 0);
  auto s = fmap.c_str();

  skip_comments(s);

  // Local cache of SNPs and RSIDs, for more locality and hence more speed. Its
  // size is somewhat arbitrary, but shouldn't be too big.
  const std::size_t BUFFER_SIZE = 200;
  RsidSNP buffer[BUFFER_SIZE];
  size_t buffer_pos = 0;

  bool internal = false; // rsid or internal id

  for ( ; *s; ++s ) {
    if ( *s == 'r' )
      internal = false;
    else if ( *s == 'i' )
      internal = true;
    else {
      skipline(s);
      continue;
    }

    RSID& rsid = buffer[buffer_pos].first;
    SNP& snp = buffer[buffer_pos].second;

    // Skip i/rs prefix and parse number
    if ( !internal )
      rsid = parse_int32(s += 2);
    else
      rsid = -parse_int32(s += 1);

    snp.chromosome = parse_chromo(skipwhite(s));
    snp.position = parse_uint32(skipwhite(s));
    snp.genotype = parse_genotype(skipwhite(s));

    genome.y_chromosome |= snp.chromosome == CHR_Y;

    // Ordinarly, we would just call `genome.insert(rsid, snp)` here, but it's
    // a tad faster to stage them in an array first, and then flush it to the
    // hash map when it's full.

    if ( ++buffer_pos == BUFFER_SIZE ) {
      buffer_pos = 0;
      for ( size_t n = 0; n < BUFFER_SIZE; ++n )
        genome.insert(buffer[n]);
    }
  }

  // Store the rest of the buffer
  for ( size_t n = 0; n < buffer_pos; ++n )
    genome.insert(buffer[n]);
}

} // namespace arv

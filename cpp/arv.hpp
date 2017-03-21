/*
 * dna-traits
 * Copyright 2014, 2016 Christian Stigen Larsen
 * Distributed under the GNU GPL v3 or later. See COPYING.
 *
 * arv
 * Copyright 2017 Christian Stigen Larsen
 * Distributed under the GNU GPL v3 or later. See COPYING.
 */

#ifndef ARV_ARV_HPP
#define ARV_ARV_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "export.hpp"

namespace arv {

typedef std::uint32_t Position;
typedef std::int32_t RSID;

enum Nucleotide {
  NONE, A, G, C, T, D, I
};

enum Chromosome {
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
  CHR_MT = 25 // Mitochondrial DNA
};

// We can get this down to a byte if we want to
#pragma pack(1)
struct DLL_PUBLIC Genotype {
  Nucleotide first  : 3;
  Nucleotide second : 3;

  Genotype();
  Genotype(const Nucleotide& a, const Nucleotide& b);

  friend Genotype operator~(const Genotype&);
  bool operator==(const Genotype& g) const;
  bool operator<(const Genotype& g) const;

  std::string to_string() const;
};

#pragma pack(1)
struct DLL_PUBLIC SNP {
  Chromosome chromosome : 5;
  Position position;
  Genotype genotype;

  SNP();
  SNP(const Chromosome&, const Position&, const Genotype&);
  SNP(const SNP&);

  SNP& operator=(const SNP&);
  bool operator!=(const SNP&) const;
  bool operator<(const SNP&) const;
  bool operator<=(const SNP&) const;
  bool operator==(const Genotype&) const;
  bool operator==(const SNP&) const;
  bool operator>(const SNP&) const;
  bool operator>=(const SNP&) const;
};

extern DLL_PUBLIC const SNP NONE_SNP;

struct DLL_LOCAL GenomeIteratorImpl;

struct DLL_PUBLIC RsidSNP {
  RSID rsid;
  SNP snp;

  bool operator==(const RsidSNP& o) const;
};

struct DLL_PUBLIC GenomeIterator {
  // todo copy ctor, assignment op, dtor
  GenomeIterator(GenomeIteratorImpl*);
  ~GenomeIterator();
  GenomeIterator(const GenomeIterator&);
  GenomeIterator& operator=(const GenomeIterator&);
  GenomeIterator& operator++();
  bool operator==(const GenomeIterator&);
  bool operator!=(const GenomeIterator&);
  const RsidSNP operator*();
private:
  GenomeIteratorImpl* DLL_LOCAL pimpl;
};

struct DLL_PUBLIC Genome {
  /*!
   * True if genome contains a Y-chromosome (with non-empty genotypes).
   */
  bool y_chromosome;

  /*!
   * Lowest RSID.
   */
  RSID first;

  /*!
   * Highest RSID.
   */
  RSID last;

  Genome();
  Genome(const std::size_t size);
  Genome(const Genome&);
  Genome& operator=(const Genome&);
  ~Genome();

  /*!
   * Access SNP. Throws on not found.
   */
  const SNP& operator[](const RSID& id) const;

  /*!
   * Checks if hash table contains given RSID.
   */
  bool has(const RSID& id) const;

  /*!
   * Add a SNP to the hash table.
   */
  void insert(const RSID& rsid, const SNP& snp);
  void insert(const std::pair<RSID, SNP>&);

  /*!
   * Underlying hash table's load factor. (For developer purposes)
   */
  double load_factor() const;

  /*!
   * Number of SNPs.
   */
  std::size_t size() const;

  /*!
   * Returns RSIDs that exist in both genomes.
   */
  std::vector<RSID> intersect_rsid(const Genome& genome) const;

  /*!
   * Returns RSID for SNPs that have the same genotype in both genomes.
   */
  std::vector<RSID> intersect_snp(const Genome& genome) const;

  /*!
   * Returns all RSIDs in this genome.
   */
  std::vector<RSID> rsids() const;

  /*!
   * Returns a copy of all SNPs.
   */
  std::vector<SNP> snps() const;

  bool operator==(const Genome&) const;
  bool operator!=(const Genome&) const;

  GenomeIterator begin() const;
  GenomeIterator end() const;

private:
  struct DLL_LOCAL GenomeImpl;
  GenomeImpl* DLL_LOCAL pimpl;
};

Nucleotide complement(const Nucleotide& n);

/*!
 * Parse a 23andMe genome text file and put contents into genome.
 */
void DLL_PUBLIC parse_file(const std::string& filename, Genome&);

Genotype DLL_PUBLIC complement(const Genotype& g);

} // namespace arv

#endif // include guard

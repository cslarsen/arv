/*
 * dna-traits
 * Copyright 2014, 2016 Christian Stigen Larsen
 * Distributed under the GNU GPL v3 or later. See COPYING.
 *
 * arv
 * Copyright 2017 Christian Stigen Larsen
 * Distributed under the GNU GPL v3 or later. See COPYING.
 */

#include "mmap.hpp"

#include <stdexcept>

namespace arv {

DLL_LOCAL MMap::MMap(void *address,
           const std::size_t length,
           const int protection_level,
           const int flags,
           const int file_descriptor,
           const off_t offset)
  : l(length),
    p(mmap(address, length, protection_level, flags, file_descriptor, offset))
{
  if ( p == reinterpret_cast<caddr_t>(-1) )
    throw std::runtime_error("mmap error");
}

DLL_LOCAL MMap::~MMap() {
  munmap(p, l);
}

} // ns arv

/*
 * dna-traits
 * Copyright 2014, 2016 Christian Stigen Larsen
 * Distributed under the GNU GPL v3 or later. See COPYING.
 *
 * arv
 * Copyright 2017 Christian Stigen Larsen
 * Distributed under the GNU GPL v3 or later. See COPYING.
 */

#ifndef ARV_MMAP_HPP
#define ARV_MMAP_HPP

#include <cstddef>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/types.h>

#ifndef BUILDING_DLL
# define BUILDING_DLL
#endif
#include "export.hpp"

namespace arv {

class DLL_LOCAL MMap {
  std::size_t l;
  void *p;
public:
  MMap(void *address,
       const std::size_t length,
       const int protection_level,
       const int flags,
       const int file_descriptor,
       const off_t offset);
  ~MMap();

  inline void* ptr() const {
    return p;
  }

  inline const char* c_str() const {
    return static_cast<const char*>(p);
  }
};

} // ns arv

#endif // guard

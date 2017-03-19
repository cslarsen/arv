/*
 * dna-traits
 * Copyright 2014, 2016 Christian Stigen Larsen
 * Distributed under the GNU GPL v3 or later. See COPYING.
 *
 * arv
 * Copyright 2017 Christian Stigen Larsen
 * Distributed under the GNU GPL v3 or later. See COPYING.
 */

#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/types.h>

#define BUILDING_DLL
#include "export.hpp"

class DLL_LOCAL MMap {
  size_t l;
  void *p;
public:
  MMap(void *address,
       size_t length,
       int protection_level,
       int flags,
       int file_descriptor,
       off_t offset);

  ~MMap();

  inline void* ptr() const {
    return p;
  }
};

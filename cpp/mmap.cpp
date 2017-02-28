/*
 * Copyright (C) 2014, 2016 Christian Stigen Larsen
 * Distributed under the GPL v3 or later. See COPYING.
 */

#include <stdexcept>
#include "mmap.hpp"

DLL_LOCAL MMap::MMap(void *address,
           size_t length,
           int protection_level,
           int flags,
           int file_descriptor,
           off_t offset)
  : l(length),
    p(mmap(address,
           length,
           protection_level,
           flags,
           file_descriptor,
           offset))
{
  if ( p == reinterpret_cast<caddr_t>(-1) )
    throw std::runtime_error("mmap error");
}

DLL_LOCAL MMap::~MMap() {
  munmap(p, l);
}

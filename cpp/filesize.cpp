/*
 * Copyright (C) 2014, 2016 Christian Stigen Larsen
 * Distributed under the GPL v3 or later. See COPYING.
 */

#include "filesize.hpp"

#include <stdexcept>
#include <string.h>
#include <sys/stat.h>

namespace arv {

std::size_t filesize(const int file_descriptor)
{
  struct stat st;
  memset(&st, 0, sizeof(struct stat));

  if ( fstat(file_descriptor, &st) < 0 )
    throw std::runtime_error("Could not stat file");

  const off_t size = st.st_size;
  return size < 0 ? 0 : static_cast<std::size_t>(size);
}

}

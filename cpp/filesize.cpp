/*
 * Copyright (C) 2014, 2016 Christian Stigen Larsen
 * Distributed under the GPL v3 or later. See COPYING.
 */

#include "filesize.hpp"
#include <stdexcept>
#include <sys/stat.h>

namespace arv {

std::size_t filesize(const int file_descriptor)
{
  struct stat stat = {0};

  if ( fstat(file_descriptor, &stat) < 0 )
    throw std::runtime_error("Could not stat file");

  const off_t size = stat.st_size;
  return size < 0 ? 0 : static_cast<std::size_t>(size);
}

}

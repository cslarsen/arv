/*
 * Copyright (C) 2014, 2016 Christian Stigen Larsen
 * Distributed under the GPL v3 or later. See COPYING.
 */

#include <stdexcept>
#include <sys/stat.h>
#include "filesize.hpp"

off_t filesize(const int file_descriptor)
{
  struct stat stat;

  if ( fstat(file_descriptor, &stat) < 0 )
    throw std::runtime_error("Could not stat file");

  return stat.st_size;
}

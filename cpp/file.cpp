/*
 * Copyright (C) 2014, 2016 Christian Stigen Larsen
 * Distributed under the GPL v3 or later. See COPYING.
 */

#include <unistd.h>
#include <fcntl.h>
#include <stdexcept>
#include <string>
#include "file.hpp"

File::File(const char* filename, const int flags):
  fd(open(filename, flags))
{
  if ( fd < 0 ) {
    std::string msg = "Could not open ";
    throw std::runtime_error((msg + filename).c_str());
  }
}

File::~File() {
  close(fd);
}

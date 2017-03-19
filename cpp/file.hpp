/*
 * dna-traits
 * Copyright 2014, 2016 Christian Stigen Larsen
 * Distributed under the GPL v3 or later. See COPYING.
 *
 * arv
 * Copyright 2017 Christian Stigen Larsen
 * Distributed under the GNU GPL v3 or later. See COPYING.
 */

#ifndef DNA_FILE_H
#define DNA_FILE_H

#define BUILDING_DLL
#include "export.hpp"

class DLL_LOCAL File {
  int fd;
public:
  File(const char* filename, const int flags);
  ~File();

  inline operator int() const {
    return fd;
  }
};

#endif

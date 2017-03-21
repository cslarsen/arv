/*
 * arv
 * Copyright 2017 Christian Stigen Larsen
 * Distributed under the GNU GPL v3 or later. See COPYING.
 */

#ifndef ARV_CYTHON_HPP
#define ARV_CYTHON_HPP

#include "Python.h"
#include "export.hpp"

// We compile globally with -fvisibility=hidden. But Cython doesn't allow us to
// put visibility specifiers in the generated code, so we just do the trick
// here. And as long as Python can find the init function, it will be able to
// fetch all the required function addresses.
extern "C" DLL_PUBLIC PyMODINIT_FUNC init_arv(void);

#endif

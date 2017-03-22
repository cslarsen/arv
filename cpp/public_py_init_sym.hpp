/*
 * arv
 * Copyright 2017 Christian Stigen Larsen
 * Distributed under the GNU GPL v3 or later. See COPYING.
 */

#ifndef ARV_CYTHON_HPP
#define ARV_CYTHON_HPP

#include <Python.h>
#include "export.hpp"

// When compiling with hidden symbols (-fvisibility=hidden), we still need to
// make the init function's symbol global (i.e. public).

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC init_arv(void) DLL_PUBLIC;
#else
PyMODINIT_FUNC PyInit__arv(void) DLL_PUBLIC;
#endif

#endif // guard

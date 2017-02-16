#ifndef TEST_H
#define TEST_H

#include "types.h"
// This header/source pair contains some tests to see if
// the python interface behaves properly.

extern "C" {

void test_types( void *px, py_int N, void *pids, void *ptypes,
                 py_int periodic, py_float *xlo );

} // extern "C"

#endif // TEST_H

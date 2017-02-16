#ifndef LINE_TENSION_H
#define LINE_TENSION_H

#include "types.h"

#include <vector>
#include <list>

extern "C" {
py_float get_line_tension( void *x, py_int N, py_int *ids,
                           py_int *types, py_int periodic,
                           py_float *xlo, py_float *xhi, py_int dims,
                           py_int tstep, const char *boxline,
                           py_int nn_expect, py_float F_per_particle );
}

double estimate_line_tension( const class block_data &edge_block,
                              const std::vector<std::list<py_int> > &neighs,
                              double F_per_particle );



#endif // LINE_TENSION_H

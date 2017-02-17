#ifndef ENSEMBLE_GENERATOR_H
#define ENSEMBLE_GENERATOR_H

#ifdef HAVE_LIB_LAMMPS

#include "types.h"
#include "lmptype.h"


extern "C" {


/*!
  @brief Creates atom positions that lie on given manifold.

  @param x         Atom position vector
  @param N         Number of atoms
  @param ids       Atom ids
  @param types     Atom types
  @param periodic  Periodic boundary settings
  @param xlo       Box lower bounds
  @param xhi       Box upper bounds
  @param manifold  Name of the manifold to use
  @param rc        Minimum distance that should be between particles.
*/
void ensemble_generate_manifold( void *x, py_int N, py_int *ids,
                                 py_int *types, py_int periodic,
                                 py_float *xlo, py_float *xhi,
                                 const char *manifold, const char *man_args,
                                 py_float rc, py_int seed );
}


void ensemble_generate_impl( arr3f &x, py_int N, arr1i &ids, arr1i &types,
                             py_int periodic, py_float *xlo, py_float *xhi,
                             const char *mname, const char *man_args,
                             py_float rc, py_int seed );

#endif // HAVE_LIB_LAMMPS

#endif // ENSEMBLE_GENERATOR_H

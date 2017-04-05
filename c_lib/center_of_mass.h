#ifndef CENTER_OF_MASS_H
#define CENTER_OF_MASS_H

/*!
  \file center_of_mass.h
  @brief Calculates the centre-of-mass of atoms

  \ingroup cpp_lib
*/

#include "types.h"

extern "C" {
/*!
   @brief Computes the RDF of species itype with jtype from atom positions.

   @param x         Double array of atom positions.
   @param N         Number of atoms
   @param ids       Array of atom ids
   @param types     Array of atom types
   @param Ngroups   Number of groups
   @param groups    Specifies which atoms to count together as a group.
   @param xlo       Lower bounds of simulation domain
   @param xhi       Upper bounds of simulation domain
   @param periodic  Is simulation domain periodic?
   @param dims      Dimension of the system
   @param com       Store the centre of mass here.
 */
void center_of_mass( void *px, py_int N, py_int *pids, py_int *ptypes,
                     py_int Ngroups, py_int *pgroups, py_float *mass,
                     py_float *xlo, py_float *xhi, py_int periodic,
                     py_int dim, void *com  );

	
}



#endif // CENTER_OF_MASS

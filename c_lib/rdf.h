#ifndef RDF_H
#define RDF_H

/*!
  \file rdf.h
  @brief Contains some functions for computing radial distribution functions.

  \ingroup  cpp_lib
*/

#include "types.h"

#include <list>

extern "C" {
/*!
   @brief Computes the RDF of species itype with jtype from atom positions.

   @param x         Double array of atom positions.
   @param N         Number of atoms
   @param ids       Array of atom ids
   @param types     Array of atom types
   @param r0        Lower bound of where RDF is computed
   @param r1        Upper bound of where RDF is computed
   @param nbins     Number of bins to use. Resolution is dr = (r1-r0)/(nbins-1)
   @param itype     Type of atom 1 to consider for computing RDF (0 for all)
   @param jtype     Type of atom 2 to consider for computing RDF (0 for all)
   @param xlo       Lower bounds of simulation domain
   @param xhi       Upper bounds of simulation domain
   @param periodic  Is simulation domain periodic?
   @param dim       Dimension of the system (used in normalization)
   @param method    Method to use for neighborizing.
   @param rdf       Array to store the RDF in
   @param coord     Array to store the coordination number in
 */
void compute_rdf( void *px, py_int N, py_int *pids, py_int *ptypes,
                  py_float r0, py_float r1, py_int nbins, py_int itype,
                  py_int jtype, py_float *xlo, py_float *xhi, py_int periodic,
                  py_int dim, py_int method, py_float *prdf, py_float *pcoord );

/*!
   @brief Computes the ADF of species itype with jtype from atom positions on sphere of radius R.

   @param x         Double array of atom positions.
   @param N         Number of atoms
   @param ids       Array of atom ids
   @param types     Array of atom types
   @param r0        Lower bound of where RDF is computed
   @param r1        Upper bound of where RDF is computed
   @param nbins     Number of bins to use. Resolution is dr = (r1-r0)/(nbins-1)
   @param itype     Type of atom 1 to consider for computing RDF (0 for all)
   @param jtype     Type of atom 2 to consider for computing RDF (0 for all)
   @param R         Radius of spherical template.
   @param method    Method to use for neighborizing.
   @param adf       Array to store the ADF in
   @param coord     Array to store the coordination number in
 */
void compute_adf( void *px, py_int N, py_int *pids, py_int *ptypes,
                  py_int nbins, py_int itype, py_int jtype,
                  py_float R, py_int method, py_float *padf, py_float *pcoord );


/*!
  @brief A C++-side test of compute_rdf.
*/
void test_rdf();

} // extern "C"


/*!
  \private
  @brief Actual implementation of the RDF computation.
*/
void compute_rdf_impl( const arr3f &x, py_int N, const arr1i &ids,
                       const arr1i &types, py_float x0, py_float x1,
                       py_int nbins, py_int itype, py_int jtype,
                       py_float *xlo, py_float *xhi,
                       py_int periodic, py_int dim, std::list<py_int> *neighs,
                       arr1f &ardf, arr1f &acoord );
/*!
  \private
  @brief Actual implementation of ADF computation
*/
void compute_adf_impl( const arr3f &x, py_int N, const arr1i &ids,
                       const arr1i &types, py_int nbins, py_int itype, py_int jtype,
                       py_float R, std::list<py_int> *neighs, arr1f &aadf, arr1f &acoord );

#endif /* RDF_H */

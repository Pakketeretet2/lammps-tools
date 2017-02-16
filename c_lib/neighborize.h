#ifndef NEIGHBORIZE_H
#define NEIGHBORIZE_H

/*!
  \file neighborize.h
  @brief Contains some functions for determining neighbor lists.

  \ingroup  cpp_lib
*/

#include "list.h"
#include "types.h"

#include <vector>

/*!
  Enum of the neighborization methods.
*/
enum NEIGHBORIZE_METHODS {
	DIST_NSQ    = 0, ///< A distance criterion, using a method slow for big data sets
	DIST_BIN    = 1, ///< A distance criterion, using a method fast for big data sets
        DELAUNAY    = 2, ///< Delaunay triangulation in 2/3D, using the CGAL lib
	CONVEX_HULL = 3  ///< Convex hull for particles on ellipsoids/spheres, using CGAL lib
};


extern "C" {

/*!
  @brief Determine some form of neighbor list of a group of particles.

  @param x         Atom positions
  @param N         Number of atoms
  @param ids       Atom ids
  @param types     Atom types
  @param rc        Maximum distance between atoms to be considered neighbors
  @param periodic  Periodic boundary settings
  @param xlo       Box lower bounds
  @param xhi       Box upper bounds
  @param dims      Box dimensions
  @param method    Neighborisation method (see NEIGHBORIZE_METHODS)
  @param pname     Name of the pipe the neighbor lists should be written to
*/
void neighborize( void *x, py_int N, py_int *ids,
                  py_int *types, py_float rc, py_int periodic,
                  py_float *xlo, py_float *xhi, py_int dims,
                  py_int method, const char *pname,
                  py_int itype, py_int jtype );


}

void neighborize_block( const class block_data &b,
                        std::vector<std::list<py_int> > &neighs );


/*!
  @brief Calls the correct neighborize function on correctly typed data.

  \private

  @param x         Atom positions
  @param N         Number of atoms
  @param ids       Atom ids
  @param types     Atom types
  @param rc        Maximum distance between atoms to be considered neighbors
  @param periodic  Periodic boundary settings
  @param xlo       Box lower bounds
  @param xhi       Box upper bounds
  @param dims      Box dimensions
  @param method    Neighborisation method (see NEIGHBORIZE_METHODS)
  @param neighs    Pointer to where the neigh list will be stored
*/
void neighborize_impl( const arr3f &x, py_int N, const arr1i &ids,
                       const arr1i &types, py_float rc, py_int periodic,
                       const py_float *xlo, const py_float *xhi, py_int dims,
                       py_int method, list *neighs,
                       py_int itype, py_int jtype );


/*!
  @brief Computes a distance-based neighbor list by computing all distances
         r_{ij} and filtering those with r_{ij} <= rc

  \private

  @param x         Atom positions
  @param N         Number of atoms
  @param ids       Atom ids
  @param types     Atom types
  @param rc        Maximum distance between atoms to be considered neighbors
  @param periodic  Periodic boundary settings
  @param xlo       Box lower bounds
  @param xhi       Box upper bounds
  @param dims      Box dimensions
  @param method    Neighborisation method (see NEIGHBORIZE_METHODS)
  @param neighs    Pointer to where the neigh list will be stored
*/
void neighborize_dist_nsq( const arr3f &x, py_int N, const arr1i &ids,
                           const arr1i &types, py_float rc, py_int periodic,
                           const py_float *xlo, const py_float *xhi, py_int dims,
                           list *neighs, py_int itype, py_int jtype  );

/*!
  @brief Computes a distance-based neighbor list by first binning all
         positions, determining which atoms are likely close and then
         computing distance r_{ij}, only for those atoms, then filtering
         those with r_{ij} <= rc. This is faster for large systems.
  \private

  @param x         Atom positions
  @param N         Number of atoms
  @param ids       Atom ids
  @param types     Atom types
  @param rc        Maximum distance between atoms to be considered neighbors
  @param periodic  Periodic boundary settings
  @param xlo       Box lower bounds
  @param xhi       Box upper bounds
  @param dims      Box dimensions
  @param neighs    Pointer to where the neigh list will be stored
*/
void neighborize_dist_bin( const arr3f &x, py_int N, const arr1i &ids,
                           const arr1i &types, py_float rc, py_int periodic,
                           const py_float *xlo, const py_float *xhi, py_int dims,
                           list *neighs, py_int itype, py_int jtype  );


/*!
  @brief Computes a neighbor list based on Delaunay triangulation using
         the CGAL library.
  \private

  @param x         Atom positions
  @param N         Number of atoms
  @param ids       Atom ids
  @param types     Atom types
  @param xlo       Box lower bounds
  @param xhi       Box upper bounds
  @param dims      Box dimensions
  @param neighs    Pointer to where the neigh list will be stored
*/
void neighborize_delaunay( const arr3f &x, py_int N, const arr1i &ids,
                           const arr1i &types, py_int periodic,
                           const py_float *xlo, const py_float *xhi, py_int dims,
                           list *neighs, py_int itype, py_int jtype );

/*!
  @brief Computes a connectivity list by computing the convex hull of the
         particles using the CGAL library. This results in Delaunay-like
         triangulation for particles on a sphere/ellipsoid.
  \private

  @param x         Atom positions
  @param N         Number of atoms
  @param ids       Atom ids
  @param types     Atom types
  @param periodic  Periodic boundary settings
  @param xlo       Box lower bounds
  @param xhi       Box upper bounds
  @param dims      Box dimensions
  @param neighs    Pointer to where the neigh list will be stored
*/
void neighborize_conv_hull( const arr3f &x, py_int N, const arr1i &ids,
                            const arr1i &types, py_int periodic,
                            const py_float *xlo, const py_float *xhi, py_int dims,
                            list *neighs, py_int itype, py_int jtype  );





#endif /* NEIGHBORIZE_H */

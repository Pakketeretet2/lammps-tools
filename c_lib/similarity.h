#ifndef SIMILARITY_H
#define SIMILARITY_H

#include "types.h"

extern "C" {

/*
 * Finds a rotation axis and angle that best maps the atoms in xi onto
 * those in xj. It stores the best matches positions in rotated, the
 * rotation axis in rot_axis and returns the rotation angle.
 */

double map_clusters( void *xi, void *xj,  py_int N,
                     py_int *ids, py_int *jds, py_int periodic,
                     void *xlo, void *xhi, py_int dims,
                     void *rotated, py_float *rot_axis );


double rotate_to_template( void *xxi, void *xxt,  py_int N,
                           py_int periodic, void *xxlo, void *xxhi,
                           py_int dims, void *mmaj_ax, void *mmin_ax,
                           void *mmaj_ax_t, void *mmin_ax_t,
                           void *xrot = nullptr );

	

} // extern "C"

double map_clusters_impl( const arr3f &xi, const arr3f &xj, py_int N,
                          py_int *ids, py_int *jds, py_int periodic,
                          py_float *xlo, py_float *xhi, py_int dims,
                          arr3f &rotated, py_float *rot_axis );

double exhaust_rotate_search( const arr3f &xi, const arr3f &xj, py_int N,
                              py_int periodic, py_float *xlo, py_float *xhi,
                              py_int dims, arr3f &rotated, py_float *rot_axis );

double rotate_to_template_impl( const arr3f &xi, const arr3f &xt, py_int N,
                                py_int periodic, py_float *xlo, py_float *xhi,
                                py_int dims, py_float *maj_ax, py_float *min_ax,
                                py_float *maj_ax_t, py_float *min_ax_t,
                                py_float *xrot );

void test_rotations();

#endif // SIMILARITY_H

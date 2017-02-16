#ifndef SKELETONIZE_H
#define SKELETONIZE_H

#include <vector>
#include <list>
#include <functional>

#include "types.h"


struct ribbon_data
{
	double R;
	
	std::size_t n_largest;
	double avg_r, var_r, rmin, rmax, max_avg, max_var;

	double area;
};


// C interface:
extern "C" {
void get_insideness( void *x, py_int N, py_int *ids,
                     py_int *types, py_int periodic,
                     py_float *xlo, py_float *xhi, py_int dims,
                     py_float *insideness );
	
void get_euclidian_distance_transform( void *x, py_int N, py_int *ids,
                                       py_int *types, py_int periodic,
                                       py_float *xlo, py_float *xhi,
                                       py_int dims, py_float R, py_float *iins,
                                       py_float *edt );
} // extern "C"



std::vector<double> get_insideness( const class block_data &b,
                                    std::vector<std::list<py_int> > *neigh_ptr = nullptr );


std::vector<double> euclidian_distance_transform( const class block_data &b,
                                                  const std::vector<double> &insideness,
                                                  double R );

void get_edge( const block_data &b, const std::vector<double> &insideness,
               block_data &b_edge );


void skeletonize_edt( const class block_data &b, std::vector<py_int> &skeleton,
                      const std::vector<double> &insideness, double R );

void skeletonize( const block_data &b, std::vector<py_int> &skeleton,
                  const std::vector<double> &distance_map, double R );


void get_local_maxima( const std::vector<std::list<py_int> > &neighs,
                       const std::vector<double> &field,
                       std::vector<py_int> &max_indices );

void get_local_maxima( const std::vector<std::list<py_int> > & neighs,
                       const std::vector<double> &field,
                       std::vector<py_int> &max_indices,
                       std::function< double(double) > );


void get_ribbon_data( class block_data &b, ribbon_data &r_data );
	
#endif // SKELETONIZE_H

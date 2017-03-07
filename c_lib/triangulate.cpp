#include "triangulate.h"
#include "neighborize.h"
#include "communication.h"
#include "util.h"
#include "dump_reader.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>


template <typename T> inline
T my_min( T a, T b )
{
	return (a < b) ? a : b;
}

template <typename T> inline
T my_max( T a, T b )
{
	return (a >= b) ? a : b;
}




void triangulate_impl(const arr3f &, py_int , py_int *, py_int *,
                      py_float , py_int, const py_float *, const py_float *,
                      py_int , py_int , const char * );

extern "C" {

void triangulate( void *x, py_int N, py_int *ids, py_int *types,
                  py_float rc, py_int periodic,
                  const py_float *xlo, const py_float *xhi, py_int dims,
                  py_int method, const char *pname )
{
	arr3f xi(x,3);
	triangulate_impl( xi, N, ids, types, rc, periodic, xlo, xhi,
	                  dims, method, pname );
}

}

template <typename T> inline
T my_max( T* x, py_int size )
{
	T x0 = x[0];
	for( py_int i = 1; i < size; ++i ){
		if( x0 < x[i] ){
			x0 = x[i];
		}
	}
	return x0;
}

double vec_dist( const arr3f &x, int i, int j )
{
	double dx = x[i][0] - x[j][0];
	double dy = x[i][1] - x[j][1];
	double dz = x[i][2] - x[j][2];
	return std::sqrt( dx*dx + dy*dy + dz*dz );
}

int insert_triangle( const arr3f &x, int i, int j, int k, int **out,
                      std::vector<triangle> &triangles, list *neighs )
{
	bool added = 0;
	if( (i < k) && (j < k) ){
		if( list_has(neighs[i], k) ){
			double rij = vec_dist( x, i, j );
			double rik = vec_dist( x, i, k );
			double rjk = vec_dist( x, j, k );
			std::cerr << "(i,j,k) = ( " << i << ", " << j
			          << ", " << k << " ).\n";
			std::cerr << "xi = ( " << x[i][0] << ", " << x[i][1]
			          << ", " << x[i][2] << " ).\n";
			std::cerr << "xi = ( " << x[j][0] << ", " << x[j][1]
			          << ", " << x[j][2] << " ).\n";
			std::cerr << "xi = ( " << x[k][0] << ", " << x[k][1]
			          << ", " << x[k][2] << " ).\n";
			
			std::cerr << "dists: " << rij << ", " << rik
			          << ", " << rjk << "...\n";
			
			out[i][j]++;
			out[i][k]++;
			out[j][k]++;
			
			triangle tri( i, j, k, x[i],
			              x[j], x[k] );
			triangles.push_back( tri );
			added++;
		}
	}
	return added ? 1 : 0;
}

void triangulate_impl( const arr3f &x, py_int N, py_int *ids, py_int *types,
                       py_float rc, py_int periodic, const py_float *xlo,
                       const py_float *xhi, py_int dims, py_int method,
                       const char *pname )
{
	// Make a neighbour list first:
	arr1i iids(ids,N);
	arr1i ttypes(types,N);

	list *neighs = new list[N];
	
	neighborize_impl( x, N, iids, ttypes, rc, periodic, xlo, xhi,
	                  dims, method, neighs, 0, 0 );



	std::vector<py_int> v_ids( ids, ids + N );
	py_int max_id = *std::max_element( v_ids.begin(), v_ids.end() );
	std::vector<py_int> id_map( max_id + 1 );

	std::vector<triangle> triangles;
	
	for( py_int i = 0; i < N; ++i ){
		id_map[ ids[i] ] = i;
	}

	int N2 = N*N;
	int *idx_out_  = new int[N*N];
	int **idx_out  = new int*[N];
	
	
	for( int i = 0; i < N; ++i ){
		idx_out[i] = idx_out_ + N*i;
		for( int j = 0; j < N; ++j ){
			idx_out_[i + j*N] = 0;
		}
	}

	int added = 0;
	for( py_int ii = 0; ii < N; ++ii ){
		// loop over all neighs of i1:

		const list &l = neighs[ii];
		for( py_int j : l ){
			if( (ii < j) && (idx_out[ii][j] < 3) ){
				// Now you got two ids... ids[i] and j. All you need
				// to know now is if neighs[j] contains ids that are
				// also in neighs[i]! However, if i and j were already
				// added twice, you can safely skip them.
				for( py_int k : neighs[j] ){
					if( (ii < k) && (j < k) &&
					    list_has( neighs[k], ii ) ){
						
						// Add triangle:
						triangles.push_back(
							triangle(ii,j,k,x[ii],
							         x[j],x[k]) );
						++added;
					}
				}
			}
			
		}
	}


	std::ofstream pout( pname );
	for( const triangle &t : triangles ){
		pout << t.i1 << " " << t.i2 << " " << t.i3 << "\n";
	}
	close_pipe( pname );

	
	delete [] idx_out_;
	delete [] idx_out;
	
	delete [] neighs;
}


void triangulate_block( block_data &b, py_float rc, py_int periodic,
                        py_int dims, py_int method, std::vector<triangle> &triangles )
{
	// Make a neighbour list first:
	std::size_t N = b.N;
	
	arr1i iids(b.ids ,N);
	arr1i ttypes(b.types ,N);
	arr3f xx(b.x, N);
	list *neighs = new list[N];
	
	neighborize_impl( xx, N, iids, ttypes, rc, periodic,
	                  b.xlo, b.xhi, dims, method, neighs, 0, 0 );
	std::vector<py_int> v_ids( b.ids, b.ids + N );
	py_int max_id = *std::max_element( v_ids.begin(), v_ids.end() );
	std::vector<py_int> id_map( max_id + 1 );
	
	for( py_int i = 0; i < N; ++i ){
		id_map[ b.ids[i] ] = i;
	}

	int N2 = N*N;
	int *idx_out_  = new int[N*N];
	int **idx_out  = new int*[N];
	
	
	for( int i = 0; i < N; ++i ){
		idx_out[i] = idx_out_ + N*i;
		for( int j = 0; j < N; ++j ){
			idx_out_[i + j*N] = 0;
		}
	}

	int added = 0;
	for( py_int ii = 0; ii < N; ++ii ){
		// loop over all neighs of i1:

		const list &l = neighs[ii];
		for( py_int j : l ){
			if( (ii < j) && (idx_out[ii][j] < 3) ){
				// Now you got two ids... ids[i] and j. All you need
				// to know now is if neighs[j] contains ids that are
				// also in neighs[i]! However, if i and j were already
				// added twice, you can safely skip them.
				for( py_int k : neighs[j] ){
					if( (ii < k) && (j < k) &&
					    list_has( neighs[k], ii ) ){
						
						// Add triangle:
						triangles.push_back(
							triangle(ii,j,k,xx[ii],
							         xx[j],xx[k]) );
						++added;
					}
				}
			}
			
		}
	}


	delete [] idx_out_;
	delete [] idx_out;
	
	delete [] neighs;
}


double triangulation_area( block_data &b, std::vector<triangle> &triangles )
{
	double A = 0.0;
	for( const triangle &t : triangles ){
		A += t.area();
	}
	return A;
}

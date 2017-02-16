#include "skeletonize.h"
#include "dump_reader.h"
#include "neighborize.h"
#include "id_map.h"
#include "triangulate.h"

#include <cmath>
#include <list>
#include <algorithm>

struct graph_vertex
{
	graph_vertex() : i(-1), pt{0,0,0}
	{}
	graph_vertex(int i, const std::array<double,3> &p) : i(i), pt(p)
	{}
	graph_vertex(int i, double x, double y, double z) : i(i), pt{x,y,z}
	{}

	std::size_t i;
	std::array<double,3> pt;
};






double sphere_dist( const std::array<double,3> &xi, const std::array<double,3> &xj,
                    double R )
{
	double R2 = R * R;
	double dot = xi[0]*xj[0] + xi[1]*xj[1] + xi[2]*xj[2];
	double theta = std::acos( dot / R2 );
	return theta*R;
	
}

std::vector<double> get_insideness( const block_data &b,
                                    std::vector<std::list<py_int> > *neigh_ptr )
{
	std::vector< std::list<py_int> > nneighs;
	std::vector< std::list<py_int> > *neigh_ptr2 = &nneighs;
	if( neigh_ptr ){
		neigh_ptr2 = neigh_ptr;
	}
	std::vector< std::list<py_int> > &neighs = *neigh_ptr2;

	if( !neigh_ptr ){
		neighborize_block( b, neighs );
	}
	
	std::vector<double> insideness( b.N );
	for( int i = 0; i < b.N; ++i ){
		insideness[i] = -1.0;
	}

	
	for( py_int i = 0; i < b.N; ++i ){
		if( neighs[i].size() < 6 ){
			insideness[i] = 0.0;
		}
	}
	
	// Stategy: Recursively loop over all points in network. Those that
	// are not yet determined to be anywhere are assigned the lowest value
	// of their neighbouring insideness plus one.
	id_map im( b.ids );	
	bool assigned_one = false;
	double current_val = 0.0;
	// std::cerr << "Finding insideness. At loop 0...";
	do{
		assigned_one = false;
		for( int i = 0; i < b.N; ++i ){
			if( insideness[i] >= 0 ){
				continue;
			}
			py_int idi = b.ids[i];

			const std::list<py_int> &ni = neighs[i];
			bool has_current_val = false;
			for( py_int idx : ni ){
				int idj = b.ids[idx];
				if( insideness[idx] == current_val ){
					has_current_val = true;
					break;
				}
			}
			if( has_current_val ){
				insideness[i] = current_val + 1.0;
				assigned_one = true;
			}
		}
		current_val += 1.0;
		// std::cerr << " " << current_val << "...";
	}while( assigned_one );

	// std::cerr << "Done assigning...\n";


	return insideness;
}


std::vector<double> euclidian_distance_transform( const class block_data &b,
                                                  const std::vector<double> &insideness,
                                                  double R )
{
	std::vector<double> edt( b.N );
	block_data b_edge;
	get_edge( b, insideness, b_edge );
	
	// Generate a separate list that contains only the edge.	
	
	b_edge.print( "edge.dump" );
	for( int i = 0; i < b.N; ++i ){
		if( insideness[i] == 0 ){
			edt[i] = 0.0;
		}else{
			const std::array<double,3> &xi = b.x[i];
			double R2 = R*R;
			double min_dist = 16*R;
			
			// Find the closest edge:
			for( const std::array<double,3> &xe : b_edge.x ){
				double dot = xi[0]*xe[0] + xi[1]*xe[1] + xi[2]*xe[2];
				double theta = std::acos( dot / R2 );
				min_dist = std::min( min_dist, theta*R );
			}
			edt[i] = min_dist;
		}
	}
	
	
	return edt;
}


void skeletonize_edt( const class block_data &b, std::vector<py_int> &skeleton,
                      const std::vector<double> &insideness, double R )
{
	std::vector<double> edt = euclidian_distance_transform( b, insideness, R );
	skeletonize( b, skeleton, edt, R );
}



void skeletonize_alg1( const block_data &b, std::vector<py_int> &skeleton,
                       const std::vector<double> &distance_map, double R )
{
	// Algorithm:
	// 0. Find point with largest distance map, i
	// 1. Find neighbour of that point j for which
	//    (DM(i) - DM(j)) / (xi - xj) is minimal.

	
	std::vector<bool> out( distance_map.size() );
	for( int i = 0; i < out.size(); ++i ){
		out[i] = false;
	}
	int i = 0;
	double d_max = 0.0;
	for( int j = 0; j < distance_map.size(); ++j ){
		if( distance_map[j] > d_max ){
			d_max = distance_map[j];
			i = j;
		}
	}
	std::cerr << "max dist = " << d_max << " for particle " << i << ".\n";
	out[i] = true;
	std::vector<std::list<py_int> > neighs;
	neighborize_block( b, neighs );
	skeleton.push_back(i);
	
	bool got_one_in = false;
	do{
		// Check the particles in the neighbourhood of i
		// that are not yet out, and find the one that has
		// the most positive gradient in EDT.
		for( int j : neighs[i] ){
			double max_grad = -1000;
			int    max_j = -1;
			got_one_in = false;
			if( !out[j] ){
				// Compute the gradient in EDT:
				double d_edt = distance_map[j] - distance_map[i];
				double d_x   = sphere_dist( b.x[i], b.x[j], R );
				double grad = d_edt / d_x;
				std::cerr << "Gradient between " << i << " and "
				          << j << " = " << grad << ".\n";
				got_one_in = true;
				if( grad > max_grad ){
					max_grad = grad;
					max_j = j;
				}
			}
			// Select j as the new i.
			skeleton.push_back(j);
			out[j] = true;
			i = j;
		}
	}while( got_one_in );
	
}

void skeletonize_alg2( const block_data &b, std::vector<py_int> &skeleton,
                       const std::vector<double> &distance_map, double R )
{
	// Another idea.
	// 1. Sort the points in terms of distance_map.
	// 2. For each point with more than one bond that is not an edge
	//    point, delete all except two bonds. Keep the bond with the
	//    largest EDT, and the one the most opposite to it.

	std::vector<std::size_t> permutation( distance_map.size() );
	for( int i = 0; i < permutation.size(); ++i ){
		permutation[i] = i;
	}
	std::sort( permutation.begin(), permutation.end(),
	           [&permutation,distance_map]( std::size_t a,
	                                        std::size_t b )
	           {
		           return distance_map[a] < distance_map[b];
	           } );
	for( std::size_t idx : permutation ){
		const std::array<double,3> &xi = b.x[idx];
		
	}
	
}

void skeletonize_cgal( const block_data &b, std::vector<py_int> &skeleton,
                       const std::vector<std::list<py_int> > &neighs )
{
	// Construct a triangulated mesh of the points:
	for( int i = 0; i < b.N; ++i ){
		
	}
}


void skeletonize( const block_data &b, std::vector<py_int> &skeleton,
                  const std::vector<double> &distance_map, double R )
{
	skeletonize_alg2( b, skeleton, distance_map, R );
}

void get_edge( const block_data &b, const std::vector<double> &insideness,
               block_data &b_edge )
{
	int N_edge = 0;
	for( double i : insideness ){
		if( i == 0 ) ++N_edge;
	}

	
	b_edge.init_empty( N_edge );

	int j = 0;
	for( int i = 0; i < insideness.size(); ++i ){
		if( insideness[i] == 0 ){
			b_edge.x[j]     = b.x[i];
			b_edge.ids[j]   = b.ids[i];
			b_edge.types[j] = b.types[i];
			++j;
		}
	}
}

void get_local_maxima( const std::vector<std::list<py_int> > &neighs,
                       const std::vector<double> &field,
                       std::vector<py_int> &max_indices )
{
	auto unit_operator = []( double x ){ return x; };
	get_local_maxima( neighs, field, max_indices, unit_operator );
}


void get_local_maxima( const std::vector<std::list<py_int> > & neighs,
                       const std::vector<double> &field,
                       std::vector<py_int> &max_indices,
                       std::function< double(double) > f )
{
	// d_out << "Indices of maxima:";
	for( std::size_t idx = 0; idx < neighs.size(); ++idx ){
		double vi = f( field[idx] );
		bool largest = true;
		for( std::size_t idj : neighs[idx] ){
			if( f( field[idj] ) >= vi ){
				largest = false;
				break;
			}
		}
		if( neighs[idx].empty() ){
			largest = false;
		}
		if( largest ){
			// d_out << " " << idx;
			max_indices.push_back( idx );
		}
	}
	// d_out << "\n";
}




template <typename T>
std::vector<std::size_t> get_sorted_permutation( const std::vector<T> &c )
{
	std::vector<std::size_t> idx( c.size() );
	
	for( std::size_t i = 0; i < c.size(); ++i ){
		idx[i] = i;
	}
	auto comp = [&]( std::size_t i, std::size_t j ) { return c[i] < c[j]; };
	std::sort( idx.begin(), idx.end(), comp );
	
	return idx;
}

template <typename container>
void mean_and_var( const container &c, double &m, double &v )
{
	m = v = 0.0;
	double n = 0.0;
	for( double value : c ){
		m += value;
		n += 1.0;
	}
	m /= n;
	for( double value : c ){
		double dx = value - m;
		v += dx*dx;
	}
	v /= (n-1);
}

void get_stats( const std::vector<double> &values,
                const std::vector<py_int> &idx,
                double &avg, double &var, double &min_val, double &max_val,
		std::vector<double> &largest )
{
	double n = 0.0;
	avg = var = min_val = max_val = 0.0;
	std::vector<double> sorted( idx.size() );
	std::size_t j = 0;
	for( std::size_t i : idx ){
		sorted[j] = values[i];
		j++;
	}

	mean_and_var( sorted, avg, var );

	std::sort( sorted.begin(), sorted.end() );
	min_val = sorted[0];
	max_val = sorted[sorted.size()-1];
	for( int i = 0; i < largest.size(); ++i ){
		largest[i] = sorted[sorted.size()-1-i];
	}
}



void get_ribbon_data( block_data &b, ribbon_data &r_data )
{
	std::cerr << "This is get_ribbon_data...\n";	
	std::vector<std::list<py_int> > neighs;
	neighborize_block( b, neighs );

	std::vector<double> insideness = get_insideness( b, &neighs );
	std::vector<double> edt = euclidian_distance_transform( b, insideness,
	                                                        r_data.R );
	std::vector<py_int> max_edts;
	get_local_maxima( neighs, edt, max_edts );

	double rmax  = 0.0;
	double rmin  = 0.0;
	double avg_r = 0.0;
	double var_r = 0.0;

	std::vector<double> largest( r_data.n_largest );
	get_stats( edt, max_edts, avg_r, var_r, rmin, rmax, largest );

	double mm, vv;
	mean_and_var( largest, mm, vv );

	r_data.avg_r = avg_r;
	r_data.var_r = var_r;
	r_data.rmin = rmin;
	r_data.rmax = rmax;
	r_data.max_avg = mm;
	r_data.max_var = vv;

	std::vector<triangle> triangles;
	triangulate_block( b, 1.3, 0, 3, 0, triangles );
	std::cerr << "Triangulated " << b.N << " particles into mesh of "
	          << triangles.size() << " triangles.\n";
	r_data.area = triangulation_area( b, triangles );
}



void get_insideness( void *x, py_int N, py_int *ids,
                     py_int *types, py_int periodic,
                     py_float *xlo, py_float *xhi, py_int dims,
                     py_float *insideness )
{
	block_data b;
	block_data_from_foreign( x, N, ids, types, periodic, xlo, xhi, dims,
	                         0, "", b );
	std::vector<double> ins = get_insideness( b );
	for( py_int i = 0; i < N; ++i ){
		insideness[i] = ins[i];
	}
}

void get_euclidian_distance_transform( void *x, py_int N, py_int *ids,
                                       py_int *types, py_int periodic,
                                       py_float *xlo, py_float *xhi,
                                       py_int dims, py_float R, py_float *iins,
                                       py_float *edt )
{
	block_data b;
	block_data_from_foreign( x, N, ids, types, periodic, xlo, xhi, dims,
	                         0, "", b );
	
	std::vector<py_float> ins( N );
	for( std::size_t i = 0; i < N; ++i ){
		ins[i] = iins[i];
	}
	std::vector<double> eedt = euclidian_distance_transform( b, ins, R );
	
	for( py_int i = 0; i < N; ++i ){
		edt[i] = eedt[i];
	}
}

#include "neighborize.h"
#include "types.h"
#include "domain.h"
#include "my_timer.hpp"
#include "my_output.hpp"
#include "dump_reader.h"


#include <iostream>
#include <fstream>
#include <list>

#include <cassert>
#include <cmath>
#include <algorithm>

static my_ostream my_out( std::cout );


extern "C" {

	
// \note: This function assumes that the neigh list is build based on
//        atom IDS, not indices!
void neighborize( void *x, py_int N, py_int *ids,
                  py_int *types, py_float rc, py_int periodic,
                  py_float *xlo, py_float *xhi, py_int dims,
                  py_int method, const char *pname,
                  py_int itype, py_int jtype  )
{
	if( !x ){
		std::cerr << "Error! x was NULL!\n";
		return;
	}
	
	arr1i tmp_ids(ids,N);
	py_int max_id = *std::max_element( tmp_ids.begin(), tmp_ids.end());
	std::list<py_int> *neigh_list = new std::list<py_int>[N];

	arr3f xx(x,N);
	arr1i iids(ids,N), ttypes(types,N);

	neighborize_impl( xx, N, iids, ttypes, 
	                  rc, periodic, xlo, xhi, dims, method, neigh_list,
	                  itype, jtype );

	
	std::size_t s = 0;
	for( py_int i = 0; i < N; ++i ){
		s += neigh_list[i].size();
	}
	s += N;
	py_int sent = 0;
	std::size_t byte_one = sizeof(py_int);
	std::size_t bytes = (s+N)*byte_one; // s ids and N separators.

	// Now send the data through the pipe in the following format:
	// ID0 ID ID ID ID ID -1 ID1 ID ID ID etc..., with
	// IDi the ID of the ith atom, and the other IDS its neighbors.
	// The -1 is used to separate them.
	std::ofstream pout( pname, std::ios::binary );

	// Needs to be a py_int, not a regular one:
	py_int sep = static_cast<py_int>(-1);
	for( py_int i = 0; i < N; ++i ){
		union py_int_char {
			py_int val;
			char   bytes[sizeof(py_int)];
		};
		py_int_char ic;
		ic.val = ids[i];

		pout.write( ic.bytes, sizeof(py_int) );
		sent++;
		// pout << ids[i];
		for( py_int j : neigh_list[i] ){
			// pout << j;
			ic.val = ids[j];
			pout.write( ic.bytes, sizeof(py_int));
			sent++;
		}
		ic.val = sep;
		pout.write( ic.bytes, sizeof(py_int) );
		sent++;
	}
	
	pout.close();
	delete [] neigh_list;
}



	



} // extern "C"



void neighborize_impl( const arr3f &x, py_int N, const arr1i &ids,
                       const arr1i &types, py_float rc, py_int periodic,
                       const py_float *xlo, const py_float *xhi, py_int dims,
                       py_int method, std::list<py_int> *neighs,
                       py_int itype, py_int jtype  )
{
	// std::cerr << "Using neighbour method " << method << "\n";
	switch(method){
		case DIST_NSQ:
			/*
			std::cerr << "For some reason method == 0 breaks. Use "
			          << "method = 1 instead for dist-based "
			          << "neighborizing!\n";
			break;
			*/
			neighborize_dist_nsq( x, N, ids, types, rc, periodic,
			                      xlo, xhi, dims, neighs, itype, jtype );
			break;
		case DIST_BIN:
			neighborize_dist_bin( x, N, ids, types, rc, periodic,
			                      xlo, xhi, dims, neighs, itype, jtype );
			break;

#ifdef HAVE_LIB_CGAL
		case DELAUNAY:
			neighborize_delaunay( x, N, ids, types, periodic,
			                      xlo, xhi, dims, neighs, itype, jtype );
			break;
		case CONVEX_HULL:
			neighborize_conv_hull( x, N, ids, types, periodic,
			                       xlo, xhi, dims, neighs, itype, jtype );
			break;
#else
		case DELAUNAY:
		case CONVEX_HULL:
			std::cerr << "Delaunay and convex hull cannot be "
			          << "used without CGAL!\n";
			break;
#endif // HAVE_LIB_CGAL
		default:
			std::cerr << "Method " << method << " not recognized!\n";
			break;
	}

}








void neighborize_dist_nsq_impl( const arr3f &x, py_int N, const arr1i &ids,
                                const arr1i &types, py_float rc, py_int periodic,
                                const py_float *xlo, const py_float *xhi, py_int dims,
                                std::list<py_int> *neighs, py_int itype, py_int jtype  )
{
	
	double rc2 = rc*rc;
	for( py_int ii = 0; ii < N; ++ii ){
		py_int i  = ids[ii];
		
		for( py_int jj = 0; jj < N; ++jj ){
			py_int j  = ids[jj];

			if( i  >= j ) continue;

			double r[3];
			
			if( periodic ){
				distance_wrap( r, x[ii], x[jj], xlo, xhi );
			}else{
				distance( r, x[ii], x[jj] );
			}
			
			double r2 = r[0]*r[0] + r[1]*r[1];
			if( dims != 2 ) r2 += r[2]*r[2];

			if( r2 > rc2 ) continue;
			/*
			my_out << "Distance between " << i << " and " << j
			<< ", (" << x[ii][0] << ", "
			<< x[ii][1] << ", " << x[ii][2] << ") and ("
			<< x[jj][0] << ", " << x[jj][1] << ", "
			<< x[jj][2] << ") = " << sqrt(r2)
			<< ", so adding...\n";
			*/
			if( itype != 0 && types[ii] != itype ){
				continue;
			}
			if( jtype != 0 && types[jj] != jtype ){
				continue;
			}

			/* Add both to list. */
			neighs[ii].push_back(jj);
			neighs[jj].push_back(ii);
		}
	}
}


void neighborize_dist_nsq( const arr3f &x, py_int N, const arr1i &ids,
                           const arr1i &types, py_float rc, py_int periodic,
                           const py_float *xlo, const py_float *xhi, py_int dims,
                           std::list<py_int> *neighs, py_int itype, py_int jtype )
{
	neighborize_dist_nsq_impl( x, N, ids, types, rc,
	                           periodic, xlo, xhi, dims, neighs, itype, jtype );
}


biguint shift_bin_index( biguint i0, py_int xinc, py_int yinc, py_int zinc,
                         py_int Nx, py_int Ny, py_int Nz, py_int periodic )
{
	py_int binx = i0 % Nx;
	py_int biny = ( (i0 - binx) % (Nx*Ny) ) / Nx;
	py_int binz = (i0 - binx - biny*Nx) / (Nx*Ny);

	binx += xinc;
	biny += yinc;
	binz += zinc;

	if( periodic ){
		if( binx < 0 )    binx += Nx;
		if( binx >=  Nx ) binx -= Nx;

		if( biny < 0 )    biny += Ny;
		if( biny >= Ny )  biny -= Ny;

		if( binz < 0 )    binz += Nz;
		if( binz >= Nz )  binz -= Nz;

		assert( (binx >= 0) && (biny >= 0) && (binz >= 0) &&
		        (binx < Nx) && (biny < Ny) && (binz < Nz) &&
		        "Illegal bin index after correction!" );
	}


	return binx + Nx*biny + Nx*Ny*binz;

}


void add_neighs_from_bin( py_int i, const arr3f &x, double rc, const arr1i &ids,
                          py_int periodic, const py_float *xlo, const py_float *xhi,
                          py_int dim, const std::list<py_int> &bin, std::list<py_int> *neighs,
                          py_int itype, py_int jtype, const arr1i &types )
{
	/*
	  Check all particles inside this bin, add them if
	  they are within rc:
	*/
	double rc2 = rc*rc;
	
	for( py_int j : bin ){
		if( j != i ){
			const py_float *xi = x[i];
			const py_float *xj = x[j];
			py_float r[3];
			if( periodic ) distance_wrap(r,xi,xj,xlo,xhi);
			else           distance(r,xi,xj);

			double r2 = r[0]*r[0] + r[1]*r[1];
			if( dim != 2 ){
				r2 += r[2]*r[2];
			}

			if( itype != 0 && types[i] != itype ){
				continue;
			}
			if( jtype != 0 && types[j] != jtype ){
				continue;
			}
			
			if( r2 < rc2 ){
				neighs[i].push_back(j);
			}
		}
	}

}

void loop_bins_make_neighs_2d( const arr3f &x, py_int N, py_float rc,
                               const arr1i &ids, py_int periodic,
                               const py_float *xlo, const py_float *xhi, py_int Nx,
                               py_int Ny, py_int Nz, biguint Nbins,
                               biguint *bin_indices, std::list<py_int> *dom_bins,
                               std::list<py_int> *neighs, py_int itype, py_int jtype,
                               const arr1i &types )
{
	biguint loop_idx[9];
	assert( loop_idx && "Failed index space allocation! Use nsq!" );
	for( py_int i = 0; i < N; ++i ){
		fprintf(stdout, "At atom id %ld\n", ids[i]);
		
		loop_idx[0] = bin_indices[i];
		loop_idx[1] = shift_bin_index(loop_idx[0],  1,  0, 0, Nx, Ny, Nz, periodic);
		loop_idx[2] = shift_bin_index(loop_idx[0], -1,  0, 0, Nx, Ny, Nz, periodic);
		loop_idx[3] = shift_bin_index(loop_idx[0],  1,  1, 0, Nx, Ny, Nz, periodic);
		loop_idx[4] = shift_bin_index(loop_idx[0], -1,  1, 0, Nx, Ny, Nz, periodic);
		loop_idx[5] = shift_bin_index(loop_idx[0],  0,  1, 0, Nx, Ny, Nz, periodic);
		loop_idx[6] = shift_bin_index(loop_idx[0],  0, -1, 0, Nx, Ny, Nz, periodic);
		loop_idx[7] = shift_bin_index(loop_idx[0],  1, -1, 0, Nx, Ny, Nz, periodic);
		loop_idx[8] = shift_bin_index(loop_idx[0], -1, -1, 0, Nx, Ny, Nz, periodic);

		for( py_int bini = 0; bini < 9; ++bini ){
			
			py_int bin_idx = loop_idx[bini];
			if( (bin_idx < 0) || (bin_idx >= Nbins) ){
				if( !periodic ){
					continue;
				}else{
					fprintf(stderr,"Illegal bin!\n");
					abort();
				}
			}
			const std::list<py_int> &bin = dom_bins[ bin_idx ];
			add_neighs_from_bin( i, x, rc, ids, periodic, xlo, xhi,
			                     2, bin, neighs, itype, jtype, types );
		}
	}
}
	

void loop_bins_make_neighs_3d( const arr3f &x, py_int N, py_float rc,
                               const arr1i &ids, py_int periodic,
                               const py_float *xlo, const py_float *xhi,
                               py_int Nx, py_int Ny, py_int Nz, biguint Nbins,
                               biguint *bin_indices, std::list<py_int> *dom_bins,
                               std::list<py_int> *neighs, py_int itype, py_int jtype,
                               const arr1i &types  )
{
	biguint loop_idx[27];
	assert( loop_idx && "Failed index space allocation! Use nsq!" );
	for( py_int i = 0; i < N; ++i ){
		loop_idx[ 0] = bin_indices[i];
		loop_idx[ 1] = shift_bin_index(loop_idx[0],  1,  0, 0, Nx, Ny, Nz, periodic);
		loop_idx[ 2] = shift_bin_index(loop_idx[0], -1,  0, 0, Nx, Ny, Nz, periodic);
		loop_idx[ 3] = shift_bin_index(loop_idx[0],  1,  1, 0, Nx, Ny, Nz, periodic);
		loop_idx[ 4] = shift_bin_index(loop_idx[0], -1,  1, 0, Nx, Ny, Nz, periodic);
		loop_idx[ 5] = shift_bin_index(loop_idx[0],  0,  1, 0, Nx, Ny, Nz, periodic);
		loop_idx[ 6] = shift_bin_index(loop_idx[0],  0, -1, 0, Nx, Ny, Nz, periodic);
		loop_idx[ 7] = shift_bin_index(loop_idx[0],  1, -1, 0, Nx, Ny, Nz, periodic);
		loop_idx[ 8] = shift_bin_index(loop_idx[0], -1, -1, 0, Nx, Ny, Nz, periodic);

		loop_idx[ 9] = shift_bin_index(loop_idx[0],  0,  0,  1, Nx, Ny, Nz, periodic);
		loop_idx[10] = shift_bin_index(loop_idx[0],  1,  0,  1, Nx, Ny, Nz, periodic);
		loop_idx[11] = shift_bin_index(loop_idx[0], -1,  0,  1, Nx, Ny, Nz, periodic);
		loop_idx[12] = shift_bin_index(loop_idx[0],  1,  1,  1, Nx, Ny, Nz, periodic);
		loop_idx[13] = shift_bin_index(loop_idx[0], -1,  1,  1, Nx, Ny, Nz, periodic);
		loop_idx[14] = shift_bin_index(loop_idx[0],  0,  1,  1, Nx, Ny, Nz, periodic);
		loop_idx[15] = shift_bin_index(loop_idx[0],  0, -1,  1, Nx, Ny, Nz, periodic);
		loop_idx[16] = shift_bin_index(loop_idx[0],  1, -1,  1, Nx, Ny, Nz, periodic);
		loop_idx[17] = shift_bin_index(loop_idx[0], -1, -1,  1, Nx, Ny, Nz, periodic);		

		loop_idx[18] = shift_bin_index(loop_idx[0],  0,  0, -1, Nx, Ny, Nz, periodic);
		loop_idx[19] = shift_bin_index(loop_idx[0],  1,  0, -1, Nx, Ny, Nz, periodic);
		loop_idx[20] = shift_bin_index(loop_idx[0], -1,  0, -1, Nx, Ny, Nz, periodic);
		loop_idx[21] = shift_bin_index(loop_idx[0],  1,  1, -1, Nx, Ny, Nz, periodic);
		loop_idx[22] = shift_bin_index(loop_idx[0], -1,  1, -1, Nx, Ny, Nz, periodic);
		loop_idx[23] = shift_bin_index(loop_idx[0],  0,  1, -1, Nx, Ny, Nz, periodic);
		loop_idx[24] = shift_bin_index(loop_idx[0],  0, -1, -1, Nx, Ny, Nz, periodic);
		loop_idx[25] = shift_bin_index(loop_idx[0],  1, -1, -1, Nx, Ny, Nz, periodic);
		loop_idx[26] = shift_bin_index(loop_idx[0], -1, -1, -1, Nx, Ny, Nz, periodic);		

		for( py_int bini = 0; bini < 27; ++bini ){
			
			py_int bin_idx = loop_idx[bini];
			if( (bin_idx < 0) || (bin_idx >= Nbins) ){
				if( !periodic ){
					continue;
				}else{
					fprintf(stderr,"Illegal bin!\n");
					abort();
				}
			}
			const std::list<py_int> &bin = dom_bins[ bin_idx ];
			add_neighs_from_bin( i, x, rc, ids, periodic, xlo, xhi,
			                     3, bin, neighs, itype, jtype, types );
		}
		
	}
}

void test_bin_shifting()
{
	/* Test on a few easy numbers. */
	biguint indices[4];
	py_int binx[4], biny[4], binz[4], Nx[4], Ny[4], Nz[4], nbins[4];

	binx[0] = 1, biny[0] = 2, binz[0] = 4, Nx[0] = 2, Ny[0] = 4, Nz[0] = 8;
	binx[1] = 3, biny[1] = 3, binz[1] = 1, Nx[1] = 2, Ny[1] = 4, Nz[1] = 3;
	binx[2] = 1, biny[2] = 1, binz[2] = 2, Nx[2] = 2, Ny[2] = 2, Nz[2] = 4;
	binx[3] = 1, biny[3] = 1, binz[3] = 2, Nx[3] = 2, Ny[3] = 4, Nz[3] = 8;

	for( py_int i = 0; i < 4; ++i ){
		nbins[i] = Nx[i]*Ny[i]*Nz[i];
		indices[i] = binx[i] + biny[i]*Nx[i] + binz[i]*Nx[i]*Ny[i];
		biguint i1 = shift_bin_index(indices[i], 0, 0, 0,
		                              Nx[i], Ny[i], Nz[i], 1 );
		
		assert( indices[i] == i1 && "wtf in test_bin_shifting" );
		i1 = shift_bin_index(indices[i], 1, 1, 1,
		                     Nx[i], Ny[i], Nz[i], 1 );
		assert( (i1 > 0) && (i1 < nbins[i]) &&
		        "WTF in test_bin_shifting" );

		biguint i2 = indices[i];
		i2 = shift_bin_index(indices[i],  1 , 1,  1, Nx[i], Ny[i], Nz[i], 1);
		i2 = shift_bin_index(i2,         -1, -1, -1, Nx[i], Ny[i], Nz[i], 1);
		assert( (i2 == indices[i]) && "wtf in test_bin_shifting" );
		
	}
}



void neighborize_dist_bin_impl( const arr3f &x, py_int N, const arr1i &ids,
                                const arr1i &types, py_float rc, py_int periodic,
                                const py_float *xlo, const py_float *xhi, py_int dims,
                                std::list<py_int> *neighs, py_int itype, py_int jtype )
{
	/* This is more tricky. You need to set up a huge amount of bins... */
	double pad = 0.3;
	double bin_r = rc + pad;
	py_int Nx = ceil( (xhi[0] - xlo[0])/bin_r );
	py_int Ny = ceil( (xhi[1] - xlo[1])/bin_r );
	py_int Nz = ceil( (xhi[2] - xlo[2])/bin_r );
	if( Nx < 3 ){
		fprintf(stderr,"!! Bins in x-direction too thin, padding...\n");
		Nx = 3;
	}
	if( Ny < 3 ){
		fprintf(stderr,"!! Bins in y-direction too thin, padding...\n");
		Ny = 3;
	}
	if( (Nz < 3) && (dims == 3) ){
		fprintf(stderr,"!! Bins in z-direction too thin, padding...\n");
		Nx = 3;
	}
	biguint Nbins = Nx*Ny*Nz;

	py_int max_id = *std::max_element(ids.begin(), ids.end());
	biguint *bin_indices = new biguint[N];
	

	/* Now you need a list for each bin... */
	std::list<py_int> *dom_bins = new std::list<py_int>[Nbins];
	assert(dom_bins && "Failed to space for bins, use neighborize nsq!");

	/* Determine for each particle the bin it needs, add it to that list. */
	for( py_int i = 0; i < N; ++i ){
		const py_float *xi = x[i];
		py_int xbin = ( xi[0] - xlo[0] ) / bin_r;
		py_int ybin = ( xi[1] - xlo[1] ) / bin_r;
		py_int zbin = ( xi[2] - xlo[2] ) / bin_r;
		assert( (xbin < Nx) && (ybin < Ny) && (zbin < Nz) &&
		        "Incorrect bin obtained!" );

		biguint bini = xbin + Nx*ybin + Nx*Ny*zbin;
		dom_bins[bini].push_back(i);
		bin_indices[i] = bini; /* Store for easy lookup */
	}

	/*
	  Here is a reminder about the current state of affairs:
	   - dom_bins contains a bunch of lists in an array the size of Nbins.
	   - bin_indices contains the bin index of each atom in the arrays.
	     Note that these are INDICES, not IDS!
	   - neighs has room for N lists of py_ints, which are to be filled
	     with atom INDICESs, not IDs!
	*/
	  
	
	if( dims == 2 ){
		loop_bins_make_neighs_2d( x, N, rc,  ids, periodic, xlo, xhi, Nx, Ny,
		                          Nz, Nbins, bin_indices, dom_bins, neighs, itype, jtype, types );
	}else{
		loop_bins_make_neighs_3d( x, N, rc, ids, periodic, xlo, xhi, Nx, Ny,
		                          Nz, Nbins, bin_indices, dom_bins, neighs, itype, jtype, types );
	}
	

	delete [] bin_indices;
	delete [] dom_bins;
}





void neighborize_dist_bin( const arr3f &x, py_int N, const arr1i &ids,
                           const arr1i &types, py_float rc, py_int periodic,
                           const py_float *xlo, const py_float *xhi, py_int dims,
                           std::list<py_int> *neighs, py_int itype, py_int jtype )
{
	test_bin_shifting();

	neighborize_dist_bin_impl( x, N, ids, types, rc, periodic,
	                           xlo, xhi, dims, neighs, itype, jtype );
}



void neighborize_block( const block_data &b,
                        std::vector<std::list<py_int> > &neighs )
{
	double *x_ = new double[3*b.N];
	py_int *i_ = new py_int[b.N];
	py_int *t_ = new py_int[b.N];
	
	arr3f x(x_, b.N);
	arr1i ids(i_,b.N);
	arr1i types(t_,b.N);
	for( int i = 0; i < b.N; ++i ){
		x[i][0]  = b.x[i][0];
		x[i][1]  = b.x[i][1];
		x[i][2]  = b.x[i][2];
		ids[i]   = b.ids[i];
		types[i] = b.types[i];
	}

	py_int method = 0;
	if( b.N > 1500 ){
		method = DIST_BIN;
	}
	
	neighs.resize( b.N );
	neighborize_impl( x, b.N, ids, types, 1.3, 0,
	                  b.xlo, b.xhi, 3, method, neighs.data(), 0, 0 );
}

#include <cmath>
#include <memory>
#include <algorithm>

#include "ensemble_generator.h"


// This software NEEDS a LAMMPS library that was compiled with
// the user-manifold package.
#ifdef HAVE_LIB_LAMMPS
#define USE_PHONY_LAMMPS
#include "manifold_factory.h"
#include <mpi.h>
#endif

#include "domain.h"
#include "random_pcg.h"


#include "my_output.hpp"
static my_ostream my_out( std::cout );


extern "C" {
void ensemble_generate_manifold( void *x, py_int N, py_int *ids,
                                 py_int *types, py_int periodic,
                                 py_float *xlo, py_float *xhi,
                                 const char *manifold, const char *man_args,
                                 py_float rc, py_int seed )
{
#ifndef HAVE_LIB_LAMMPS
	std::cerr << "ensemble_generate_manifold can only be used if compiled "
	          << "with LAMMPS support!\n";
#else
	arr1i t_ids(ids,N);
	arr1i t_types(ids,N);
	arr3f t_x(x,N);

	ensemble_generate_impl( t_x, N, t_ids, t_types, periodic,
	                        xlo, xhi, manifold, man_args, rc, seed );
#endif	
}

}

#ifdef HAVE_LIB_LAMMPS
int line_search( double x[3], int id, double tol, std::size_t maxit,
                 LAMMPS_NS::user_manifold::manifold *ptr_m,
                 py_float *xlo, py_float *xhi, std::size_t &iters )
{
	double res = std::abs( ptr_m->g(x) );	
	
	for( iters = 0; iters < maxit; ++iters ){
		// The only thing I can think of is doing a line search:
		// g(x + alpha*n) = 0 -->
		// g(x + alpha*n) ~= g(x) + alpha * n * dg/dx =
		// g(x) + alpha * n^2 = 0, so alpha = -g(x) / n^2
		// Therefore, as an update, do
		// x = x - g(x) * n / n^2.
		double n[3];
		ptr_m->n(x,n);
		double g = ptr_m->g(x);
		double n2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
		double c = g / n2;
		double l = 1.0 / sqrt(1.0 + res*res);
		x[0] -= l*c * n[0];
		x[1] -= l*c * n[1];
		x[2] -= l*c * n[2];

		res = std::abs( g );
		if( res < tol ){
			break;
		}
	}
	if( res >= tol ){
		return 1;
	}else{
		// check if position was valid:
		return (x[0] < xlo[0]) || (x[0] > xhi[0]) || (x[1] < xlo[1]) ||
			(x[1] > xhi[1]) || (x[2] < xlo[2]) || (x[2] > xhi[2]);
	}
}





void parse_man_args( const char *man_args, int &man_narg, char **&man_arg_arr )
{
	// Read through string and count separators:
	std::string s( man_args );
	man_narg = std::count(s.begin(), s.end(), ' ');

	man_arg_arr = new char*[man_narg];

	if( !man_arg_arr ) return;

	int L = strlen(man_args);
	
	int start, end = 0;
	start = end = 0;
	int i = 0;
	int curr = 0;
	while( i < L ){
		std::cerr << "At character '" << man_args[i] << "'.\n";
		if( man_args[i] == ' ' ){
			end = i;
			int size_needed = end - start;
			
			std::cerr << "Allocating " << size_needed << " chars...\n";
			man_arg_arr[curr] = new char[size_needed+1];
			if( !man_arg_arr[curr] ){
				// Problem.
				std::cerr << "WTF?!\n";
				
			}
			

			strncpy( man_arg_arr[curr], man_args + start, size_needed );
			man_arg_arr[curr][size_needed] = '\0';
			std::cerr << "Arg " << curr << " is now " << man_arg_arr[curr] << ".\n";
			++curr;
			start = end+1;
			// ++i; // To skip over the space
		}
		++i;
	}

	// Test:
	std::cerr << "Found " << man_narg << " manifold args:\n";
	for( int i = 0; i < man_narg; ++i ){
		std::cerr << i+1 << ": " << man_arg_arr[i] << "\n";
	}
	std::cerr << "String passed from python was " << man_args << ".\n";
	
}

void ensemble_generate_impl( arr3f &x, py_int N, arr1i &ids, arr1i &types,
                             py_int periodic, py_float *xlo, py_float *xhi,
                             const char *mname, const char *man_args, py_float rc,
                             py_int seed )
{
	using manifold = LAMMPS_NS::user_manifold::manifold;
	using man_ptr = std::unique_ptr<manifold>;
	using LAMMPS_NS::user_manifold::create_manifold;

	

	char **man_arg_arr = nullptr;
	int man_narg = 0;

	parse_man_args( man_args, man_narg, man_arg_arr );

	man_ptr ptr_m;

	std::cerr << "Attempting to create manifold of type " << mname << ".\n";
	LAMMPS_NS::user_manifold::manifold *manifold_ptr = 
		create_manifold( mname, NULL, man_narg, man_arg_arr );
	std::cerr << "Manifold_ptr is " << manifold_ptr << "...\n";
	ptr_m = man_ptr( manifold_ptr );

	if( !ptr_m ){
		std::cerr << "Allocation of manifold " << mname << " failed!\n";
		return;
	}else{
		std::cerr << "Succesfully allocated manifold of type " << mname << ".\n";
	}

	// Init the manifold now.
	ptr_m->params = new double[ptr_m->nparams()];
	for( int i = 0; i < man_narg; ++i ){
		ptr_m->params[i] = std::stof( man_arg_arr[i] );
	}
	ptr_m->post_param_init();
	

	// Clean up the manifold args:
	std::cerr << "Cleaning up man_arg_arr... \n";
	for( int i = 0; i < man_narg; ++i ){
		std::cerr << "Deleting arg " << i << ": ";
		std::cerr << man_arg_arr[i] << "\n";
		delete [] man_arg_arr[i];
	}
	delete [] man_arg_arr;

	std::cerr << "Generating " << N << " points on " << mname << " in "
		"domain [ " << xlo[0] << ", " << xhi[0] << " ] x [ " << xlo[1]
	          << ", " << xhi[1] << " ] x [ " << xlo[2] << ", " << xhi[2] << " ].\n";
	py_int i = 0;
	RanPCG ran_pcg(seed);
	double rc2 = rc*rc;
	while( i < N ){
		double xi[3];
		xi[0] = ran_pcg.uniform_lo_hi(xlo[0], xhi[0]);
		xi[1] = ran_pcg.uniform_lo_hi(xlo[1], xhi[1]);
		xi[2] = ran_pcg.uniform_lo_hi(xlo[2], xhi[2]);
		std::cerr << i << ": Attempting to start with point ( " << xi[0]
		          << ", " << xi[1] << ", " << xi[2] << " )...";


		std::size_t iters = 0;
		int status = line_search( xi, i+1, 1e-8, 1000, ptr_m.get(), xlo, xhi, iters );
		if( status ){
			std::cerr << "Something went wrong!\n";
			continue;
		}
		// Check if distance is OK (slow though, be careful!)
		bool reject = false;
		for( py_int j = 0; j < i; ++j ){
			double r[3];
			distance_wrap( r, xi, x[j], xlo, xhi, periodic );
			double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
			if( r2 < rc2 ){
				// Reject!
				reject = true;
				std::cerr << " Rejected!\n";
				break;
			}
		}
		if( !reject ){
			std::cerr << " Accepted!\n";
			if( i && (i % 10000) == 0 ){
				std::cerr << "At point " << i << "...\n";
			}
			// Push point to x:
			x[i][0]  = xi[0];
			x[i][1]  = xi[1];
			x[i][2]  = xi[2];
			ids[i]   = i+1;
		        
			++i;
		}
	}

	
	
}
#endif // HAVE_LIB_LAMMPS

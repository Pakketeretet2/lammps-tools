#include "line_tension.h"
#include "dump_reader.h"
#include "neighborize.h"

double estimate_line_tension( const block_data &edge_block,
                              const std::vector<std::list<py_int> > &neighs,
                              std::size_t nn_expect, double F_per_particle )
{
	double avg_gamma = 0.0;
	double count = 0.0;
	for( py_int i = 0; i < edge_block.N; ++i ){
		int nn = neighs[i].size();
		if( nn == 0 ) continue;
		
		int deficit = nn_expect - nn;
		if( deficit != 0 ){
			double def_fac = static_cast<double>( deficit ) / nn;
			avg_gamma += def_fac * F_per_particle;
			count += 1.0;
			/*
			std::cerr << "deficit = " << deficit
			          << ", def_fac = " << def_fac
			          << ", nn = " << nn << ", n_expect = "
			          << nn_expect << "\n";
			*/
		}
	}
	avg_gamma /= count;
	/*
	std::cerr << "Counted " << count << " edge particles, gamma ~= "
	          << avg_gamma << "\n";

	if( avg_gamma < -100 ){
		std::terminate();
	}
	*/
	return avg_gamma;
}


extern "C" {
py_float get_line_tension( void *x, py_int N, py_int *ids,
                           py_int *types, py_int periodic,
                           py_float *xlo, py_float *xhi, py_int dims,
                           py_int tstep, const char *boxline,
                           py_int nn_expect, py_float F_per_particle )
{
	block_data b;
	block_data_from_foreign( x, N, ids, types, periodic, xlo, xhi, dims,
	                         tstep, boxline, b );
	std::vector<std::list<py_int> > neighs;
	neighborize_block( b, neighs );
	py_float gamma = estimate_line_tension( b, neighs,
	                                        nn_expect,
	                                        F_per_particle );
	return gamma;
	
}

} // extern "C"

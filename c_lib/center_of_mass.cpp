#include "center_of_mass.h"

#include <vector>
#include <cmath>




void com_periodic( const arr3f &x,  const arr1i &ids, const arr1i &types,
                   const arr1f &mass, double xlo, double xhi, int dim,
                   const arr1i &groups, arr3f &com )
{
	for( py_int i = 0; i < com.size(); ++i ){
		com[i][dim] = 0.0;
	}

	py_int Ngroups = com.size() - 1;
	
	std::vector<double> xi_avg(Ngroups+1);
	std::vector<double> zeta_avg(Ngroups+1);
	std::vector<double> M(Ngroups+1);
	double two_pi = 2.0*math_const::pi;

	double Lx = xhi - xlo;

	for( py_int i = 0; i < x.size(); ++i ){
		py_int id = ids[i];
		py_int group_id = groups[i];

		if( group_id == 0 ) continue; // Use 0 to flag do not use.

		double m  = mass[i];
		double xx = (x[i][dim] - xlo) / Lx;

		int group_idx = group_id;
		
		// This hack is required because atoms might be slightly
		// outside of the box due to infrequent rewrapping.
		if( xx < 0.0 ){
			xx += 1.0;
		}else if( xx > 1.0 ){
			xx -= 1.0;
		}
		M[group_idx] += m;
		double theta = two_pi*xx;
		
		xi_avg[group_idx]   += m * std::cos( theta );
		zeta_avg[group_idx] += m * std::sin( theta );

	}

	for( int i = 1; i <= Ngroups; ++i ){
		xi_avg[i]   /= M[i];
		zeta_avg[i] /= M[i];
		double theta_avg = std::atan2( -zeta_avg[i], -xi_avg[i] );
		
		theta_avg += math_const::pi;
		com[i][dim] = theta_avg * Lx / two_pi + xlo;

	}
}


void com_nonperiodic( const arr3f &x,  const arr1i &ids, const arr1i &types,
                      const arr1f &mass, double xlo, double xhi, int dim,
                      const arr1i &groups, arr3f &com )
{
	for( py_int i = 0; i < com.size(); ++i ){
		com[i][dim] = 0.0;
	}

	py_int Ngroups = com.size() - 1;
	std::vector<double> x_avg(Ngroups+1);
	std::vector<double> M(Ngroups+1);
	
	double two_pi = 2.0*math_const::pi;

	double Lx = xhi - xlo;
	
	for( py_int i = 0; i < x.size(); ++i ){
		py_int id = ids[i];
		py_int group_id = groups[i];

		if( group_id == 0 ) continue; // Use 0 to flag do not use.

		double m  = mass[i];
		int group_idx = group_id;

		M[group_idx]     += m;
		x_avg[group_idx] += m*x[i][dim];
	}

	for( int i = 1; i <= Ngroups; ++i ){
		x_avg[i]   /= M[i];
		com[i][dim] = x_avg[i];
	}
}

extern "C" {

void center_of_mass( void *px, py_int N, py_int *pids, py_int *ptypes,
                     py_int Ngroups, py_int *pgroups, py_float *pmass,
                     py_float *xlo, py_float *xhi, py_int periodic,
                     py_int dims, void *pcom  )
{
	// Assumes com is already allocated.
	arr3f x( px, N );
	arr1i ids( pids, N );
	arr1i types( pids, N );
	arr1i groups( pgroups, N );
	arr1f mass( pmass, N );
	
	arr3f com( pcom, Ngroups + 1 );

	int dim = 0;
	if( periodic & 1 ){
		com_periodic( x, ids, types, mass, xlo[dim], xhi[dim],
		              dim, groups, com );
	}else{
		com_nonperiodic( x, ids, types, mass, xlo[dim], xhi[dim],
		                 dim, groups, com );		
	}

	dim = 1;
	if( periodic & 2 ){
		com_periodic( x, ids, types, mass, xlo[dim], xhi[dim],
		              dim, groups, com );
	}else{
		com_nonperiodic( x, ids, types, mass, xlo[dim], xhi[dim],
		                 dim, groups, com );		
	}

	dim = 2;
	if( periodic & 4 ){
		com_periodic( x, ids, types, mass, xlo[dim], xhi[dim],
		              dim, groups, com );
	}else{
		com_nonperiodic( x, ids, types, mass, xlo[dim], xhi[dim],
		                 dim, groups, com );		
	}
	
	
}

} // extern "C"

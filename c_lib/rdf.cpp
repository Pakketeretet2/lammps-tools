#include "rdf.h"
#include "neighborize.h"
#include "id_map.h"
#include "domain.h"
#include "my_output.hpp"

#include <algorithm>
#include <cmath>
#include <list>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <fstream>


static my_ostream my_out( std::cerr );

void compute_rdf_impl( const arr3f &x, py_int N, const arr1i &ids,
                       const arr1i &types, py_float x0, py_float x1,
                       py_int nbins, py_int itype, py_int jtype,
                       py_float *xlo, py_float *xhi,
                       py_int periodic, py_int dim, std::list<py_int> *neighs,
                       arr1f &ardf, arr1f &acoord )
{
	double dr  = (x1 - x0) / ( nbins - 1 );
	arr1f &rdf = ardf;
	arr1f &coord = acoord;
	
	for( py_int i = 0; i < nbins; ++i ){
		rdf[i]   = 0.0;
		coord[i] = 0.0;
	}

	// Count the total possible number of interactions:
	double count_i = 0.0, count_j = 0.0;
	for( py_int i = 0; i < N; ++i ){
		count_i += ( (types[i] == itype) + ( itype == 0) );
		count_j += ( (types[i] == jtype) + ( jtype == 0) );
	}
	
	std::ofstream out("test_rdf.dat");
	my_out << "Computing rdf for " << N << " particles over "
	          << nbins << " bins between types " << itype
	          << " (Ni = " << count_i << " ) and " << jtype << " (Nj = "
	          << count_j << " )...\n";
	my_out << "Binning between " << x0 << " and "
	         << x1 << " with a resolution " << dr << "\n";
	my_out << "dimensions = " << dim << " box is [ " << xlo[0] << ", "
	          << xhi[0] << " ] x [ " << xlo[1] << ", " << xhi[1]
	          << " ] x [ " << xlo[2] << ", " << xhi[2] << " ].\n";
	id_map id_to_idx(ids);


	py_int adds = 0;
	double Lx = xhi[0] - xlo[0];
	double Ly = xhi[1] - xlo[1];
	double Lz = xhi[2] - xlo[2];
	my_out << "Box lenghs: " << Lx << " by " << Ly << " by " << Lz << "\n";
	

	// Loop over the neighbor list of each particle.
	for( py_int i = 0; i < N; ++i ){
		// if( types[i] != itype ) continue;
		if( itype && (types[i] != itype) ) continue;
		
		py_int id = ids[i];
		const std::list<py_int> &curr = neighs[i];
		
		for ( py_int j : curr ){
			// if( jd <= id ) continue;
			if( ids[j] >= ids[i] ) continue;
			if( jtype && (types[j] != jtype) ) continue;

			const py_float *xi = x[i];
			const py_float *xj = x[j];
			py_float r[3];

			distance_wrap( r, xi, xj, xlo, xhi );

			double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];

			// determine bin
			double rr = sqrt(r2);
			py_int bini = (rr - x0) / dr;
			if( bini < 0 || bini >= nbins ) continue;

			rdf[bini] += 2.0;
			++adds;
		}
	}
	
	
	my_out << "Binned " << adds << " interactions.\n";

	//  Now normalize RDF and create CDF.
	// Most code borrowed from LAMMPS' compute_rdf.

	double V = Lx*Ly;
	if( dim != 2 ) V *= Lz;

	double factor;
	if( dim != 2 ){
		factor = 4.0*math_const::pi / (3.0 * V);
	}else{
		factor = math_const::pi / V;	
	}
	// factor *= count_j;
	coord[0] = 0.0;
	// V is already set, it is total volume for 3D or area for 2D.
	for( py_int bin = 0; bin < nbins; ++bin ){
		double binr  = x0 + dr*bin;
		double binrp = binr + dr;
		double nideal = factor*(binrp*binrp*binrp - binr*binr*binr)*count_j;

		rdf[bin] /= (nideal*count_i);

		if( bin == 0 ){
			coord[bin] = 0.0;
		}else{
			coord[bin] = coord[bin-1] + rdf[bin]*nideal;
		}
	}

	for( py_int i = 0; i < nbins; ++i ){
		double r = x0 + (i+0.5)*dr;
		out << i << " " << r << " " << rdf[i]
		    << " " << coord[i] << "\n";
	}
}


void compute_adf_impl( const arr3f &x, py_int N, const arr1i &ids,
                       const arr1i &types, py_int nbins, py_int itype, py_int jtype,
                       py_float R, std::list<py_int> *neighs, arr1f &aadf, arr1f &acoord )
{
	// First check if the atom positions are all on the sphere:
	const double R2 = R*R;
	bool mismatch_radius = false;
	for( py_int i = 0; (i < x.size()) && (!mismatch_radius); ++i ){
		const py_float *xi = x[i];
		double r2 = xi[0]*xi[0] + xi[1]*xi[1] + xi[2]*xi[2];

		if( fabs( r2 - R2 ) > 1e-4 ){
			fprintf(stderr, "WARNING! Position of atom %ld does not "
			        " appear to not be on sphere! |x| = %f"
			        " Are you sure data is correct?\n", ids[i], sqrt(r2));
			mismatch_radius = true;
		}
	}

	auto dist_to_angle = [R2,mismatch_radius] (const py_float *x1, const py_float *x2)
		{
			const double dot = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
			if( mismatch_radius ){
				return acos( dot / R2 );
			}else{
				const double r12 = x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2];
				const double r22 = x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2];
				
				const double R2_alt = sqrt(r12)*sqrt(r22);
				return acos( dot / R2_alt );
			}
		};

	// Determine bin width in terms of angles:
	double dtheta = math_const::pi / (nbins - 1 );
	py_float weight = 1.0 / N;

	// Put all bins to 0:
	for( py_int bin = 0; bin < nbins; ++bin ){
		acoord[bin] = aadf[bin] = 0.0;
	}
	
	for( py_int i = 0; i < x.size(); ++i ){
		const py_float *xi = x[i];
		if( i == x.size() ){
			std::cerr << "Error! index problem!\n";
			return;
		}

		for( py_int j : neighs[i] ){
			if( j == x.size() ){
				std::cerr << "Error! index problem!\n";
				return;
			}
			const py_float *xj = x[j];
			double theta_ij = dist_to_angle( xi, xj );
			py_int bin = static_cast<py_int>( theta_ij / dtheta );
			aadf[bin] += weight;
		}
	}

	// First determine coordination numbers:
	acoord[0] = 0.0;
	for( py_int bin = 1; bin < nbins; ++bin ){
		acoord[bin] = aadf[bin] + acoord[bin-1];
	}

	// Now normalize the ADF:
	if( aadf[0] != 0.0 ){
		// This is most likely an error.
		fprintf(stderr, "WARNING! adf at theta = 0 was not 0 but %f!\n",
			aadf[0]);
	}
	double As   = 4 * math_const::pi * R * R;
	double rho  = x.size() / As;
	double norm = 1.0 / ( rho * 2*math_const::pi * R * R * dtheta );
	for( py_int bin = 1; bin < nbins; ++bin ){
		py_float theta0 = (bin-1)*dtheta;
		py_float theta1 = bin*dtheta;
		py_float theta  = theta0 + 0.5*(theta1 - theta0);
		
		aadf[bin] *= norm / sin(theta);
	}
	std::ofstream out_adf( "adf_test.dat" );
	for( py_int i = 0; i < nbins; ++i ){
		out_adf << i << " " << i*dtheta << " " << aadf[i] << " " << acoord[i] << "\n";
	}
}




extern "C" {


void compute_rdf( void *px, py_int N, py_int *pids, py_int *ptypes,
                  py_float x0, py_float x1, py_int nbins, py_int itype,
                  py_int jtype, py_float *xlo, py_float *xhi, py_int periodic,
                  py_int dim, py_int method, py_float *prdf, py_float *pcoord )
{
	arr1i ids(pids,N);
	arr1i types(ptypes,N);
	arr3f x(px,N);
	arr1f rdf  ( prdf,   nbins );
	arr1f coord( pcoord, nbins );

	// Make neigh list first:
	py_int max_id = *std::max_element( ids.begin(), ids.end() );
	std::list<py_int> *neighs = new std::list<py_int>[max_id+1];

	neighborize_impl( x, N, ids, types, x1+0.1, periodic,
	                  xlo, xhi, dim, method, neighs, 0, 0 );

	compute_rdf_impl( x, N, ids, types, x0, x1, nbins, itype, jtype,
	                  xlo, xhi, periodic, dim, neighs, rdf, coord );

	delete [] neighs;
}



void test_rdf()
{
	py_int N = 5;
	py_float xlo[3] = { 0, 0, 0 }, xhi[3] = { 6, 6, 6 };
	py_int nbins = 101;

	py_int iids[5] = {1,2,3,4,5}, ttypes[5] = {1,1,1,1,1};

	arr1f rdf(NULL, nbins), coord(NULL, nbins);
	arr1i ids(iids,5), types(ttypes,5);

	py_float *rdf_data   = new py_float[nbins];
	py_float *coord_data = new py_float[nbins];
	py_float **x_data    = new py_float*[N];
	for( py_int i = 0; i < N; ++i ){
		x_data[i] = new py_float[3];
	}
	
	for( uint i = 0; i < nbins; ++i ){
		rdf_data[i] = coord_data[i] = 0.0;
	}
	rdf.set_data(rdf_data);
	coord.set_data(coord_data);

	assert( rdf.raw_data()   && "RDF data pointer was NULL!" );
	assert( coord.raw_data() && "coord data pointer was NULL!" );

	arr3f x(x_data,N);


	x[0][0] = 1.0;
	x[0][1] = 0.0;
	x[0][2] = 0.0;

	x[1][0] = 1.0;
	x[1][1] = 0.0;
	x[1][2] = 0.0;

	x[2][0] = 0.0;
	x[2][1] = 1.0;
	x[2][2] = 0.0;
	
	x[3][0] = 5.0;
	x[3][1] = 0.0;
	x[3][2] = 0.0;

	x[4][0] = 1.0;
	x[4][1] = 0.0;
	x[4][2] = 5.0;
	



	py_float r0 = 0, r1 = 5.0;
	py_float dr = (r1 - r0)/(nbins-1);

	py_int max_id = *std::max_element(ids.begin(), ids.end());
	std::list<py_int> *neighs = new std::list<py_int>[max_id+1];

	neighborize_impl( x, N, ids, types, r1+0.1, 1,
	                  xlo, xhi, 3, DIST_BIN, neighs, 0, 0 );
	
	fprintf( stdout, "Neighbors:\n");
	for( py_int i = 0; i < max_id+1; ++i ){
		my_out << "neighs[ " << i << " ]: ";
		for ( std::list<py_int>::const_iterator j = neighs[i].begin();
		      j != neighs[i].end(); ++j ){
			my_out << *j;
			std::list<py_int>::const_iterator next = j;
			next++;
			if( next != neighs[i].end() ){
				my_out << " --> ";
			}
		}
		my_out << "\n";
	}
	
	compute_rdf_impl( x, N, ids, types, r0, r1,
	                  nbins, 1, 1, xlo, xhi, 1, 3, neighs,
	                  rdf, coord );

	fprintf(stdout, "RDF and coordination:\n");
	for( py_int i = 0; i < nbins; ++i ){
		py_float r = r0 + i*dr;
		fprintf(stdout,"%ld  %f   %f   %f\n",i, r, rdf[i], coord[i]);
	}

	delete [] neighs;
	delete [] rdf_data;
	delete [] coord_data;
	for( py_int i = 0; i < N; ++i ){
		delete [] x_data[i];
	}
	delete [] x_data;
}





void compute_adf( void *px, py_int N, py_int *pids, py_int *ptypes,
                  py_int nbins, py_int itype, py_int jtype, py_float R,
                  py_int method, py_float *padf, py_float *pcoord )
{
	if( method == DELAUNAY ){
		fprintf(stderr, "Cannot use Delaunay with spherical data!\n");
		return;
	}
	
	arr1i ids(pids,N);
	arr1i types(ptypes,N);
	arr3f x(px,N);
	arr1f adf  ( padf,   nbins );
	arr1f coord( pcoord, nbins );

	// Make neigh list first:
	py_int max_id = *std::max_element( ids.begin(), ids.end() );
	std::list<py_int> *neighs = new std::list<py_int>[max_id+1];

	py_float *xlo = nullptr, *xhi = nullptr;
	
	if( method == DIST_BIN ){
		fprintf( stderr, "Using bins for ADF does not make sense, "
		         "using nsq instead!" );
		method = DIST_NSQ;
	}

	// 2R + 0.2 should be enough to include ALL neighbours.
	neighborize_impl( x, N, ids, types, 2*R + 0.2, 0,
	                  xlo, xhi, 3, method, neighs, 0, 0 );

	compute_adf_impl( x, N, ids, types, nbins, itype,
	                  jtype, R, neighs, adf, coord );
	
	delete [] neighs;	
}
	

}// extern "C"

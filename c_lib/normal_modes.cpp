#include "normal_modes.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <array>
#include <armadillo>

#include "id_map.h"
#include "dump_reader.h"

typedef std::vector<std::array<double,3> > coord_block;


void normal_mode_analysis( dump_reader &r,  py_int N, void *eigenvalues,
                           void *eigenvectors )
{
	std::list< coord_block > blocks;
	block_data b;
	
	while( r.next_block( b ) ){
		id_map im( b.ids );
		py_int Nparticles = b.N;
		
		coord_block block( Nparticles );
		for( int i = 0; i < Nparticles; ++i ){
			int id  = i + 1;
			int idx = im[id];

			block[i][0] = b.x[i][0];
			block[i][1] = b.x[i][1];
			block[i][2] = b.x[i][2];
		}

		blocks.push_back( block );
	}
	std::cerr << "Done reading dump file...\n";

	get_normal_modes( blocks, N, eigenvalues, eigenvectors );

}


void read_blocks_from_pipe( const char *pname,
                            std::list< coord_block > &bs )
{
	std::ifstream in(pname);
	std::string line;
	std::ofstream out("pipe_test2");
	
	/*
	  Expected format:
	  First line: Number of particles
	  N times 3 coordinates, ordered along id.
	*/

	int block_idx = 1;
	
	// Make sure you exhaust the pipe:
	while( in ){
		std::getline( in, line );
		int Nparticles = 0;
		std::stringstream ss(line);
		ss >> Nparticles;

		if( Nparticles <= 0 ) break;
		
		coord_block block(Nparticles);

		
		for( int i = 0; i < Nparticles; ++i ){
			std::getline( in, line );
			ss.str("");
			ss.clear();
			ss << line;
			ss  >> block[i][0] >> block[i][1] >> block[i][2];
		}
		block_idx++;
		bs.push_back( block );
	}
	std::cerr << "Grabbed " << bs.size() << " blocks in total.\n";
}


void average_positions( const std::list< coord_block > &block_data, py_int N,
                        std::vector<double> &X_avg )
{
	int Nparticles = block_data.begin()->size();
	for( int i = 0; i < N; ++i ) X_avg[i] = 0.0;

	std::cerr << "Averaging positions...\n";
	for( const coord_block &b : block_data ){
		for( int i = 0; i < Nparticles; ++i ){
			X_avg[3*i + 0] += b[i][0];
			X_avg[3*i + 1] += b[i][1];
			X_avg[3*i + 2] += b[i][2];
		}
	}
}


void get_normal_modes( const std::list< coord_block > &block_data, py_int N,
                       void *eigenvalues, void *eigenvectors,
                       std::vector<double> *X_avg_ptr )
{
	// Perform the normal mode analysis...
	int Nparticles = block_data.begin()->size();
	double c = 1.0 / static_cast<double>( block_data.size() );
	if( 3*Nparticles != N ){
		std::cerr << "Something fishy is going on...\n";
		std::cerr << "Nparticles = " << Nparticles << ", 3*Nparticles = "
		          << 3*Nparticles << " != N = " << N << "!\n";
		return;
	}

	std::vector<double> X_avg;
	std::cerr << "Averaging positions...\n";
	if( !X_avg_ptr ){
		X_avg.resize(N);
		X_avg_ptr = &X_avg;
		average_positions( block_data, N, *X_avg_ptr );
	}else{
		X_avg = *X_avg_ptr;
		average_positions( block_data, N, *X_avg_ptr );
	}
	
	// Extract covariance matrix:
	std::cerr << "Constructing covariance matrix...\n";
	arma::mat covar(N,N);
	
	for( int i = 0; i < N; ++i ){
		for( int j = 0; j < N; ++j ){
			covar(i, j) = 0.0;
		}
	}
	
	for( const coord_block &b : block_data ){
		for( int i = 0; i < Nparticles; ++i ){
			double dxi = b[i][0] - X_avg[3*i];
			double dyi = b[i][1] - X_avg[3*i+1];
			double dzi = b[i][2] - X_avg[3*i+2];

			for( int j = 0; j < Nparticles; ++j ){
				double dxj = b[j][0] - X_avg[3*j];
				double dyj = b[j][1] - X_avg[3*j+1];
				double dzj = b[j][2] - X_avg[3*j+2];

				covar(3*i  , 3*j  ) += dxi*dxj*c;
				covar(3*i  , 3*j+1) += dxi*dyj*c;
				covar(3*i  , 3*j+2) += dxi*dzj*c;

				covar(3*i+1, 3*j  ) += dyi*dxj*c;
				covar(3*i+1, 3*j+1) += dyi*dyj*c;
				covar(3*i+1, 3*j+2) += dyi*dzj*c;

				covar(3*i+2, 3*j  ) += dzi*dxj*c;
				covar(3*i+2, 3*j+1) += dzi*dyj*c;
				covar(3*i+2, 3*j+2) += dzi*dzj*c;
			}
		}
	}

	// Perform the decomposition:
	std::cerr << "Performing decomposition in normal modes...\n";
	const char *method = "std";
	arma::vec eigval(N);
	arma::mat eigvec(N,N);
	arma::eig_sym( eigval, eigvec, covar, method );
	
	std::cerr << "Writing " << N << " eigenvalues to " << eigenvalues
	          << " and eigenvectors to " << eigenvectors << "\n";

	double **eigenvecs = static_cast<double**>(eigenvectors);
	double *eigenvals  = static_cast<double*> (eigenvalues);

	for( int i = 0; i < N; ++i ){
		eigenvals[i] = eigval(i);
		for( int j = 0; j < N; ++j ){
			eigenvecs[i][j] = eigvec(i,j);
		}
	}
	std::cerr << "Done!\n";

}


void normal_mode_analysis( const char *pname, void *eigenvalues,
                           void *eigenvectors, py_int N  )
{

	std::list< coord_block > block_data;

	read_blocks_from_pipe( pname, block_data );
	
	get_normal_modes( block_data, N, eigenvalues, eigenvectors );
}


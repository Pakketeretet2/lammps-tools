#ifndef BLOCK_DATA_H
#define BLOCK_DATA_H

#include "types.h"

#include <algorithm>
#include <iosfwd>
#include <vector>
#include <list>
#include <string>



extern "C"{

struct dump_col
{
	std::string header;
	std::vector<double> data;

	void resize( py_int M )
	{
		data.resize( M );
	}
	
};

enum atom_styles {
	ATOMIC,
	MOLECULAR
};


struct block_data
{
	// Needs a custom assignment operator and copy constructor to work.
	py_float **x;
	py_float *x_;
	py_int *ids, *types;
	py_int *mol;
	py_int N, tstep;

	py_float xlo[3], xhi[3];
	py_int periodic;
	py_int atom_style;

	std::string boxline;
	std::vector<dump_col> other_cols;
	

	block_data();
	block_data( int N );
	~block_data();

	// Deliberately void to prevent assignment chaining.
	void operator=( const block_data &o );
	block_data( const block_data &o );
	void swap( block_data &o ) throw(); // For copy-and-swap idiom.
	
	void resize( int N );
	void init( int N );

private:
	void delete_members();
};


// void copy( block_data &b, const block_data &source );
	
}


void block_data_from_foreign( const void *x, py_int N, const py_int *ids,
                              const py_int *types, const py_int *mol,
                              py_int periodic,
                              const py_float *xlo, const py_float *xhi,
                              py_int dims, py_int tstep, const char *box_line,
                              block_data &b );

block_data block_data_from_data_file( const char *fname, int &status );



void print_block_data_lmp( const block_data &b, std::ostream &o );
void print_block_data_lmp( const block_data &b, const std::string &s );
void print_block_data_lmp( const block_data &b, const std::string &s,
                           std::ios_base::openmode mode );


template <typename container>
block_data filter_block( block_data b, const container &ids )
{
	block_data b_filter;
	py_int M = 0;
	b_filter = b;
	

	for( int i = 0; i < b.N; ++i ){
		auto idx = std::find( ids.begin(), ids.end(), b.ids[i] );
		if( idx != ids.end() ){
			b_filter.x[M][0] = b.x[i][0];
			b_filter.x[M][1] = b.x[i][1];
			b_filter.x[M][2] = b.x[i][2];

			b_filter.ids[M]  = b.ids[i];
			b_filter.types[M]  = b.types[i];
			if( b.mol ){
				b_filter.mol[M]  = b.mol[i];
			}
			++M;
		}
	}
	b_filter.resize( M );
	return b_filter;
}

template <typename container>
block_data filter_block_by_indices( block_data b, const container &indices )
{
	block_data b_filter;
	py_int M = 0;
	b_filter = b;
	

	for( py_int i : indices ){
		b_filter.x[M][0] = b.x[i][0];
		b_filter.x[M][1] = b.x[i][1];
		b_filter.x[M][2] = b.x[i][2];

		b_filter.ids[M]  = b.ids[i];
		b_filter.types[M]  = b.types[i];
		if( b.mol ){
			b_filter.mol[M]  = b.mol[i];
		}
		++M;
	}
	b_filter.resize( M );
	return b_filter;
}



// For some python functions that want to pass around a C++-type
// block_data directly to C++ functions.
extern "C" {

block_data *new_block_data();
void set_block_data( block_data *b,
                     py_int N, py_int t, py_float *x, py_int *ids,
                     py_int *types, py_int *mol,
                     py_float *xlo, py_float *xhi, py_int periodic,
                     const char* boxline );
void free_block_data( block_data * );

}





#endif // BLOCK_DATA_H

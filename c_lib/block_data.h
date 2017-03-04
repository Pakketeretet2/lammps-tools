#ifndef BLOCK_DATA_H
#define BLOCK_DATA_H

#include "types.h"

#include <iosfwd>
#include <vector>
#include <list>
#include <string>


extern "C"{
	
struct block_data
{
	py_float **x;
	py_float *x_;
	py_int *ids, *types;
	py_int *mol;
	py_int N, tstep;

	double xlo[3], xhi[3];
	py_int periodic;
	std::string boxline;

	enum atom_styles {
		ATOMIC,
		MOLECULAR
	};
	py_int atom_style;
};


block_data *new_block_data();
void        delete_block_data(block_data*);


}

void resize( block_data &b, int N );
void init( block_data &b, int N );


block_data filter_block( block_data b, const std::list<long int> &ids );
block_data filter_block( block_data b, const std::vector<long int> &ids );

block_data filter_block_by_indices( block_data b, const std::vector<long int> &indices );
block_data filter_block_by_indices( block_data b, const std::list<long int> &indices );

void block_data_from_foreign( void *x, py_int N, py_int *ids,
                              py_int *types, py_int periodic,
                              py_float *xlo, py_float *xhi, py_int dims,
                              py_int tstep, const char *box_line,
                              block_data &b );

block_data block_data_from_data_file( const char *fname, int &status );



void print_block_data_lmp( const block_data &b, std::ostream &o );
void print_block_data_lmp( const block_data &b, const std::string &s );
void print_block_data_lmp( const block_data &b, const std::string &s,
                           std::ios_base::openmode mode );


/*
class block_data
{
public:
	block_data();
	~block_data();

	block_data &operator=( const block_data &b );

	int alloc( int size );
	void clear();
	double memory();

	void init_empty( int size );
	void copy_meta( const block_data &b );

	void set_xyz_tags( const char *tx, const char *ty, const char *tz );
	void set_id_tag( const char *tid );
	void set_type_tag( const char *ttype );

	std::vector<std::array<double,3> > x;
	std::vector<long int> ids, types;

	std::string boxline;
	double xlo[3], xhi[3];
	int N;
	int tstep;

	std::vector<std::vector<double> > other_cols;
	std::vector<std::string> other_col_headers;

	int N_other_cols;
	int columns[5]; // Lists the index of id, type, x, y, z column, resp.
	std::vector<int> ocols;
	std::string idt, typet, xt, yt, zt;

	
	void print( std::ostream &out, std::ostream &err = std::cerr );
	void print( const std::string &out_name, std::ostream &err = std::cerr );
	void print_data( std::ostream &out, std::ostream &err = std::cerr );

	void set_column_indices( const std::string &line );		
	void word_to_array( const std::string &w, int idx, int i );
	int  get_column_index( const std::string &header );
	
	void append_other_col( std::string header, std::vector<double> data );
};
*/


#endif // BLOCK_DATA_H

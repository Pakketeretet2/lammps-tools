#ifndef DUMP_READER_H
#define DUMP_READER_H

#include "util.h"
#include <string>
#include <vector>
#include <fstream>
#include <array>

#define HAVE_BOOST_GZIP
#ifdef HAVE_BOOST_GZIP
#  include <boost/iostreams/filter/zlib.hpp>
#  include <boost/iostreams/filtering_stream.hpp>
#endif


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


class dump_reader
{
public:
	dump_reader( const std::string &fname );
	~dump_reader();

	bool next_block( block_data &block );
	bool last_block( block_data &block );

	operator bool() const
	{
		return !at_eof;
	}

	bool skip_blocks( int Nblocks );
	bool skip_block ( );

	std::size_t block_count();
	void rewind();
	
private:
	enum FILE_FORMATS {
		PLAIN = 0,
		GZIP  = 1
	};

	
	bool getline( std::string &line, int mode = 0 );
	void init_infile( const std::string &fname );

	bool at_eof;	
	int file_format;

	std::string infile_name;	
	std::ifstream infile;
	// For the gzip case:
#ifdef HAVE_BOOST_GZIP
	boost::iostreams::filtering_istream infile_filt;
#endif
	
};


bool get_dump_line( std::string &line, std::istream &in );
bool next_block_from_istream( block_data &b, std::istream &in );
bool last_block_from_istream( block_data &block, std::istream &in );
bool skip_block_from_istream( std::istream &in );

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



#endif // DUMP_READER_H

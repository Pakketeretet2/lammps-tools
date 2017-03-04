#include "dump_reader.h"
#include "id_map.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>

#include <boost/iostreams/filter/gzip.hpp>

struct dummy_stream
{
	template <typename T>
	dummy_stream &operator<<( const T &t )
	{
		return *this;
	}

	void flush(){}
};

struct debug_stream
{
	debug_stream( std::ostream &out ) : o(out){}
		
	template <typename T>
	debug_stream &operator<<( const T &t )
	{
		o << t;
		o.flush();
		return *this;
	}

	void flush()
	{
		o.flush();
	}
	
	std::ostream &o;
};

// static dummy_stream d_out;
//static std::ofstream d_out( "/dev/stderr" );
static debug_stream d_out( std::cerr );



dump_reader::dump_reader( const std::string &fname ) : at_eof(false),
                                                       infile_name(fname)
{
	init_infile(fname);
}



dump_reader::~dump_reader()
{
	
}

// Skips the next Nblocks blocks.
bool dump_reader::skip_blocks( int Nblocks )
{
	for( int i = 0; i < Nblocks; ++i ){
		bool success = skip_block();
		if( !success ){
			d_out << "Block skip failed, returning...\n";
			return false;
		}
	}
	return true;
}

bool dump_reader::skip_block()
{
	switch( file_format ){
#ifdef HAVE_BOOST_GZIP
		case GZIP:
			return skip_block_from_istream( infile_filt );
#endif
		case PLAIN:
		default:
			return skip_block_from_istream( infile );
	}	
	
	return true;
}



bool dump_reader::next_block( block_data &block )
{
	if( at_eof ) return false;

	switch( file_format ){
#ifdef HAVE_BOOST_GZIP
		case GZIP:
			return next_block_from_istream( block, infile_filt );
#endif
		case PLAIN:
		default:
			return next_block_from_istream( block, infile );
	}
}

bool dump_reader::last_block( block_data &block )
{
	if( at_eof ) return false;


	switch( file_format ){
#ifdef HAVE_BOOST_GZIP
		case GZIP:
			return last_block_from_istream( block, infile_filt );
#endif
		case PLAIN:
		default:
			return last_block_from_istream( block, infile );
	}
}



void dump_reader::init_infile( const std::string &fname )
{
	if( ends_with( fname, ".gz" ) ){
#ifdef HAVE_BOOST_GZIP
		infile = std::ifstream( fname, std::ios_base::in | std::ios_base::binary );
		file_format = GZIP;
		infile_filt.push(boost::iostreams::gzip_decompressor());
		infile_filt.push(infile);
#else
		assert(false && "BOOST GZIP LIBRARY NOT INSTALLED!");
#endif
	}else{
		file_format = PLAIN;
		infile = std::ifstream( fname, std::ios_base::in );
	}

	
}

bool dump_reader::getline( std::string &line, int mode )
{
	if( mode == 1 ){
		std::getline( std::cin, line );
		return false; // stuff;
	}
		
	switch( file_format ){
#ifdef HAVE_BOOST_GZIP
		case GZIP:
			if( std::getline( infile_filt, line ) ){
				return true;
			}else{
				return false;
			}
			break;
#endif
		case PLAIN:
		default:
			if( std::getline( infile, line ) ){
				return true;
			}else{
				return false;
			}
			break;
	}
}




std::size_t dump_reader::block_count()
{
	std::size_t n_blocks = 0;
	block_data b;
	while( next_block(b) ){
		n_blocks++;
		if( n_blocks % 100 == 0 ){
			std::cerr << "At block " << n_blocks << "...\n";
		}
	}
	
	at_eof = false;
	rewind();
	
	return n_blocks;
}



void dump_reader::rewind()
{
	if( infile ){
		infile.close();
	}
#ifdef HAVE_BOOST_GZIP
	if( infile_filt ){
		infile_filt.clear();
	}
#endif // HAVE_BOOST_GZIP
	
	init_infile( infile_name );
}


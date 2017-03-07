#include "reader_core_gzip.h"

#include <fstream>

#ifdef HAVE_BOOST_GZIP
#  include <boost/iostreams/filtering_stream.hpp>
#  include <boost/iostreams/filter/zlib.hpp>
#  include <boost/iostreams/filter/gzip.hpp>
#endif // HAVE_BOOST_GZIP

reader_core_gzip::reader_core_gzip( const std::string &fname )
	: infile( fname, std::ios_base::in | std::ios_base::binary )
{
	in.push( boost::iostreams::gzip_decompressor() );
	in.push( infile );
}

reader_core_gzip::~reader_core_gzip()
{}

bool reader_core_gzip::getline( std::string &line )
{
	if( std::getline( in, line ) ){
		return true;
	}else{
		return false;
	}
}

void reader_core_gzip::rewind()
{
	
}

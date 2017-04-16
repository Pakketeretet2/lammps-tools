#include "text_readers.h"

#include <iostream>
#include <fstream>



text_reader_plain::text_reader_plain( const std::string &fname )
	: text_reader( fname ), in(nullptr)
{
	in = new std::ifstream(fname);
}

text_reader_plain::~text_reader_plain()
{
	if( in ) delete in;
}

bool text_reader_plain::getline( std::string &line )
{
	/*
	if( in ){
		std::cerr << "in exists.\n";
	}else{
		std::cerr << "in does not exist.\n";
	}
	*/
	if( std::getline( *in, line ) ){
		return true;
	}else{
		return false;
	}
}

int text_reader_plain::peek()
{
	return in->peek();
}

bool text_reader_plain::eof() const
{
	if( in ) return in->eof();
	else     return false;
}

bool text_reader_plain::good() const
{
	if( in ) return in->good();
	else     return false;
}


text_reader_gzip::text_reader_gzip( const std::string &fname )
	: text_reader(fname), in_file(nullptr)
{
#ifndef HAVE_BOOST_GZIP
	std::cerr << "Cannot read gzip files without boost support!\n";
	std::terminate();
#endif
	in_file = new std::ifstream( fname, std::ios_base::in );
	if( !in_file ){
		std::cerr << "Failed to open file " << fname << "!\n";
		std::terminate();
	}
#ifdef HAVE_BOOST_GZIP
	in.push( boost::iostreams::gzip_decompressor() );
	in.push( *in_file );
#endif // HAVE_BOOST_GZIP
}

text_reader_gzip::~text_reader_gzip()
{
	if( in_file ) delete in_file;
}


bool text_reader_gzip::getline( std::string &line )
{
#ifdef HAVE_BOOST_GZIP
	return static_cast<bool>( std::getline(in,line) );
#else
	return false;
#endif // HAVE_BOOST_GZIP
}

int text_reader_gzip::peek()
{
#ifdef HAVE_BOOST_GZIP
	return in.peek();
#else
	return 0;
#endif
}

bool text_reader_gzip::eof()  const
{
#ifdef HAVE_BOOST_GZIP
	return in.eof();
#else
	return false;
#endif
}

bool text_reader_gzip::good() const
{
#ifdef HAVE_BOOST_GZIP
	return in.good();
#else
	return false;
#endif
}

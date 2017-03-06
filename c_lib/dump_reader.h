#ifndef DUMP_READER_H
#define DUMP_READER_H

#include "util.h"
#include "block_data.h"

#include <string>
#include <vector>
#include <fstream>
#include <array>

#define HAVE_BOOST_GZIP
#ifdef HAVE_BOOST_GZIP
#  include <boost/iostreams/filter/zlib.hpp>
#  include <boost/iostreams/filtering_stream.hpp>
#endif


class dump_reader
{
public:
	dump_reader( const std::string &fname );
	
	
	virtual ~dump_reader();

	virtual bool next_block( block_data &block );
	virtual bool last_block( block_data &block );

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
		PLAIN   = 0,
		GZIP    = 1,
		ISTREAM = 2
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

#endif // DUMP_READER_H

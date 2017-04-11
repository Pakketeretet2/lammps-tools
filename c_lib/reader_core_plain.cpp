#include "reader_core_plain.h"

#include <fstream>

reader_core_plain::reader_core_plain( const std::string &fname ) : got_file(true),
                                                                   lc(0), debug(false)
{
	in = new std::ifstream( fname );
	if( !in ){
		std::cerr << "Failed to open input file!\n";
		std::terminate();
	}
}

reader_core_plain::reader_core_plain( std::istream &istream ) : got_file(false),
                                                                lc(0), debug(false)
{
	in = &istream;
}

reader_core_plain::~reader_core_plain()
{
	if( got_file && in ) delete in;
}

bool reader_core_plain::getline( std::string &line )
{
	if( debug ) std::cerr << "Attempting getline...\n";
	if( std::getline( *in, line ) ){
		lc ++;
		if( debug ){
			std::cerr << "  Got line " << lc << ": \"" << line << "\".\n";
		}
		return true;
	}else{
		if( debug ){
			std::cerr << "  Failed to get line " << lc << "!\n";
		}

		return false;
	}
}

void reader_core_plain::rewind()
{
	// Rewind only makes sense on a file.
	if( got_file ){
		in->clear();
		in->seekg(0);
		lc = 0;
	}else{
		std::cerr << "Cannot rewind stream! Ignoring call...\n";
	}
}

#include "reader_core_plain.h"

#include <fstream>

reader_core_plain::reader_core_plain( const std::string &fname ) : got_file(true),
                                                                   lc(0)
{
	in = new std::ifstream( fname );
}

reader_core_plain::reader_core_plain( std::istream &istream ) : got_file(false),
                                                                lc(0)
{
	in = &istream;
}

reader_core_plain::~reader_core_plain()
{
	if( got_file && in ) delete in;
}

bool reader_core_plain::getline( std::string &line )
{
	if( std::getline( *in, line ) ){
		lc ++;
		return true;
	}else{
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

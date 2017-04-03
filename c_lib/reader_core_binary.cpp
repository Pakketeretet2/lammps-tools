#include "reader_core_binary.h"
#include <fstream>

reader_core_bin::reader_core_bin( const std::string &fname ) :
	in(nullptr), got_external_stream(false), current_byte(0)
{
	in = new std::ifstream( fname, std::ios_base::binary |
	                        std::ios_base::in );
	if( !in ){
		std::cerr << "Failed to open file " << fname << "!\n";
		std::terminate();
	}
}

reader_core_bin::reader_core_bin( std::istream &istream ) :
	in( &istream ), got_external_stream(true), current_byte(0)
{
	if( !in ){
		std::cerr << "Faulty input stream!\n";
		std::terminate();
	}
}

reader_core_bin::~reader_core_bin()
{
	if( !got_external_stream ) delete in;
}

void reader_core_bin::get_char_data( char *dest, uint n )
{
	in->read(dest,n);
	current_byte += n;
}

unsigned long int reader_core_bin::get_current_byte()
{
	return current_byte;
}

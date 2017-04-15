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
	return static_cast<bool>(std::getline( *in, line ));
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
	

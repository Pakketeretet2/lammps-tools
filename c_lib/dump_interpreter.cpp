#include "dump_interpreter.h"

dump_interpreter::dump_interpreter( const std::string &dname )
	: dname_(dname)
{}

dump_interpreter::dump_interpreter( std::istream &in )
	: dname_("<input stream>")
{}
	

int dump_interpreter::next_block( block_data &block )
{
	std::cerr << "Do not use dump_interpreter directly! "
	          << "Use a derived class!\n";
	std::terminate();
	return false;
}

int dump_interpreter::next_block_meta( block_data &block )
{
	std::cerr << "Do not use dump_interpreter directly! "
	          << "Use a derived class!\n";
	std::terminate();
	return false;
}

int dump_interpreter::next_block_body( block_data &block )
{
	std::cerr << "Do not use dump_interpreter directly! "
	          << "Use a derived class!\n";
	std::terminate();
	return false;
}


int dump_interpreter::last_block( block_data &block )
{
	block_data last_block;
	int status = next_block( block );
	if( !status ){
		do {
			last_block = block;
		} while ( next_block( block ) );
		return 0;
	}else{
		return -1;
	}
}


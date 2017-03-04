#include "dump_reader_lammps.h"

dump_reader_lammps::dump_reader_lammps( const std::string &fname )
{
	
}

dump_reader_lammps::dump_reader_lammps( std::istream &in )
{
}
		
dump_reader_lammps::~dump_reader_lammps()
{
}

bool dump_reader_lammps::next_block( block_data &block )
{
	bool success = false;
	return success;	
}

bool dump_reader_lammps::last_block( block_data &block )
{
	bool success = false;
	return success;
}



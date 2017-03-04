#ifndef DUMP_READER_LAMMPS_H
#define DUMP_READER_LAMMPS_H

#include <string>

#include "block_data.h"

class dump_reader_lammps
{
public:
	dump_reader_lammps( const std::string &fname );
	dump_reader_lammps( std::istream &in );
		
	virtual ~dump_reader_lammps();

	virtual bool next_block( block_data &block );
	virtual bool last_block( block_data &block );


private:
	
};



#endif // DUMP_READER_LAMMPS_H

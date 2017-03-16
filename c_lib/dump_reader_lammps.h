#ifndef DUMP_READER_LAMMPS_H
#define DUMP_READER_LAMMPS_H

#include <string>
#include "block_data.h"

/*!
  @file dump_reader_lamms.h
  @brief Reads in LAMMPS dump files. Supports plain text
         and gzip (if compiled with boost support)

  \ingroup cpp_lib
*/



class dump_reader_lammps
{
public:
	/// Read in dump from file.
	dump_reader_lammps( const std::string &fname );
	/// Read in dump from /dev/stdin
	dump_reader_lammps( std::istream &in );
		
	virtual ~dump_reader_lammps();

	virtual bool next_block( block_data &block );
	virtual bool last_block( block_data &block );


private:
	
};



#endif // DUMP_READER_LAMMPS_H

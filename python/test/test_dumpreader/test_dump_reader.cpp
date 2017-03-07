#include "dump_reader.h"
#include "util.h"
#include "block_data_writers.h"

#include <iostream>

int main( int argc, char **argv )
{
	if( argc < 2 ){
		std::cerr << "Pass a dump file!\n";
		return -1;
	}
	
	std::string fname = argv[1];
	if( ends_with( fname, ".dump" ) ){
		std::cerr << "Dump file ends with .dump!\n";
	}else{
		std::cerr << "Dump file does not end with .dump!\n";		
	}

	dump_reader d( fname );
	block_data b;

	while( d.next_block(b) ){
		// write_block_lammps_dump( b, std::cout );
	}



	return 0;
}

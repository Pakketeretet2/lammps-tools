#include "dump_reader.h"
#include "block_data_writers.h"

#include <fstream>

int main( int argc, char **argv )
{
	dump_reader d( "triangles.gsd" );
	block_data b;
	
	bool success = d.next_block( b );
	while( success ){
		std::cerr << "Read block at t = " << b.tstep << ".\n";
		std::cerr << "x[2] = " << b.x[2][0] << ", " << b.x[2][1] << ", " << b.x[2][2] << "\n";
		success = d.next_block( b );
	}

	// Attempt to write the block to a lammps and a gsd format:
	
	dump_reader d2( "triangles.gsd" );
	d2.next_block( b );
	d2.next_block( b );
	d2.next_block( b );
	d2.next_block( b );
	
	write_block_hoomd_gsd( b, "test1.gsd" );

	std::ofstream lmp_file( "test1.dump" );
	write_block_lammps_dump( b, lmp_file );
	

	return 0;
}

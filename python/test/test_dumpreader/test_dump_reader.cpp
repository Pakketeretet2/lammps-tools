#include "dump_reader_lammps.h"
#include "dump_reader.h"


int main( int argc, char **argv )
{
	std::string fname = argv[1];
	dump_reader d( fname );
	
	return 0;
}

#include "dump_interpreter_dcd.h"




int dump_interpreter_dcd::next_block( reader_core *r, block_data &block )
{
	bool success = false;

	// Reinterpret r as a binary reader.
	reader_core_bin *rb = dynamic_cast<reader_core_bin*>(r);
	if( !rb ){
		std::cerr << "Dynamic cast of reader_core to reader_core_bin!\n";
	}
	if( rb->peek() == EOF ){
		if( rb ){
			// Assume at EOF?
			return 1;
		}
	}


	// First four bytes are garbage.
	char c[4];
	rb->get_char_data( c, 4 );
	char cord[5];
	rb->get_char_data( cord, 4 );
	cord[4] = '\0';
	if( std::string( cord ) != "CORD" ){
		std::cerr << "NAMD DCD format inconsistent! " << cord
		          << " != " << "CORD!\n";
		std::terminate();
	}

	
	std::cerr << "Sizes in byes of various data types:\n"
	          << "   char    " << sizeof(char) << "\n"
	          << "   short   " << sizeof(short) << "\n"
	          << "   int     " << sizeof(int) << "\n"
	          << "   long    " << sizeof(long) << "\n"
	          << "   float   " << sizeof(float) << "\n"
	          << "   double  " << sizeof(double) << "\n"
	          << "   ldouble " << sizeof(long double) << "\n";


	int nset = 0;
	int istrt, nsavc, natom_nfreat;
	
	nset = rb->get_data<int>();
	std::cerr << "NSET = " << nset << "\n";
	istrt = rb->get_data<int>();
	std::cerr << "ISTRT = " << istrt << "\n";
	nsavc = rb->get_data<int>();
	std::cerr << "NSAVC = " << nsavc << "\n";
	natom_nfreat = rb->get_data<int>();
	std::cerr << "NATOM-NFREAT = " << natom_nfreat << "\n";
	
	int five_zeros[5];
	rb->get_data( five_zeros, 5 );
	five_zeros[5] = '\0';
	std::cerr << "5 zeros: ";
	for( int i = 0; i < 5; ++i ){
		std::cerr << five_zeros[i] << " ";
	}
	std::cerr << "\n";

	float dt;
	dt = rb->get_data<float>(  ) * TIMEFACTOR;
	std::cerr << "DELTA = " << dt << " fs\n";
	int nine_zeros[9];
	rb->get_data( nine_zeros, 9 );
	std::cerr << "9 zeros: ";
	for( int i = 0; i < 9; ++i ){
		std::cerr << nine_zeros[i] << " ";
	}
	std::cerr << "\n";
	int ntitle = rb->get_data<int>();
	std::cerr << "NTITLE = " << ntitle << "\n";

	char title[33];
	rb->get_char_data( title, 32 );
	title[32] = '\0';
	std::cerr << "Title: " << title << "\n";

	int natoms = rb->get_data<int>();
	std::cerr << "NATOM = " << natoms << "\n";

	int fast_forward = 144 - rb->get_current_byte();
	std::cerr << "Reading " << fast_forward << " bytes of info.\n";
	char *remarks = new char[fast_forward+1];
	rb->get_char_data( remarks, fast_forward );
	remarks[fast_forward] = '\0';
	std::cerr << "REMARKS = " << remarks << "\n";

	std::cerr << "Read " << rb->get_current_byte() << " bytes.\n";
	for( int i = 0; i < 4; ++i ){
		float x, y, z;
		x = rb->get_data<float>();
		y = rb->get_data<float>();
		z = rb->get_data<float>();
		std::cerr << "(x,y,z) = ( " << x << ", " << y << ", "
		          << z << " ).\n";
	}
	return 0;
}
	

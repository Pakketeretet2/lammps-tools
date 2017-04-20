#include <iostream>
#include <fstream>
#include <string>

#include "block_data.h"

constexpr const int Nkeys = 5;
const char *tests[Nkeys] = { "new_block_data",
                             "set_block_data",
                             "get_block_data",
                             "free_block_data",
                             "block_data" };

void test_block_data()
{
	block_data b;

	block_data bb( 123 );
}

void test_new_block_data()
{
	block_data *b = new_block_data();
	
}
void test_get_block_data()
{
	
}
void test_set_block_data()
{
	
}
void test_free_block_data()
{
	block_data *b = new_block_data();
	free_block_data(b);
}



void print_usage()
{
	std::cerr << "Pass one or more of these keywords from the "
		"command line:\n";
	for( int i = 0; i < Nkeys; ++i ){
		std::cerr << tests[i] << " ";
	}
	std::cerr << "\n";
}


// Tests the interfaces for leaks.
int main( int argc, char **argv )
{
	if( argc <= 1 ){
		std::cerr << "Pass the interfaces you want to test for leaks!\n";
		print_usage();
	}

	for( int i = 1; i < argc; ++i ){
		std::string key = argv[i];
		if( key == tests[0] ){
			test_new_block_data();
		}else if( key == tests[1] ){
			test_set_block_data();
		}else if( key == tests[2] ){
			test_get_block_data();
		}else if( key == tests[3] ){
			test_free_block_data();
		}else if( key == tests[4] ){
			test_block_data();
		}else{
			std::cerr << "Unknown keyword.\n";
		}
	}
	
}

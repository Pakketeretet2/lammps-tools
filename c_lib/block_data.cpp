#include "block_data.h"

#include <iostream>
#include <fstream>

extern "C" {
	
block_data *new_block_data()
{
	block_data *b = new block_data;

	b->x  = nullptr;
	b->x_ = nullptr;
	b->ids = b->types = b->mol = nullptr;
	b->N = b->tstep = 0;
	b->boxline = "";
	b->atom_style = block_data::ATOMIC;
	
	return b;
}

void delete_members( block_data &b )
{
	delete [] b.x_;
	delete [] b.x;
	delete [] b.ids;
	delete [] b.types;
	delete [] b.mol;
}

void delete_block_data(block_data* b)
{
	delete_members( *b );
	delete b;
}


}



void resize( block_data &b, int N )
{
	if( N < b.N ){
		b.N = N;
	}else{
		delete_members(b);
		init( b, N );
	}
}
	
void init( block_data &b, int N )
{
	b.x  = new py_float*[N];
	b.x_ = new py_float[3*N];

	b.ids   = new py_int[N];
	b.types = new py_int[N];
	b.mol   = new py_int[N];
	b.N = N;

	
}




void print_block_data_lmp( const block_data &b, const std::string &s )
{
	print_block_data_lmp( b, s, std::ofstream::out );
}
                           
void print_block_data_lmp( const block_data &b, const std::string &s,
                           std::ios_base::openmode mode )
{
	std::ofstream o( s, mode );
	print_block_data_lmp( b, o );
}

void print_block_data_lmp( const block_data &b, std::ostream &o )
{
	
}

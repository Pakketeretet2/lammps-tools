#include "block_data.h"
#include <iostream>
#include <fstream>
#include <algorithm>

extern "C" {


block_data::block_data() : x(nullptr), x_(nullptr), ids(nullptr),
                           types(nullptr), mol(nullptr), N(0),
                           tstep(0), xlo{0,0,0}, xhi{0,0,0},
                           periodic(0), atom_style(atom_styles::ATOMIC),
                           boxline("       ")
{}

block_data::block_data( int N ) : x(nullptr), x_(nullptr), ids(nullptr),
                                  types(nullptr), mol(nullptr), N(0),
                                  tstep(0), xlo{0,0,0}, xhi{0,0,0},
                                  periodic(0), atom_style(atom_styles::ATOMIC),
                                  boxline("       ")
{
	init(N);
}


void block_data::delete_members()
{
	delete [] x_;
	delete [] x;
	delete [] ids;
	delete [] types;
	delete [] mol;
}




block_data::~block_data()
{
	delete_members();
}

void block_data::resize( int NN )
{
	if( NN < N ){
		N = NN;
	}else{
		delete_members();
		init( NN );
	}
}
	
void block_data::init( int NN )
{
	N = NN;
	x  = new py_float*[N];
	x_ = new py_float[3*N];

	ids   = new py_int[N];
	types = new py_int[N];
	mol   = new py_int[N];

	for( py_int i = 0; i < N; ++i ){
		x[i] = x_ + 3*i;
	}
}

} // extern "C";


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
	o << "ITEM: TIMESTEP\n" << b.tstep << "\n";
	o << "ITEM: NUMBER OF ATOMS\n" << b.N << "\n";
	o << b.boxline << "\n";
	o << b.xlo[0] << " " << b.xhi[0] << "\n";
	o << b.xlo[1] << " " << b.xhi[1] << "\n";
	o << b.xlo[2] << " " << b.xhi[2] << "\n";

	if( b.atom_style == atom_styles::MOLECULAR ){
		o << "ITEM: ATOMS id mol type x y z\n";
		for( py_int i = 0; i < b.N; ++i ){
			o << b.ids[i] << " " << b.mol[i] << " " << b.types[i]
			  << " " << b.x[i][0] << " " << b.x[i][1] << " "
			  << b.x[i][2] << "\n";
		}
	}else{
		o << "ITEM: ATOMS id type x y z\n";
		for( py_int i = 0; i < b.N; ++i ){
			o << b.ids[i] << " " << b.types[i]
			  << " " << b.x[i][0] << " " << b.x[i][1] << " "
			  << b.x[i][2] << "\n";
		}
	}
	
}


void block_data_from_foreign( const void *x, py_int N, const py_int *ids,
                              const py_int *types, const py_int *mol,
                              py_int periodic,
                              const py_float *xlo, const py_float *xhi,
                              py_int dims, py_int tstep, const char *box_line,
                              block_data &b )
{
	b.resize( N );
	const py_float *xx = static_cast< const py_float* const>(x);
	for( int i = 0; i < N; ++i ){
		b.x[i][0] = xx[3*i];
		b.x[i][1] = xx[3*i+1];
		b.x[i][2] = xx[3*i+2];
		
	}

	std::copy( ids, ids + N, b.ids );
	std::copy( types, types + N, b.types );
	if( mol ){
		std::copy( mol, mol + N, b.mol );
		b.atom_style = atom_styles::MOLECULAR;
	}else{
		b.atom_style = atom_styles::ATOMIC;
	}

	b.periodic = periodic;
	b.tstep = tstep;
	
	b.xlo[0] = xlo[0];
	b.xlo[1] = xlo[1];
	b.xlo[2] = xlo[2];
	b.xhi[0] = xhi[0];
	b.xhi[1] = xhi[1];
	b.xhi[2] = xhi[2];
	
	b.boxline = box_line;
}

void copy( block_data &b, const block_data &source )
{
	b.resize( source.N );
	
	b.xlo[0] = source.xlo[0];
	b.xlo[1] = source.xlo[1];
	b.xlo[2] = source.xlo[2];
	
	b.xhi[0] = source.xhi[0];
	b.xhi[1] = source.xhi[1];
	b.xhi[2] = source.xhi[2];
	
	b.tstep = source.tstep;
	b.atom_style = source.atom_style;
	b.boxline = source.boxline;
	b.periodic = source.periodic;

	/*
	std::cerr << "Copying " << source.other_cols.size()
	          << " other cols as well...\n";
	*/
	
	b.other_cols.resize( source.other_cols.size() );
	for( int i = 0; i < b.other_cols.size(); ++i ){
		b.other_cols[i].data.resize( source.N );
		b.other_cols[i].header = source.other_cols[i].header;
	}
	
	for( int i = 0; i < source.N; ++i ){
		b.x[i][0]  = source.x[i][0];
		b.x[i][1]  = source.x[i][1];
		b.x[i][2]  = source.x[i][2];

		b.ids[i]   = source.ids[i];
		b.types[i] = source.types[i];

		if( b.atom_style == atom_styles::MOLECULAR ){
			b.mol[i] = source.mol[i];
		}

		for( std::size_t j = 0; j < source.other_cols.size(); ++j ){
			b.other_cols[j].data[i] = source.other_cols[j].data[i];
		}
	}	
}


extern "C" {

block_data *new_block_data()
{
	block_data *b = new block_data();
	std::cerr << "Made block data at " << b << ".\n";
	return b;
}
	
void set_block_data( block_data *b,
                     py_int N, py_int t, py_float *x, py_int *ids,
                     py_int *types, py_int *mol,
                     py_float *xlo, py_float *xhi, py_int periodic,
                     const char* boxline )
{
	std::cerr << "Setting block data on block at " << b << ".\n";
	b->N = N;
	b->tstep = t;

	b->xlo[0] = xlo[0];
	b->xlo[1] = xlo[1];
	b->xlo[2] = xlo[2];

	b->xhi[0] = xhi[0];
	b->xhi[1] = xhi[1];
	b->xhi[2] = xhi[2];

	b->boxline = boxline;
	b->resize(N);
	
	std::cerr << "Boxline is " << boxline << "\n";
	std::cerr << "Resized block data at " << b << " to N = " << N << ".\n";

	for( int i = 0; i < N; ++i ){
		
		b->x[i][0]  = x[ 3*i + 0 ];
		b->x[i][1]  = x[ 3*i + 1 ];
		b->x[i][2]  = x[ 3*i + 2 ];
		
		b->ids[i]   = ids[i];
		b->types[i] = types[i];

		if( mol ){
			b->mol[i] = mol[i];
		}
	}

	std::cerr << "After setting, this is block_data:\n";
	print_block_data_lmp( *b, std::cerr );
	
	
}
	
void free_block_data( block_data *b )
{
	std::cerr << "Freeing block data at " << b << ".\n";
	delete b;
}

}

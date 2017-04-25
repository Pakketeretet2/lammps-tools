#include "block_data.h"
#include <iostream>
#include <fstream>
#include <algorithm>

extern "C" {


block_data::block_data() : x(nullptr), x_(nullptr), ids(nullptr),
                           types(nullptr), mol(nullptr), N(0),
                           tstep(0), xlo{0,0,0}, xhi{0,0,0},
                           periodic(0), atom_style(atom_styles::ATOMIC),
                           boxline("       "), Ntypes(0), mass(nullptr)
{}

block_data::block_data( int N ) : x(nullptr), x_(nullptr), ids(nullptr),
                                  types(nullptr), mol(nullptr), N(0),
                                  tstep(0), xlo{0,0,0}, xhi{0,0,0},
                                  periodic(0), atom_style(atom_styles::ATOMIC),
                                  boxline("       "), Ntypes(0), mass(nullptr)
{
	init(N);
}

	

block_data::block_data( const block_data &o )
{
	*this = o;
}


// TODO: Rewrite this into copy & swap:
void block_data::operator=( const block_data &o )
{
	resize( o.N );
	// This needs to copy _all_ the fields of block_data, _except_ x
	// because x refers to the memory of x_!
	std::copy( o.x_,    o.x_    + 3*N, x_ );
	std::copy( o.ids,   o.ids   + N,   ids );
	std::copy( o.types, o.types + N,   types );
	if( o.atom_style == atom_styles::MOLECULAR ){
		if( !o.mol ){
			std::cerr << "Copying block dat with atom_style "
			          << "molecular but it's missing mol array!\n";
			std::terminate();
		}
		std::copy( o.mol, o.mol + N, mol );
	}

	std::copy( o.xlo, o.xlo + 3, xlo );
	std::copy( o.xhi, o.xhi + 3, xhi );

	tstep      = o.tstep;
	periodic   = o.periodic;
	atom_style = o.atom_style;
	boxline    = o.boxline;

	// Copy the other cols:
	int N_other_cols = o.other_cols.size();
	other_cols.resize( N_other_cols );
	
	for( int i = 0; i < o.other_cols.size(); ++i ){
		
		int M = o.other_cols[i].data.size();
		other_cols[i].resize( M );
		other_cols[i].header = o.other_cols[i].header;
		for( int j = 0; j < M; ++i ){
			other_cols[i].data[j] = o.other_cols[i].data[j];
		}
	}
}

/*
void block_data::operator=( const block_data &o )
{
	block_data tmp(o);
	tmp.swap(*this);
}

void block_data::swap( block_data &o ) throw()
{
	std::swap( *this, o );
}
*/


void block_data::copy_meta( const block_data &o )
{
	this->atom_style = o.atom_style;
	this->boxline    = o.boxline;
	this->tstep      = o.tstep;
	
	this->xlo[0] = o.xlo[0];
	this->xlo[1] = o.xlo[1];
	this->xlo[2] = o.xlo[2];
		
	this->xhi[0] = o.xhi[0];
	this->xhi[1] = o.xhi[1];
	this->xhi[2] = o.xhi[2];
	this->N      = o.N;

}

void block_data::delete_members()
{
	delete [] x_;
	delete [] x;
	delete [] ids;
	delete [] types;
	delete [] mol;

	if( mass ) delete [] mass;
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

void block_data::init_per_type_arrays( int NNtypes )
{
	Ntypes = NNtypes;
	mass = new double[NNtypes+1];
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
	
	/*
	  std::cerr << "Boxline is " << boxline << "\n";
	  std::cerr << "Resized block data at " << b << " to N = " << N << ".\n";
	*/

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

	/*
	  std::cerr << "After setting, this is block_data:\n";
	  print_block_data_lmp( *b, std::cerr );
	*/
	
	
}


void get_block_data( block_data *b,
                     py_int *N, py_int *t, py_float *x, py_int *ids,
                     py_int *types, py_int *mol,
                     py_float *xlo, py_float *xhi, py_int *periodic,
                     char *boxline )
{
	std::cerr << "Setting block data on block at " << b << ".\n";
	*N = b->N;
	*t = b->tstep;

	xlo[0] = b->xlo[0];
	xlo[1] = b->xlo[1];
	xlo[2] = b->xlo[2];

	xhi[0] = b->xhi[0];
	xhi[1] = b->xhi[1];
	xhi[2] = b->xhi[2];

	b->boxline.copy( boxline, 0, b->boxline.size() );
	boxline[b->boxline.size()] = '\0';
	// b->resize(N);
	
	for( int i = 0; i < b->N; ++i ){
		
		x[ 3*i + 0 ] = b->x[i][0];
		x[ 3*i + 1 ] = b->x[i][1];
		x[ 3*i + 2 ] = b->x[i][2];
		             
		ids[i]   = b->ids[i];
		types[i] = b->types[i];

		if( b->mol && mol ){
			mol[i] = b->mol[i];
		}
	}
	
	/*
	  std::cerr << "After setting, this is block_data:\n";
	  print_block_data_lmp( *b, std::cerr );
	*/
	
	
}

	
void free_block_data( block_data *b )
{
	std::cerr << "Freeing block data at " << b << ".\n";
	delete b;
}

}

#include "dump_reader.h"
#include "id_map.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>

#include <boost/iostreams/filter/gzip.hpp>

#include "reader_core_plain.h"
#include "reader_core_gzip.h"
#include "reader_core_binary.h"

#include "dump_interpreter_lammps.h"
#include "dump_interpreter_dcd.h"
#include "dump_interpreter_gsd.h"


struct dummy_stream
{
	template <typename T>
	dummy_stream &operator<<( const T &t )
	{
		return *this;
	}

	void flush(){}
};

struct debug_stream
{
	debug_stream( std::ostream &out ) : o(out){}
		
	template <typename T>
	debug_stream &operator<<( const T &t )
	{
		o << t;
		o.flush();
		return *this;
	}

	void flush()
	{
		o.flush();
	}
	
	std::ostream &o;
};



bool reader_core::getline( std::string &line )
{
	std::cerr << "Do not use reader_core::getline. Use a derived class!\n";
	std::terminate();
	return false;
}

void reader_core::rewind()
{
	std::cerr << "Do not use reader_core::rewind. Use a derived class!\n";
	std::terminate();
}


bool dump_interpreter::next_block( reader_core *r, block_data &block )
{
	std::cerr << "Do not use dump_interpreter directly! "
	          << "Use a derived class!\n";
	std::terminate();
	return false;
}

bool dump_interpreter::next_block_meta( reader_core *r, block_data &block )
{
	std::cerr << "Do not use dump_interpreter directly! "
	          << "Use a derived class!\n";
	std::terminate();
	return false;
}

bool dump_interpreter::next_block_body( reader_core *r, block_data &block )
{
	std::cerr << "Do not use dump_interpreter directly! "
	          << "Use a derived class!\n";
	std::terminate();
	return false;
}


bool dump_interpreter::last_block( reader_core *r, block_data &block )
{
	block_data last_block;
	bool success = next_block( r, block );
	if( success ){
		do{
			last_block = block;
		}while( next_block( r, block ) );
		return success;
	}else{
		return false;
	}
}




// ***************** DUMPREADER STUFF ********************
dump_reader::dump_reader( std::istream &stream, int dformat ) :
	dump_format( dformat ), at_eof(false)
{
	file_format = ISTREAM;
	dump_format = dformat;
	setup_interpreter();
	setup_reader( stream );
}



dump_reader::dump_reader( const std::string &fname, int fformat, int dformat )
	: dump_format( dformat ), file_format( fformat )
{
	if( dump_format < 0 ) guess_dump_type( fname );
	if( file_format < 0 ) guess_file_type( fname );

	post_constructor( fname );
}

void dump_reader::guess_file_type( const std::string &fname )
{
	if( ends_with( fname, ".gz" ) ){
		file_format = GZIP;
	}else{
		std::cerr << "File format assumed to be plain text.\n";
		file_format = PLAIN;
	}
}

void dump_reader::guess_dump_type( const std::string &fname )
{
	if( ends_with( fname, ".gsd" ) ){
		dump_format = GSD;
	}else if( ends_with( fname, ".dump" ) ){
		dump_format = LAMMPS;
	}else if( ends_with( fname, ".dump.gz" ) ){
		dump_format = LAMMPS;
	}else if( ends_with( fname, ".dcd" ) ){
		dump_format = DCD;
		file_format = BIN;
	}else{
		std::cerr << "Extension of dump file " << fname
		          << " not recognized! If you are certain of the "
		          << "format, explicitly state the dump format!\n";
		std::terminate();
	}
}

void dump_reader::post_constructor( const std::string &fname )
{
	if( dump_format != GSD ){
		setup_interpreter();
		setup_reader( fname );
	}else if( dump_format == GSD ){
		setup_interpreter_gsd( fname );
		reader = nullptr; // gsd uses it's own reader.
	}else{
		std::cerr << "Dunno what to do with dump format "
		          << dump_format << "!\n";
		std::terminate();
	}
}




dump_reader::~dump_reader()
{
	if( reader ) delete reader;
	if( interp ) delete interp;
}


void dump_reader::setup_reader( std::istream &stream )
{
	reader = new reader_core_plain( stream );
	if( !reader ){
		std::cerr << "Failure setting up reader!\n";
		std::terminate();
	}
}

void dump_reader::setup_reader( const std::string &fname )
{
	std::cerr << "Setting up reader for format " << file_format << ".\n";
	if( !file_exists( fname ) ){
		std::cerr << "File " << fname << " doesn't exist!\n";
		std::terminate();
	}
	
	if( file_format == GZIP ){
		reader = new reader_core_gzip( fname );
	}else if( file_format == PLAIN ){
		reader = new reader_core_plain( fname );
	}else if( file_format == BIN ){
		reader = new reader_core_bin( fname );
	}else{
		std::cerr << "File format for reader not recognized!\n";
		std::terminate();
	}

	if( !reader || !(*reader) ){
		std::cerr << "Failure setting up reader! Perhaps file doesn't exist?\n";
		std::terminate();
	}
}

void dump_reader::setup_interpreter()
{
	if( dump_format == LAMMPS ){
		interp = new dump_interpreter_lammps();
	}else if( dump_format == DCD ){
		std::cerr << "Creating dcd interpreter.\n";
		interp = new dump_interpreter_dcd();
	}else{
		std::cerr << "Unrecognized dump format!\n";
		std::terminate();
	}
}

void dump_reader::setup_interpreter_gsd( const std::string &fname )
{
	interp = new dump_interpreter_gsd( fname );
}



bool dump_reader::next_block( block_data &block )
{
	return interp->next_block( reader, block );
}

bool dump_reader::last_block( block_data &block )
{
	bool success = interp->next_block( reader, block );
	if( !success ) return false;
	
	while( success ){
		success = interp->next_block( reader, block );
	}
	return true;
}



bool dump_reader::skip_blocks( int Nblocks )
{
	bool success = false;
	for( int i = 0; i < Nblocks; ++i ){
		block_data b;
		success = interp->next_block_meta( reader, b );
		if( !success ){
			return success;
		}
		// line is at "ITEM: TIMESTEP now. Skip b.N + 1 line."
		std::string tmp;
		for( int i = 0; i < b.N; ++i ){
			reader->getline( tmp );
		}
	}
	return success;
}


bool dump_reader::skip_block ( )
{
	block_data b;
	return next_block( b );
}



std::size_t dump_reader::block_count()
{
	std::size_t i = 0;
	bool success = false;
	while( true ){
		block_data b;
		success = next_block( b );
		if( success ){
			++i;
		}else{
			return i;
		}
	}
}


// **************** For interfacing with foreign languages *********************
extern "C" {

dump_reader_handle *get_dump_reader_handle( const char *dname,
                                            py_int dformat, py_int fformat )
{
	std::cerr << "dname = " << dname << "\n";
	dump_reader_handle *dh = new dump_reader_handle;
	dh->last_block = nullptr;
	dh->reader     = nullptr;
	std::string name( dname );
	dh->reader = new dump_reader( name, fformat, dformat );
	dh->last_block = new block_data;
	
	return dh;
}

void release_dump_reader_handle( dump_reader_handle *dh )
{
	if( dh ){
		if( dh->reader ) delete dh->reader;
		if( dh->last_block ) delete dh->last_block;
		delete dh;
	}
}

bool dump_reader_next_block( dump_reader_handle *dh )
{
	block_data block;
	std::cerr << "Calling next_block on dh @ " << dh << "\n";
	bool status = dh->reader->next_block( block );
	if( !status ){
		std::cerr << "Failed to get next block\n";
		return false;
	}
	
	// std::cerr << "Got next block of " << block.N
	//           << " particles at time " << block.tstep << ".\n";
	block_data *address = dh->last_block;
	std::cerr << "Got next block, it's @ " << address << "\n";
	*address = block;
	return true;
}

void dump_reader_get_block_meta( dump_reader_handle *dh,
                                 py_int *tstep, py_int *N,
                                 py_float *xlo, py_float *xhi,
                                 py_int *periodic, char *boxline,
                                 py_int *atom_style )
{
	block_data *lb = dh->last_block;
	*tstep = lb->tstep;
        *N     = lb->N;

	xlo[0] = lb->xlo[0];
	xlo[1] = lb->xlo[1];
	xlo[2] = lb->xlo[2];

	xhi[0] = lb->xhi[0];
	xhi[1] = lb->xhi[1];
	xhi[2] = lb->xhi[2];

	*periodic   = lb->periodic;
	*atom_style = lb->atom_style;
	
	// std::cerr << "Atom style in C++ part is " << lb->atom_style << "\n";
	
	if( boxline ){
		// std::cerr << "Before writing, external boxline is ";
		for( int i = 0; i < lb->boxline.size(); ++i ){
			// std::cerr << static_cast<char>(boxline[i]);
		}
		// std::cerr << "\n";
		for( int i = 0; i < lb->boxline.size(); ++i ){
			boxline[i] = static_cast<char>( lb->boxline[i] );
		}
		
		// boxline[lb->boxline.size()] = '\0';
		// std::cerr << "Internal boxline is ";
		for( int i = 0; i < lb->boxline.size(); ++i ){
			// std::cerr << lb->boxline[i];
		}
		// std::cerr << "\n";
		// std::cerr << "After writing, external boxline is ";
		for( int i = 0; i < lb->boxline.size(); ++i ){
			// std::cerr << static_cast<char>(boxline[i]);
		}
		// std::cerr << "\n";
	}
}
	
void dump_reader_get_block_data( dump_reader_handle *dh,
                                 py_int N, py_float *x, py_int *ids,
                                 py_int *types, py_int *mol )
{
	// std::cerr << "Writing " << N << " particles.\n";
	if( !x || !ids || !types ){
		std::cerr << "One of the essential arrays not "
		          << "allocated properly.\n";
	}

	for( py_int i = 0; i < N; ++i ){
		x[3*i + 0] = dh->last_block->x[i][0];
		x[3*i + 1] = dh->last_block->x[i][1];
		x[3*i + 2] = dh->last_block->x[i][2];

		ids[i] = dh->last_block->ids[i];
		types[i] = dh->last_block->types[i];
	}

	
	if( mol ){
		if( dh->last_block->atom_style == atom_styles::ATOMIC ){
			/*
			std::cerr << "WARNING: Mol requested but atom style "
			          << "does not support it! Ignoring...\n";
			*/
		}else if( !dh->last_block->mol ){
			std::cerr << "WARNING: Mol requested but block does "
			          << "not have mol array! Ignoring...\n";
		}else{
			for( int i = 0; i < N; ++i ){

				mol[i] = dh->last_block->mol[i];
			}
		}
	}

}

bool dump_reader_fast_forward( dump_reader_handle *dh,
                               py_int N_blocks )
{
	return dh->reader->skip_blocks( N_blocks );
}

	
} // extern "C"

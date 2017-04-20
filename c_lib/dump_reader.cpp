#include "dump_reader.h"
#include "id_map.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>

#include <boost/iostreams/filter/gzip.hpp>

#include "dump_interpreter_lammps.h"
#include "dump_interpreter_dcd.h"
#include "dump_interpreter_gsd.h"


dump_reader::dump_reader( const std::string &fname )
	: dump_format(-1), file_format(-1), interp(nullptr)
{
	guess_dump_type( fname );
	guess_file_type( fname );

	setup_interpreter( fname );
}



dump_reader::dump_reader( const std::string &fname, int dformat, int fformat )
	: dump_format( dformat ), file_format( fformat ), interp(nullptr)
{
	if( dump_format < 0 ) guess_dump_type( fname );
	if( file_format < 0 ) guess_file_type( fname );

	setup_interpreter( fname );
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
		dump_format = NAMD;
		file_format = BIN;
	}else{
		std::cerr << "Extension of dump file " << fname
		          << " not recognized! If you are certain of the "
		          << "format, explicitly state the dump format!\n";
		std::terminate();
	}
}



dump_reader::~dump_reader()
{
	if( interp ) delete interp;
}

void dump_reader::setup_interpreter( const std::string &fname )
{
	if( dump_format == LAMMPS ){
		if( file_format == PLAIN ){
			interp = new dump_interpreter_lammps( fname );
		}else{
			interp = new dump_interpreter_lammps( fname, 1 );
		}
	}else if( dump_format == GSD ){
		interp = new dump_interpreter_gsd( fname );
	}else if( dump_format == NAMD ){
		
	}
}


int dump_reader::next_block( block_data &block )
{
	return interp->next_block( block );
}

int dump_reader::last_block( block_data &block )
{
	int status = interp->next_block( block );
	while( !status ){
		status = interp->next_block( block );
	}
	return status;
}



int dump_reader::skip_blocks( int Nblocks )
{
	int status = 0;
	for( int i = 0; i < Nblocks; ++i ){
		status = skip_block();
		if( !status ){
			return status;
		}
	}
	return status;
}


int dump_reader::skip_block ( )
{
	return interp->skip_block();
}



std::size_t dump_reader::block_count()
{
	std::size_t i = 0;
	while( true ){
		block_data b;
	        int status = next_block( b );
		if( !status ){
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
        std::string name( dname );
        dh->reader = new dump_reader( name, dformat, fformat );
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

int dump_reader_next_block( dump_reader_handle *dh )
{
	// std::cerr << "Calling next_block on dh @ " << dh << "\n";
	block_data &block = *dh->last_block;
	int status = dh->reader->next_block( block );
	if( status ){
		if( status > 0 ){
			std::cerr << "Dump file at EOF.\n";
		}else{
			std::cerr << "Failed to get next block\n";
		}
		return status;
	}
        
	return status;
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
		boxline[lb->boxline.size()] = '\0';
		
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

int dump_reader_fast_forward( dump_reader_handle *dh,
                               py_int N_blocks )
{
	return dh->reader->skip_blocks( N_blocks );
}

	
} // extern "C"

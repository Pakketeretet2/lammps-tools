#include "dump_interpreter_gsd.h"

#include <iostream>

#define MY_CERR std::cerr << status << " ( " << gh << " ) "


dump_interpreter_gsd::dump_interpreter_gsd( const std::string &fname )
	: status( 0 ), gh( nullptr ), current_frame( -1 )
{
	MY_CERR << "Attempting to open gsd file " << fname << ".\n";
	gh = new gsd_handle;
	status = gsd_open( gh, fname.c_str(), GSD_OPEN_READONLY );
	MY_CERR << "Attempted to open.\n";
	if( status ){
		MY_CERR << "Error occured opening file!\n";
		std::terminate();
	}

	/*
	MY_CERR << "There are " << gsd_get_nframes( gh )
	        << " frames in this file for a total of "
	        << gh->file_size << " bytes.\n";
	MY_CERR << "File was created by " << gh->header.application
	        << " and uses " << gh->header.schema << " schema.\n";
	*/
	
}

dump_interpreter_gsd::~dump_interpreter_gsd()
{
	if( gh ){
		status = gsd_close( gh );
		MY_CERR << "Closed handle.\n";
		delete gh;
	}
	
}

int dump_interpreter_gsd::get_chunk_data( const std::string &name, void *dest )
{
	const gsd_index_entry *idx  = gsd_find_chunk( gh, current_frame,
	                                              name.c_str() );
	if( !idx ){
		/*
		MY_CERR << "find_chunk failed to find " << name
		        << " for frame " << current_frame << "\n";
		*/
		return 1;
	}

	status = gsd_read_chunk( gh, dest, idx );
	
	// MY_CERR << "Read chunk " << name << " at " << idx << ".\n";
	return status;
}


bool dump_interpreter_gsd::next_block( reader_core *, block_data &b )
{
	++current_frame;
	const gsd_index_entry *entry;

	// Read out all the chunks you want to add to your dumpfile.
	uint64_t tstep;
	uint8_t  dims;
	float    box[6];

	uint32_t N;

	bool default_type = false;
	
	
	status = get_chunk_data( "configuration/step", &tstep );
	status = get_chunk_data( "configuration/dimensions", &dims );
	status = get_chunk_data( "configuration/box", box );
	status = get_chunk_data( "particles/N", &N );

	if( status ){
		MY_CERR << "An error occured!\n";
		return false;
	}

	b.tstep = tstep;
	b.xlo[0] = -0.5*box[0];
	b.xlo[1] = -0.5*box[1];
	b.xlo[2] = -0.5*box[2];
	b.xhi[0] =  0.5*box[0];
	b.xhi[1] =  0.5*box[1];
	b.xhi[2] =  0.5*box[2];
	
	b.resize(N);

	float *x        = new float[3*N];
	uint32_t *types = new uint32_t[N];

	status = get_chunk_data( "particles/position", x );
	status = get_chunk_data( "particles/typeid",   types );
	if( status == 1 ){
		// This means types was not in the file. Probably because
		// everything is the default 1.
		std::cerr << "Didn't find particles/typeid!\n";
		default_type = true;
	}

	// MY_CERR << "Read positions and types.\n";
	for( int i = 0; i < N; ++i ){
		b.x[i][0] = x[3*i  ];
		b.x[i][1] = x[3*i+1];
		b.x[i][2] = x[3*i+2];

		b.ids[i]   = i+1;
		
	}

	if( !default_type ){
		for( int i = 0; i < N; ++i ){
			b.types[i] = types[i];
		}
	}else{
		for( int i = 0; i < N; ++i ){
			b.types[i] = 1;
		}
	}
	b.mol = nullptr;
	b.atom_style = atom_styles::ATOMIC;
	b.boxline = "ITEM: BOX BOUNDS pp pp pp";

	delete [] x;
	delete [] types;
	
	return true;
	
}

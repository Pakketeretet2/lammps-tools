#include "dump_interpreter_gsd.h"
#include "util.h"

#include <iostream>


#define MY_CERR std::cerr << status << " ( " << gh << " ) "


dump_interpreter_gsd::dump_interpreter_gsd( const std::string &fname )
	: dump_interpreter( fname ), status( 0 ),
	  gh( nullptr ), current_frame( -1 ), eof_( true ), good_( false )
{
	MY_CERR << "Attempting to open gsd file " << fname << ".\n";
	gh = new gsd_handle;
	status = gsd_open( gh, fname.c_str(), GSD_OPEN_READONLY );
	MY_CERR << "Attempted to open.\n";
	if( status ){
		MY_CERR << "Error occured opening file!\n";
		std::terminate();
	}else{
		good_ = true;
		eof_ = false;
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


int dump_interpreter_gsd::next_block( block_data &b )
{
	++current_frame;
	const gsd_index_entry *entry;

	// Read out all the chunks you want to add to your dumpfile.
	uint64_t tstep;
	uint8_t  dims;
	float    box[6];

	uint32_t N;

	bool default_type = false;
	
	if( get_chunk_data( "configuration/step", &tstep ) ) {
		// Assume EOF, because this would be first block
		eof_ = true;
		return 0;
	}
	if ( get_chunk_data( "configuration/dimensions", &dims ) ) {
		// This is not good.
		good_ = false;
		return -1;
	}
	if( get_chunk_data( "configuration/box", box ) ){
		good_ = false;
		return -1;
	}
	if( get_chunk_data( "particles/N", &N ) ){
		good_ = false;
		return -1;
	}

	b.tstep = tstep;
	b.xlo[0] = -0.5*box[0];
	b.xlo[1] = -0.5*box[1];
	b.xlo[2] = -0.5*box[2];
	b.xhi[0] =  0.5*box[0];
	b.xhi[1] =  0.5*box[1];
	b.xhi[2] =  0.5*box[2];
	
	b.resize(N);

	float *x           = new float[3*N];
	uint32_t *type_ids = new uint32_t[N];
	

	status = get_chunk_data( "particles/position", x );
	if( status == 1 ){
		// This means types was not in the file. Probably because
		// everything is the default 1.
		std::cerr << "Didn't find particles/position!\n";
		return -2;
	}
	
	status = get_chunk_data( "particles/typeid",  type_ids );
	if( status == 1 ){
		// This means types was not in the file. Probably because
		// everything is the default 1.
		std::cerr << "Didn't find particles/typeid!\n";
		default_type = true;
	}

	if( !default_type ){
		// At this point you know the type ids. Count them:
		int n_types = 0;
		for( int i = 0; i < N; ++i ){
			int type_id_inc = type_ids[i]+1;
			if( type_id_inc > n_types ) n_types = type_id_inc;
		}
		/*
		std::cerr << "Found " << n_types << " particle types. "
		          << "Reading in type names...\n";
		*/
		constexpr int buff_size = TYPE_BUFFER_SIZE;
		gsd_utf8_char *type_names = new gsd_utf8_char[ n_types * buff_size ];
		
		status = get_chunk_data( "particles/types", type_names );
		for( int i = 0; i < n_types; ++i ){
			gsd_utf8_char *name_i = type_names + i*buff_size;
			char name_buffer[buff_size];
			int j = 0;
			while( name_i[j].c[0] ){
				name_buffer[j] = name_i[j].c[0];
				++j;
			}
			name_buffer[j] = '\0';
			std::string name = name_buffer;
			
			//std::cerr << "Name " << i+1 << " = " << name << "\n";
		}

		

		delete [] type_names;
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
			// Reconstruct the type as an integer rather than
			// the named type. I am not happy about this loss
			// of info, but ok.
			b.types[i] = type_ids[i] + 1;
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
	delete [] type_ids;
	
	return 0;
	
}

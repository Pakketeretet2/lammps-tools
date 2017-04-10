#include "block_data_writers.h"
#include "gsd/gsd.h"

#include <fstream>


void write_block_lammps_dump( const block_data &b, const std::string &fname )
{
	std::ofstream o( fname );
	write_block_lammps_dump( b, o );
}


void write_block_lammps_dump( const block_data &b, std::ostream &o )
{
	std::string header_line = "ITEM: ATOMS id type x y z";
	if( b.atom_style == atom_styles::MOLECULAR ){
		header_line = "ITEM: ATOMS id mol type x y z";
	}else if( b.atom_style == atom_styles::ATOMIC ){
		
	}else{
		std::cerr << "Unknown atom style!\n";
		std::terminate();
	}
	
	o << "ITEM: TIMESTEP\n" << b.tstep << "\nITEM: NUMBER OF ATOMS\n";
	o << b.N << "\nITEM: BOX BOUNDS " << b.boxline << "\n";
	o << b.xlo[0] << " " << b.xhi[0] << "\n";
	o << b.xlo[1] << " " << b.xhi[1] << "\n";
	o << b.xlo[2] << " " << b.xhi[2] << "\n";
	o << header_line;
	if( b.other_cols.size() > 0 ){
		for( std::size_t i = 0; i < b.other_cols.size(); ++i ){
			o << " " << b.other_cols[i].header;
		}
	}
	o << "\n";
	for( py_int i = 0; i < b.N; ++i ){
		if( b.atom_style == atom_styles::ATOMIC ){
			o << b.ids[i] << " " << b.types[i] << " " << b.x[i][0]
			  << " " << b.x[i][1] << " " << b.x[i][2];
		}else if( b.atom_style == atom_styles::MOLECULAR ){
			o << b.ids[i] << " " << b.mol[i] << " " << b.types[i]
			  << " " << b.x[i][0] << " " << b.x[i][1] << " "
			  << b.x[i][2];
		}
		for( py_int j = 0; j < b.other_cols.size(); ++j ){
			o << " " << b.other_cols[j].data[i];
		}
		o << "\n";
	}
}

void write_block_hoomd_gsd( const block_data &b, const std::string &fname )
{
	gsd_handle gh;
	uint32_t schema_version = gsd_make_version( 1, 1 );
	
	int status = gsd_create_and_open( &gh, fname.c_str(), "lammps-tools",
	                                  "hoomd", schema_version,
	                                  GSD_OPEN_READWRITE, 0 );

	if( status ){

		std::cerr << "An error occured creating GSD file!\n";
	}else{
		write_block_hoomd_gsd( b, &gh );
	}

	gsd_close( &gh );
}

void write_block_hoomd_gsd( const block_data &b, gsd_handle *gh )
{
	int status;
	uint64_t step = b.tstep;
	uint8_t  dims = 3;
	uint32_t N    = b.N;
	double L[3];
	float box[6];
	L[0] = b.xhi[0] - b.xlo[0];
	L[1] = b.xhi[1] - b.xlo[1];
	L[2] = b.xhi[2] - b.xlo[2];
	box[3] = box[4] = box[5] = 0.0;

	box[0] = L[0];
	box[1] = L[1];
	box[2] = L[2];


	status = gsd_write_chunk( gh, "configuration/step", GSD_TYPE_UINT64,
	                          1, 1, 0, &step );
	std::cerr << "( " << status << " ): Wrote step to file.\n";
	status = gsd_write_chunk( gh, "configuration/dimensions", GSD_TYPE_UINT8,
	                          1, 1, 0, &dims );
	std::cerr << "( " << status << " ): Wrote dims to file.\n";
	status = gsd_write_chunk( gh, "configuration/box", GSD_TYPE_FLOAT,
	                          6, 1, 0, box );
	std::cerr << "( " << status << " ): Wrote box to file.\n";
	status = gsd_write_chunk( gh, "particles/N", GSD_TYPE_UINT32,
	                          1, 1, 0, &N );
	std::cerr << "( " << status << " ): Wrote N to file.\n";


	// float    *xx     = new float[3*N];
	float    *x     = new float[3*N];
	
	uint32_t *types	= new uint32_t[N];

	
	
	for( int i = 0; i < N; ++i ){
		// Sort them along id:
		int j = b.ids[i] - 1;

		// Remap the positions to -0.5L and 0.5L.
		double xi[3];
		for( int d = 0; d < 3; ++d ){
			xi[d] = b.x[i][d] - b.xlo[d];
			xi[d] -= 0.5*L[d];

			// Check box bounds:
			if( xi[d] > 0.5*L[d] || xi[d] < -0.5*L[d] ){
				std::cerr << "Particle " << b.ids[i] << " is "
				          << "out of box bound in dim " << d
				          << " ( " << -0.5*L[d] << ", "
				          << 0.5*L[d] << " ) with x = "
				          << xi[d] << ".\n";
				std::terminate();
			}
			
			x[3*j+d] = xi[d];
		}
		types[j] = b.types[i];
	}
	status = gsd_write_chunk( gh, "particles/position", GSD_TYPE_FLOAT,
	                          N, 3, 0, x );
	std::cerr << "( " << status << " ): Wrote x to file.\n";
	status = gsd_write_chunk( gh, "particles/typeid", GSD_TYPE_UINT32,
	                          N, 1, 0, types );
	std::cerr << "( " << status << " ): Wrote types to file.\n";

	status = gsd_end_frame( gh );

	delete [] x;
	delete [] types;

}


extern "C" {
	
void write_block_to_file( const block_data *bh, const char *fname,
                          const char *fformat, const char *dformat )
{
	std::string data_format( dformat );
	std::string file_format( fformat );

	std::cerr << "Writing block data at " << bh << " to " << fname << "\n";
	std::cerr << "File format is " << file_format << " and Data format is "
	          << data_format << ".\n";

	if( data_format == "LAMMPS" ){
		if( file_format == "GZIP" ){
			std::cerr << "GZIP not supported for LAMMPS yet!\n";
		}else if( file_format == "BIN" ){
			std::cerr << "BIN not supported for LAMMPS yet!\n";
		}else{
			write_block_lammps_dump( *bh, fname );
		}
		
	}else if( data_format == "HOOMD" ){
		if( file_format == "GZIP" || file_format == "PLAIN" ){
			std::cerr << file_format
			          << " not supported for HOOMD data!\n";
		}else{
			write_block_hoomd_gsd( *bh, fname );
		}
	}else{
		std::cerr << "Data format " << data_format
		          << " not recognized!\n";
	}
	
}

} // extern "C"

#include "block_data_writers.h"
#include "gsd/gsd.h"
#include "dump_interpreter_gsd.h"
#include "id_map.h"

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


void write_block_lammps_data( const block_data &b, const std::string &fname )
{
	std::ofstream o(fname);
	write_block_lammps_data( b, o );
		
}

void write_block_lammps_data( const block_data &b, std::ostream &o )
{
	o << "LAMMPS data file via lammps-tools\n\n";
	o << b.N << " atoms\n";
	py_int n_types = *(std::max_element( b.types, b.types + b.N ));
	o << n_types << " atom types\n\n";
	const char *words[3][2] = { { "xlo", "xhi"}, {"ylo", "yhi"}, {"zlo", "zhi"} };
	for( int dim = 0; dim < 3; ++dim ){
		o << b.xlo[dim] << " " << b.xhi[dim]
		  << " " << words[dim][0] << " " << words[dim][1] << "\n";
	}
	o << "\n";
	
	std::string atom_style = "atomic";
	if( b.atom_style == MOLECULAR ) atom_style = "molecular";

	o << "Atoms # " << atom_style << "\n\n";
	for( py_int i = 0; i < b.N; ++i ){
		o << b.ids[i];
		if( b.atom_style == MOLECULAR ){
			o << " " << b.mol[i];
		}
		o << " " << b.types[i] << " " << b.x[i][0]
		  << " " << b.x[i][1] << " " << b.x[i][2];

		// If you ever support image flags or something, add here.
		
		
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
	//std::cerr << "( " << status << " ): Wrote step to file.\n";
	status = gsd_write_chunk( gh, "configuration/dimensions", GSD_TYPE_UINT8,
	                          1, 1, 0, &dims );
	//std::cerr << "( " << status << " ): Wrote dims to file.\n";
	status = gsd_write_chunk( gh, "configuration/box", GSD_TYPE_FLOAT,
	                          6, 1, 0, box );
	//std::cerr << "( " << status << " ): Wrote box to file.\n";
	status = gsd_write_chunk( gh, "particles/N", GSD_TYPE_UINT32,
	                          1, 1, 0, &N );
	//std::cerr << "( " << status << " ): Wrote N to file.\n";


	// float    *xx     = new float[3*N];
	float    *x     = new float[3*N];
	uint32_t *types	= new uint32_t[N];
	int n_types = 0;
	
	
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
		int current_type = b.types[i];
		types[j] = current_type - 1;
		if( current_type > n_types ) n_types = current_type;
		
	}
	status = gsd_write_chunk( gh, "particles/position", GSD_TYPE_FLOAT,
	                          N, 3, 0, x );
	//std::cerr << "( " << status << " ): Wrote x to file.\n";
	status = gsd_write_chunk( gh, "particles/typeid", GSD_TYPE_UINT32,
	                          N, 1, 0, types );
	//std::cerr << "( " << status << " ): Wrote type ids to file.\n";
	
	// Write the actual type names.
	const uint buff_size = dump_interpreter_gsd::TYPE_BUFFER_SIZE;
	gsd_utf8_char *type_names = new gsd_utf8_char[buff_size * n_types]();
	//std::cerr << "gsd utf-8 char is " << sizeof(gsd_utf8_char) << " bytes.\n";
	
	for( int t = 0; t < n_types; ++t ){
		gsd_utf8_char * current_name = type_names + t*buff_size;
		std::string number = std::to_string( t+1 );
		int idx;
		for( idx = 0; idx < number.length(); ++idx ){
			current_name[idx].c[0] = number[idx];
		}
		current_name[idx].c[0] = '\0';
		/*
		std::cerr << "Setting name " << t+1 << " to ";
		for( int i = 0; i < idx; ++i ){
			std::cerr << current_name[i].c[0];
		}
		std::cerr << "\n";
		*/
	}
	status = gsd_write_chunk( gh, "particles/types", GSD_TYPE_UINT8,
	                          n_types, buff_size, 0, type_names );
	status = gsd_end_frame( gh );

	delete [] x;
	delete [] types;
	delete [] type_names;
}

void read_block_lammps_data_body( std::istream &in, std::string &line,
                                  block_data &b )
{
	std::vector<std::string> words = split(line);
	
	bool warned_image_flags = false;
	bool warned_velocities = false;
	bool warned_masses = false;
	bool warned_bonds = false;
	bool warned_angles = false;
	bool warned_dihedrals = false;
	bool warned_impropers = false;

	if( words[0] == "Atoms" ){
		if( words.size() > 3 && words[1] == "#" ){
			// Contains descriptor of atom_style:
			if( words[2] == "atomic" ){
				b.atom_style = ATOMIC;
			}else if( words[2] == "molecular" ){
				b.atom_style = MOLECULAR;
			}else{
				std::cerr << "WARNING: atom style " << words[2]
				          << " not known! Assuming atomic!\n";
				b.atom_style = ATOMIC;
			}
		}
		int n_expect = 5 + 1*(b.atom_style == MOLECULAR);
		std::getline( in, line );
		for( py_int i = 0; i < b.N; ++i ){
			std::getline( in, line );
			words = split(line);
			int n_entries = words.size();
			if( n_entries != n_expect &&
			    n_entries != n_expect + 3 ){
				std::cerr << "Unexpected number of columns!\n";
				std::cerr << "i = " << i << ", b.N = " << b.N << "\n";
				std::terminate();
			}

			std::stringstream ss(line);
			ss >> b.ids[i];
			if( b.atom_style == MOLECULAR ){
				ss >> b.mol[i];
			}
			ss >> b.types[i];
			ss >> b.x[i][0] >> b.x[i][1] >> b.x[i][2];
			std::cerr << "x[i] = " << b.x[i][0] << " " << b.x[i][1]
			          << " " << b.x[i][2] << "\n";
			if( !warned_image_flags && (n_entries == n_expect+3) ){
				std::cerr << "WARNING: Ignoring image flags.\n";
				warned_image_flags = true;
			}
		}
	}else if( words[0] == "Velocities" ){
		int n_expect = 4;
		std::getline( in, line );
		id_map im( b.ids, b.N );

		// Gotta add the other cols now.
		dump_col vx( "vx", b.N );
		dump_col vy( "vy", b.N );
		dump_col vz( "vz", b.N );
		
		
		
		for( py_int i = 0; i < b.N; ++i ){
			std::getline( in, line );
			std::cerr << "line: " << line << "\n";
			words = split(line);
			int n_entries = words.size();
			if( n_entries != n_expect ){
				std::cerr << "Unexpected number of columns!\n";
				std::cerr << "i = " << i << ", b.N = " << b.N << "\n";
				std::terminate();
			}

			std::stringstream ss(line);
			py_int idi;
			ss >> idi;
			py_int idx = im[idi];
			ss >> vx.data[idx] >> vy.data[idx] >> vz.data[idx];
		}

	}else if( words[0] == "Masses" ){
		std::getline( in, line );
		std::getline( in, line );
		for( int i = 0; i < b.Ntypes; ++i ){
			
		}
	}else if( words[0] == "Bonds" ){
		std::cerr << "WARNING: Ignoring bonds.\n";
		warned_velocities = true;
	}else if( words[0] == "Impropers" ){
		std::cerr << "WARNING: Ignoring impropers.\n";
		warned_velocities = true;
	}else if( words[0] == "Dihedrals" ){
		std::cerr << "WARNING: Ignoring dihedrals.\n";
		warned_velocities = true;
	}else if( words[0] == "Angles" ){
		std::cerr << "WARNING: Ignoring angles.\n";
		warned_velocities = true;
	}else{
		std::cerr << "Word " << words[0] << " not recognized!\n";
		std::terminate();
	}

	while( std::getline(in,line) ){
		if( !line.empty() ){
			std::cerr << "Got line " << line << " before recursion.\n";
			read_block_lammps_data_body( in, line, b );
		}
	}
}

block_data read_block_lammps_data( const std::string &fname )
{
	std::ifstream in(fname);
	std::string line = "";

	// Ignore first two lines:
	std::getline(in,line);
	std::getline(in,line);

	// Read the header first.
	py_int N = 0;
	py_int Ntypes = 0;
	py_int N_angles = 0;
	py_int N_angle_types = 0;
	py_int N_bonds = 0;
	py_int N_bond_types = 0;

	py_float xlo[3], xhi[3];

	std::vector<std::string> body_headers =
		{ "Atoms", "Velocities", "Masses", "Bonds",
		  "Angles", "Impropers", "Dihedrals" };

	block_data b;
	
	while( std::getline(in,line) ){
		if( line.empty() ) continue;
		
		std::vector<std::string> words = split( line );
		if( (words.size() < 2) &&
		    std::find( body_headers.begin(), body_headers.end(),
		               words[0] ) != body_headers.end() ){
			b.top.N_bonds  = N_bonds;
			b.top.N_angles = N_angles;
			
			b.init_topology();
			std::cerr << "Hit " << words[0] << ", gonna read body.\n";
			
			read_block_lammps_data_body( in, line, b );
			break;
		}
		
		if( words[1] == "atoms" ){
			N = std::stoi( words[0] );
			b.init(N);
		}else if( words[1] == "atom" && words[2] == "types" ){
			Ntypes = std::stoi( words[0] );
			b.init_per_type_arrays(Ntypes);
		}else if( words[1] == "bonds" ){
			N_bonds = std::stoi( words[0] );
			if( b.atom_style == ATOMIC ){
				std::cerr << "Bonds encountered for "
					"atomic data!\n";
			}
			
		}else if( words[1] == "bond" && words[2] == "types" ){
			N_bond_types = std::stoi( words[0] );
			if( b.atom_style == ATOMIC ){
				std::cerr << "Bonds encountered for "
					"atomic data!\n";
			}
			
		}else if( words[1] == "angles" ){
			N_angles = std::stoi( words[0] );
			if( b.atom_style == ATOMIC ){
				std::cerr << "Angles encountered for "
					"atomic data!\n";
			}
						
		}else if( words[1] == "angle" && words[2] == "types" ){
			N_angle_types = std::stoi( words[0] );
			if( b.atom_style == ATOMIC ){
				std::cerr << "Angles encountered for "
					"atomic data!\n";
			}
			
		}else if( words.size() > 3 ){
			if( words[2] == "xlo" && words[3] == "xhi" ){
				xlo[0] = std::stof( words[0] );
				xhi[0] = std::stof( words[1] );
			}
			if( words[2] == "ylo" && words[3] == "yhi" ){
				xlo[1] = std::stof( words[0] );
				xhi[1] = std::stof( words[1] );
			}
			if( words[2] == "zlo" && words[3] == "zhi" ){
				xlo[2] = std::stof( words[0] );
				xhi[2] = std::stof( words[1] );
			}
			
		}else{
			std::cerr << "Unrecognized line \"" << line
			          << "\" encountered in header!\n";
		}
	}
	
	
	return b;
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

	
	if( data_format == "LAMMPS" || data_format == "LAMMPS_DUMP" ){
		if( file_format == "GZIP" ){
			std::cerr << "GZIP not supported for LAMMPS yet!\n";
		}else if( file_format == "BIN" ){
			std::cerr << "BIN not supported for LAMMPS yet!\n";
		}else{
			if( std::string(fname) == "-" ){
				write_block_lammps_dump( *bh, std::cout );
			} else {
				write_block_lammps_dump( *bh, fname );
			}
		}
	}else if( data_format == "LAMMPS_DATA" ){
		if( file_format == "GZIP" ){
			std::cerr << "GZIP not supported for LAMMPS yet!\n";
		}else if( file_format == "BIN" ){
			std::cerr << "BIN not supported for LAMMPS yet!\n";
		}else{
			if( std::string(fname) == "-" ){
				write_block_lammps_data( *bh, std::cout );
			} else {
				write_block_lammps_dump( *bh, fname );
			}
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


py_int read_block_from_data_file( block_data *bh, const char *fname,
                                  const char *dformat )
{
	std::string dump_format( dformat );
	if( dump_format != "LAMMPS_DATA" ){
		std::cerr << "Don't know how to read in format "
		          << dformat << "\n";
		return 0;
	}
	block_data tmp = read_block_lammps_data( fname );
	// Copy the temporary over the pointed one.
	// bh->resize( tmp.N );
	*bh = tmp;
	
	return bh->N;
}

} // extern "C"

#include "dump_interpreter_lammps.h"
#include "util.h"

#include <fstream>


#ifdef VERBOSE_LIB
#define MY_CERR std::cerr << __FILE__ << "(" << __LINE__ << "): "
#else
static std::ofstream gobble("/dev/null");
#define MY_CERR gobble
#endif


dump_interpreter_lammps::dump_interpreter_lammps( const std::string &dname,
                                                  int file_type )
	: dump_interpreter(dname), r(nullptr), id_idx(-1), type_idx(-1),
	  x_idx(-1), y_idx(-1), z_idx(-1), mol_idx(-1),
	  n_cols(0), atom_style( atom_styles::ATOMIC ),
	  scaled(0), last_line("")
{
	if( file_type == dump_reader::PLAIN ){
		r = new text_reader_plain( dname );
	}else if( file_type == dump_reader::GZIP ){
		r = new text_reader_gzip( dname );
	}
	if( !r ){
		std::cerr << "Failed to set up text reader for file " << dname
		          << ", file format " << file_type << "!\n";
		std::terminate();
	}
}

dump_interpreter_lammps::~dump_interpreter_lammps()
{
	if( r ) delete r;
}


int dump_interpreter_lammps::next_block_meta( block_data &block )
{
	std::string line;

	while( r->getline( line ) ){
		MY_CERR << "Line is \"" << line << "\".\n";
		if( line == "ITEM: TIMESTEP" ){
			r->getline( line );
			last_meta.tstep = std::stoul( line );
		}else if( line == "ITEM: NUMBER OF ATOMS" ){
			r->getline( line );
			last_meta.N = std::stoul( line );
			MY_CERR << "N = " << last_meta.N << "\n";
		}else if( starts_with( line, "ITEM: BOX BOUNDS " ) ){
			last_meta.boxline = line;
			r->getline( line );
			std::stringstream dims( line );
			dims >> last_meta.xlo[0] >> last_meta.xhi[0];
			
			r->getline( line );
			dims.str(""); dims.clear();
			dims.str( line );
			dims >> last_meta.xlo[1] >> last_meta.xhi[1];

			r->getline( line );
			dims.str(""); dims.clear();
			dims.str( line );
			dims >> last_meta.xlo[2] >> last_meta.xhi[2];
		}else if( starts_with( line, "ITEM: ATOMS" ) ){
			// Stop there.
			last_line = line;
			return 0;
		}else{
			std::cerr << "Encountered line '" << line
			          << "' and have no clue what to do!\n";
			return -1;
		}
	}
	return -1;
}

int dump_interpreter_lammps::next_block_body( block_data &block )
{
	std::string line = last_line;
	
	block.copy_meta( last_meta );
	block.resize( block.N );

	MY_CERR << "Checking line \"" << line << "\"\n";
	
	if( starts_with( line, "ITEM: ATOMS" ) ){
		// Figure out which column maps which.
		// Read out the next block.N lines.
		if( headers.empty() ){
			if( line == "ITEM: ATOMS" ){
				dump_style = ATOMIC;
				set_headers( "" );
				n_cols = 5;
			}else{
				dump_style = CUSTOM;
				set_headers( line.substr( 12 ) );
			}
		}
		
		double Lx = block.xhi[0] - block.xlo[0];
		double Ly = block.xhi[1] - block.xlo[1];
		double Lz = block.xhi[2] - block.xlo[2];
		
		block.other_cols.resize( other_cols.size() );
		for( int i = 0; i < other_cols.size(); ++i ){
			block.other_cols[i].header = other_col_headers[i];
			// std::cerr << "Resizing other cols to " << block.N << "\n";
			block.other_cols[i].resize( block.N );
			// std::cerr << block.other_cols[i].header << " was the header->\n";
		}
		// std::cerr << "Grabbing " << block.N << " atoms.\n";
		for( int i = 0; i < block.N; ++i ){
			r->getline( line );
				
			std::stringstream ss( line );
			std::vector<std::string> words( n_cols );
			std::string w;
			std::size_t j = 0;
			while( ss >> w ){
				// std::cerr << "Word = " << w << "\n";
				words[j] = w;
				++j;
			}
			
			block.ids[i]   = std::stoul( words[id_idx] );
			block.types[i] = std::stoul( words[type_idx] );
			block.x[i][0]  = std::stof( words[x_idx] );
			block.x[i][1]  = std::stof( words[y_idx] );
			block.x[i][2]  = std::stof( words[z_idx] );
			block.atom_style = atom_style;
			
			
			if( is_bit<BIT_X>(scaled) ){
				block.x[i][0] *= Lx;
				block.x[i][0] += block.xlo[0];
			}
			
			if( is_bit<BIT_Y>(scaled) ){
				block.x[i][1] *= Ly;
				block.x[i][1] += block.xlo[1];
			}
			
			if( is_bit<BIT_Z>( scaled ) ){
				block.x[i][2] *= Lz;
				block.x[i][2] += block.xlo[2];
			}
			
			if( atom_style == atom_styles::MOLECULAR ){
				block.mol[i] = std::stoul( words[mol_idx] );
				block.atom_style = atom_style;
			}

			
			std::size_t k = 0;
			for( int j : other_cols ){
				double val = std::stof( words[j] );
				block.other_cols[k].data[i] = val;
			}
		}

		// std::cerr << "Block has atom_style " << block.atom_style << ".\n";

		
		return 0;
		
	}else{
		// std::cerr << "Encountered unknown header!\n";
		// std::cerr << line << "\n";
		return -1;
	}
}

       



int dump_interpreter_lammps::next_block( block_data &block )
{
	if( r->peek() == EOF ){
		if( r ){
			// Assume at EOF?
			return 1;
		}
	}
	int status = next_block_meta( block );
	if( !status ){
		MY_CERR << "Reading in block of size " << last_meta.N << "\n";
		return next_block_body( block );		
	}else{
		std::cerr << "Failed to get meta!\n";
		return -1;
	}
}

void dump_interpreter_lammps::set_headers( const std::string &h_line )
{
	// std::cerr << "Figuring out atom columns...\n";
	// std::cerr << "Headers: " << h_line << "\n";
	if( h_line.empty() ){
		// This happens if the dump style was atomic. In that case
		// the headers are preset.
		std::cerr << "Headers missing, assuming file was dump atom.\n";
		id_idx = 0;
		type_idx = 1;
		x_idx = 2;
		y_idx = 3;
		z_idx = 4;

		atom_style = atom_styles::ATOMIC;
		last_meta.atom_style = atom_style;
		return;
	}
	
	std::stringstream h( h_line );
	int header_idx = 0;

	if( !other_col_headers.empty() ){
		// std::cerr << "This is strange...\n";
	}
	
	do{
		std::string hh;
		if( h >> hh ){
			n_cols++;
			headers.push_back( hh );
			if( hh == "id" ){
				id_idx = header_idx;
			}else if( hh == "mol" ){
				mol_idx = header_idx;
			}else if( hh == "type" ){
				type_idx = header_idx;
			}else if( starts_with( hh, "x" ) ){
				x_idx = header_idx;
				if( hh == "xs" ){
					set_bit< BIT_X >( scaled );
				}else{
					clear_bit< BIT_X >( scaled );
				}
			}else if( starts_with( hh, "y" ) ){
				y_idx = header_idx;
				if( hh == "ys" ){
					set_bit< BIT_Y >( scaled );
				}else{
					clear_bit< BIT_Y >( scaled );
				}
			}else if( starts_with( hh, "z" ) ){
				z_idx = header_idx;
				if( hh == "zs" ){
					set_bit< BIT_Z >( scaled );
				}else{
					clear_bit< BIT_Z >( scaled );
				}
			}else{
				other_col_headers.push_back( hh );
				other_cols.push_back( header_idx );
				// std::cerr << "Got other col: " << hh
				//           << ", at index " << header_idx << "\n";
			}
			++header_idx;
		}
	}while( h );

	if( id_idx < 0 || type_idx < 0 || x_idx < 0 || y_idx < 0 || z_idx < 0 ){
		std::cerr << "Column mapping failed!\n";
		std::cerr << "id mol type x y z = " << id_idx << " " << mol_idx
		          << " " << type_idx << " " << x_idx << " " << y_idx
		          << " " << z_idx << "\n";
		std::terminate();
	}
	if( mol_idx > 0 ){
		atom_style = atom_styles::MOLECULAR;
	}else{
		atom_style = atom_styles::ATOMIC;
	}

	// std::cerr << "Other cols:";
	for( int i = 0; i < other_cols.size(); ++i ){
		// std::cerr << " " << other_col_headers[i] << "( "
		//           << other_cols[i] << " )";
	}
	// std::cerr << "\n";
}


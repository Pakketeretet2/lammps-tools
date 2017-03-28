#include "dump_interpreter_lammps.h"


bool dump_interpreter_lammps::next_block( reader_core *r, block_data &block )
{
	bool success = false;
	std::string line;

	py_int N     = 0;
	py_int tstep = 0;
	std::string boxline = "";
	py_float xlo[3], xhi[3];
	
	while( r->getline( line ) ){
		if( line == "ITEM: TIMESTEP" ){
			r->getline( line );
			tstep = std::stoul( line );
		}else if( line == "ITEM: NUMBER OF ATOMS" ){
			r->getline( line );
			N = std::stoul( line );
		}else if( starts_with( line, "ITEM: BOX BOUNDS " ) ){
			boxline = line.substr( 17 );
			
			r->getline( line );
			std::stringstream dims( line );
			dims >> xlo[0] >> xhi[0];
			
			r->getline( line );
			dims.str(""); dims.clear();
			dims.str( line );
			dims >> xlo[1] >> xhi[1];

			r->getline( line );
			dims.str(""); dims.clear();
			dims.str( line );
			dims >> xlo[2] >> xhi[2];

		}else if( starts_with( line, "ITEM: ATOMS " ) ){
			// Figure out which column maps which.
			// Read out the next b.N lines.
			block_data b(N);
			if( headers.empty() ){
				set_headers( line.substr( 12 ) );
			}
			
			
			b.boxline = boxline;
			b.tstep = tstep;

			b.xlo[0] = xlo[0];
			b.xlo[1] = xlo[1];
			b.xlo[2] = xlo[2];
			
			b.xhi[0] = xhi[0];
			b.xhi[1] = xhi[1];
			b.xhi[2] = xhi[2];

			double Lx = b.xhi[0] - b.xlo[0];
			double Ly = b.xhi[1] - b.xlo[1];
			double Lz = b.xhi[2] - b.xlo[2];
			

			b.other_cols.resize( other_cols.size() );
			for( int i = 0; i < other_cols.size(); ++i ){
				b.other_cols[i].header = other_col_headers[i];
				// std::cerr << "Resizing other cols to " << b.N << "\n";
				b.other_cols[i].resize( b.N );
				// std::cerr << b.other_cols[i].header << " was the header.\n";
			}
			// std::cerr << "Grabbing " << b.N << " atoms.\n";
			for( int i = 0; i < b.N; ++i ){
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
				
				b.ids[i]   = std::stoul( words[id_idx] );
				b.types[i] = std::stoul( words[type_idx] );
				b.x[i][0]  = std::stof( words[x_idx] );
				b.x[i][1]  = std::stof( words[y_idx] );
				b.x[i][2]  = std::stof( words[z_idx] );

				if( scaled | BIT_X ){
					b.x[i][0] *= Lx;
					b.x[i][0] += b.xlo[0];
				}
				
				if( scaled | BIT_Y ){
					b.x[i][1] *= Ly;
					b.x[i][1] += b.xlo[1];
				}
				
				if( scaled | BIT_Z ){
					b.x[i][2] *= Lz;
					b.x[i][2] += b.xlo[2];
				}

				if( atom_style == block_data::MOLECULAR ){
					b.mol[i] = std::stoul( words[mol_idx] );
				}

				std::size_t k = 0;
				for( int j : other_cols ){
					double val = std::stof( words[j] );
					b.other_cols[k].data[i] = val;
				}
			}
			// At this point, b contains everything you want.
			// Copy it to block and return.
			// std::cerr << "Done with reading block.\n";
			copy( block, b );
			
			return true;
			
		}else{
			// std::cerr << "Encountered unknown header!\n";
			// std::cerr << line << "\n";
			return false;
		}
	}
	return false;
}

void dump_interpreter_lammps::set_headers( const std::string &h_line )
{
	// std::cerr << "Figuring out atom columns...\n";
	// std::cerr << "Headers: " << h_line << "\n";

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
					scaled += BIT_X;
				}
			}else if( starts_with( hh, "y" ) ){
				y_idx = header_idx;
				if( hh == "ys" ){
					scaled += BIT_Y;
				}
			}else if( starts_with( hh, "z" ) ){
				z_idx = header_idx;
				if( hh == "zs" ){
					scaled += BIT_Z;
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
		// std::cerr << "Column mapping failed!\n";
		// std::cerr << "id mol type x y z = " << id_idx << " " << mol_idx
		//           << " " << type_idx << " " << x_idx << " " << y_idx
		//           << " " << z_idx << "\n";
		std::terminate();
	}
	if( mol_idx > 0 ){
		atom_style = block_data::MOLECULAR;
	}else{
		atom_style = block_data::ATOMIC;
	}

	// std::cerr << "Other cols:";
	for( int i = 0; i < other_cols.size(); ++i ){
		// std::cerr << " " << other_col_headers[i] << "( "
		//           << other_cols[i] << " )";
	}
	// std::cerr << "\n";
}


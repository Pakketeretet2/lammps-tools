#include "block_data_writers.h"



void write_block_lammps_dump( const block_data &b, std::ostream &o )
{
	std::string header_line = "ITEM: ATOMS id type x y z";
	if( b.atom_style == block_data::MOLECULAR ){
		header_line = "ITEM: ATOMS id mol type x y z";
	}else if( b.atom_style == block_data::ATOMIC ){
		
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
		if( b.atom_style == block_data::ATOMIC ){
			o << b.ids[i] << " " << b.types[i] << " " << b.x[i][0]
			  << " " << b.x[i][1] << " " << b.x[i][2];
		}else if( b.atom_style == block_data::MOLECULAR ){
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
	
}

void write_block_hoomd_gsd( const block_data &b, gsd_handle *h )
{
	
}

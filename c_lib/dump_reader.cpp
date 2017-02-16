#include "dump_reader.h"
#include "id_map.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>

#include <boost/iostreams/filter/gzip.hpp>

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

// static dummy_stream d_out;
//static std::ofstream d_out( "/dev/stderr" );
static debug_stream d_out( std::cerr );


block_data::block_data() : N(-1), idt("id"), typet("type"), xt("x"), yt("y"),
                           zt("z"), N_other_cols(-1)
{
}


block_data::~block_data()
{
}


// Copy constructor to make sure all members are properly set up.
block_data &block_data::operator=( const block_data &b )
{
	N = b.N;
	N_other_cols = b.N_other_cols;

	copy_meta( b );
	
	alloc( N );

	// After allocing, copy all data:
	std::copy( b.columns, b.columns + 5, columns );
	
	std::copy( b.x.begin(),     b.x.end()    , x.begin()     );
	std::copy( b.ids.begin(),   b.ids.end()  , ids.begin()   );
	std::copy( b.types.begin(), b.types.end(), types.begin() );
	
	std::copy( b.other_col_headers.begin(), b.other_col_headers.end(),
	           other_col_headers.begin() );
	std::copy( b.other_cols.begin(), b.other_cols.end(),
	           other_cols.begin() );
	if( ocols.size() > 0 ){
		std::copy( b.ocols.begin(), b.ocols.end(), ocols.begin() );
	}
	return *this;
}

int block_data::alloc( int size )
{
	d_out << "Allocking for size " << size << "...\n";
	x.resize( size );
	ids.resize( size );
	types.resize( size );
	other_cols.resize( N_other_cols );
	other_col_headers.resize( N_other_cols );
	ocols.resize(N_other_cols);
	
	for( int i = 0; i < N_other_cols; ++i ){
		other_cols[i].resize(size);
	}

	

	return 0;
}
	
void block_data::clear()
{
	
}

double block_data::memory()
{
	return N*6*sizeof(double);
}

void block_data::word_to_array( const std::string &w, int idx, int i )
{
	if( idx == columns[0] )      ids[i]   = std::stod( w );
	else if( idx == columns[1] ) types[i] = std::stod( w );
	else if( idx == columns[2] ) x[i][0]  = std::stod( w );
	else if( idx == columns[3] ) x[i][1]  = std::stod( w );
	else if( idx == columns[4] ) x[i][2]  = std::stod( w );
	else{
		for( int j = 0; j < N_other_cols; ++j ){
			if( idx == ocols[j] ){
				other_cols[j][i] = std::stod( w );
			}
		}
	}
				
}


void block_data::set_xyz_tags( const char *tx, const char *ty, const char *tz )
{
	xt = tx;
	yt = ty;
	zt = tz;
	
}

void block_data::set_id_tag( const char *tid )
{
	idt = tid;
}

void block_data::set_type_tag( const char *ttype )
{
	typet = ttype;
}

void block_data::set_column_indices( const std::string &line )
{
	// Read out the column headers and map the indices
	std::stringstream ss(line);
	std::string word;
	
	int idx = 0;
	int ocol_idx = 0;
	N_other_cols = 0;
	
	// Format of line: ITEM: ATOMS ...
	ss >> word >> word;

	std::vector<std::string> col_headers;
	
	while( ss >> word ){
		col_headers.push_back(word);

		d_out << "Trying to match word \"" << word << "\"...";
		
		if( word == idt ){
			columns[0] = idx;
			d_out << " Matched id (" << idx << ")!\n";
		}else if( word == typet ){
			columns[1] = idx;
			d_out << " Matched type (" << idx << ")!\n";
		}else if( word == xt ){
			columns[2] = idx;
			d_out << " Matched x (" << idx << ")!\n";
		}else if( word == yt ){
			columns[3] = idx;
			d_out << " Matched y (" << idx << ")!\n";
		}else if( word == zt ){
			columns[4] = idx;
			d_out << " Matched z (" << idx << ")!\n";
		}else{
			d_out << " Was other col!\n";
			++N_other_cols;
		}
		++idx;
	}

	// Need a second run to assign the other cols properly:
	d_out << "Resizing ocols to " << N_other_cols << " to fit other cols.\n";
	ocols.resize(N_other_cols);
	other_col_headers.resize(N_other_cols);
	d_out << "ocols now has size " << ocols.size() << ".\n";
	
	idx = 0;
	for( const std::string &header : col_headers ){
		if( header != idt && header != typet && header != xt &&
		    header != yt && header != zt ){
			ocols[ocol_idx] = idx;
			other_col_headers[ocol_idx] = header;
			d_out << " Matched other col (" << idx << ") to "
			      << ocol_idx+1 << "!\n";
			++ocol_idx;
		}
		++idx;
	}
	d_out << "Other cols are: ";
	for( int i = 0; i < N_other_cols; ++i ){
		d_out << other_col_headers[i] << " ";
	}
	d_out << "\n";
	d_out.flush();
}


void block_data::init_empty( int size )
{
	N = size;
	N_other_cols = 0;

	alloc(size);
	
	
}

void block_data::copy_meta( const block_data &b )
{
	idt = b.idt;
	typet = b.typet;
	xt = b.xt;
	yt = b.yt;
	zt = b.zt;

	std::copy( b.xlo, b.xlo + 3, xlo );
	std::copy( b.xhi, b.xhi + 3, xhi );
	tstep = b.tstep;
}


dump_reader::dump_reader( const std::string &fname ) : at_eof(false),
                                                       infile_name(fname)
{
	init_infile(fname);
}

dump_reader::~dump_reader()
{
}

// Skips the next Nblocks blocks.
bool dump_reader::skip_blocks( int Nblocks )
{
	for( int i = 0; i < Nblocks; ++i ){
		bool success = skip_block();
		if( !success ){
			d_out << "Block skip failed, returning...\n";
			return false;
		}
	}
	return true;
}

bool dump_reader::skip_block()
{
	switch( file_format ){
#ifdef HAVE_BOOST_GZIP
		case GZIP:
			return skip_block_from_istream( infile_filt );
#endif
		case PLAIN:
		default:
			return skip_block_from_istream( infile );
	}	
	
	return true;
}



bool dump_reader::next_block( block_data &block )
{
	if( at_eof ) return false;

	switch( file_format ){
#ifdef HAVE_BOOST_GZIP
		case GZIP:
			return next_block_from_istream( block, infile_filt );
#endif
		case PLAIN:
		default:
			return next_block_from_istream( block, infile );
	}
}

bool dump_reader::last_block( block_data &block )
{
	if( at_eof ) return false;


	switch( file_format ){
#ifdef HAVE_BOOST_GZIP
		case GZIP:
			return last_block_from_istream( block, infile_filt );
#endif
		case PLAIN:
		default:
			return last_block_from_istream( block, infile );
	}
}



void dump_reader::init_infile( const std::string &fname )
{
	if( ends_with( fname, ".gz" ) ){
#ifdef HAVE_BOOST_GZIP
		infile = std::ifstream( fname, std::ios_base::in | std::ios_base::binary );
		file_format = GZIP;
		infile_filt.push(boost::iostreams::gzip_decompressor());
		infile_filt.push(infile);
#else
		assert(false && "BOOST GZIP LIBRARY NOT INSTALLED!");
#endif
	}else{
		file_format = PLAIN;
		infile = std::ifstream( fname, std::ios_base::in );
	}

	
}

bool dump_reader::getline( std::string &line, int mode )
{
	if( mode == 1 ){
		std::getline( std::cin, line );
		return false; // stuff;
	}
		
	switch( file_format ){
#ifdef HAVE_BOOST_GZIP
		case GZIP:
			if( std::getline( infile_filt, line ) ){
				return true;
			}else{
				return false;
			}
			break;
#endif
		case PLAIN:
		default:
			if( std::getline( infile, line ) ){
				return true;
			}else{
				return false;
			}
			break;
	}
}




void block_data::print( const std::string &out_name, std::ostream &err )
{
	std::ofstream tmp(out_name);
	print( tmp, err );
}


void block_data::print_data( std::ostream &out, std::ostream &err )
{
	long int max_type = *std::max_element( types.begin(), types.end() );
	out << "# LAMMPS data file made by dump_reader.\n\n";
	out << N << " atoms\n" << max_type << " atom types\n\n";
	out << xlo[0] << " " << xhi[0] << " xlo xhi\n";
	out << xlo[1] << " " << xhi[1] << " ylo yhi\n";
	out << xlo[2] << " " << xhi[2] << " zlo zhi\n\n";
	out << "Masses\n\n";
	for( long int i = 1; i <= max_type; ++i ){
		out << i << " " << i << "\n";
	}
	out << "\n";
	out << "Atoms\n\n";
	for( std::size_t i = 0; i < N; ++i ){
		out << ids[i] << " " << types[i] << " " << x[i][0]
		       << " " << x[i][1] << " " << x[i][2] << "\n";
	}
}


void block_data::print( std::ostream &out, std::ostream &err )
{
	if( N < 0 ){
		err << "Careful! Block seems empty!\n";
		return;
	}
	out << "ITEM: TIMESTEP\n";
	out << tstep << "\nITEM: NUMBER OF ATOMS\n" << N << "\n";
	out << "ITEM: BOX BOUNDS " << boxline << "\n";
	out << xlo[0] << " " << xhi[0] << "\n";
	out << xlo[1] << " " << xhi[1] << "\n";
	out << xlo[2] << " " << xhi[2] << "\n";
	out << "ITEM: ATOMS " << idt << " " << typet << " " << xt << " " << yt
	    << " " << zt;
	for( int i = 0; i < ocols.size(); ++i ){
		out << " " << other_col_headers[i];
	}
	out << "\n";
	for( int i = 0; i < N; ++i ){
		out << ids[i] << " " << types[i] << " " << x[i][0] << " "
		    << x[i][1] << " " << x[i][2];
		for( int j = 0; j < ocols.size(); ++j ){
			out << " " << other_cols[j][i];
		}
		out << "\n";
	}
}


void block_data::append_other_col( std::string header, std::vector<double> data )
{
	d_out << "ocols now has size " << ocols.size() << ".\n";
	int old_ncols = N_other_cols;
	other_col_headers.push_back( header );
	d_out << "Other headers: ";
	for( const std::string &s : other_col_headers ){
		d_out << s << " ";
	}
	d_out << "\n";
	other_cols.push_back( data );
	ocols.push_back( 5 + old_ncols );
	N_other_cols++;
	d_out << "After appending column, there are now " << ocols.size()
	      << " other columns.\n";
}

int block_data::get_column_index( const std::string &header )
{
	int i = 0;
	while( other_col_headers[i] != header ) ++i;
	return i;
}

block_data filter_block( block_data b, const std::vector<long int> &ids )
{
	std::list<long int> l(ids.begin(), ids.end());
	return filter_block( b, l );
}

block_data filter_block_by_indices( block_data b,
                                    const std::vector<long int> &indices )
{
	std::list<long int> l(indices.begin(), indices.end());
	return filter_block_by_indices( b, l );	
}



block_data filter_block_by_indices( block_data b,
                                    const std::list<long int> &indices )
{
	block_data b2 = b;

	// Copy only the relevant ids.
	int N = indices.size();
	b2.x.resize( N );
	b2.ids.resize( N );
	b2.types.resize( N );
	b2.N = N;
	
	id_map im( b.ids );
	int j = 0;
	for( const long int &i : indices ){
		b2.x[j]     = b.x[i];
		b2.ids[j]   = b.ids[i];
		b2.types[j] = b.types[i];

		for( int k = 0; k < b.N_other_cols; ++k ){
			b2.other_cols[k][j] = b.other_cols[k][i];
		}
		++j;
	}

	return b2;
}



block_data filter_block( block_data b, const std::list<long int> &ids )
{
	block_data b2 = b;
	std::ofstream out( "id.dat" );

	// Copy only the relevant ids.
	int N = ids.size();
	b2.x.resize( N );
	b2.ids.resize( N );
	b2.types.resize( N );
	b2.N = N;
	
	id_map im( b.ids );
	int j = 0;
	for( const long int &id : ids ){
		int i = im[ id ];
		b2.x[j]     = b.x[i];
		b2.ids[j]   = b.ids[i];
		b2.types[j] = b.types[i];

		out << "Got match! " << b.ids[i] << " == " << id << ".\n";

		for( int k = 0; k < b.N_other_cols; ++k ){
			b2.other_cols[k][j] = b.other_cols[k][i];
		}
		++j;
	}

	out << "Ids after filtering:\n";
	int i = 0;
	for( const long int &id : ids ){
		out << id << " " << b2.ids[i] << "\n";
		++i;
	}

	return b2;
}



std::size_t dump_reader::block_count()
{
	std::size_t n_blocks = 0;
	block_data b;
	while( next_block(b) ){
		n_blocks++;
		if( n_blocks % 100 == 0 ){
			std::cerr << "At block " << n_blocks << "...\n";
		}
	}
	
	at_eof = false;
	rewind();
	
	return n_blocks;
}



void dump_reader::rewind()
{
	if( infile ){
		infile.close();
	}
#ifdef HAVE_BOOST_GZIP
	if( infile_filt ){
		infile_filt.clear();
	}
#endif // HAVE_BOOST_GZIP
	
	init_infile( infile_name );
}



bool next_block_from_istream( block_data &block, std::istream &in )
{
	if( !in ) return false;

	std::string line;
	std::string garbage;

	block_data bb;
	block_data *b = nullptr;

	int tstep;
	int N;

	bool need_alloc = true;
	bool success = false;

	while( std::getline( in, line ) ){
		d_out << "line: \"" << line << "\".\n";
		d_out.flush();
		if( line == "ITEM: TIMESTEP" ){
			std::getline( in, line );
			tstep = std::stoi(line);
			d_out << "At t = " << tstep << "...\n";
		}else if( starts_with( line, "ITEM: BOX BOUNDS" ) ){
			std::stringstream ss( line );
			ss >> garbage >> garbage >> garbage;
			ss >> garbage;
			b->boxline = garbage;
			b->boxline += " ";
			ss >> garbage;
			b->boxline = garbage;
			b->boxline += " ";
			ss >> garbage;
			b->boxline = garbage;
			std::getline( in, line );
			ss.str(""); ss.clear(); ss.str(line);
			ss >> b->xlo[0] >> b->xhi[0];
			std::getline( in, line );
			ss.str(""); ss.clear(); ss.str(line);
			ss >> b->xlo[1] >> b->xhi[1];
			std::getline( in, line );
			ss.str(""); ss.clear(); ss.str(line);
			ss >> b->xlo[2] >> b->xhi[2];
			d_out << "xlo, xhi = [ " << b->xlo[0] << " "
			      << b->xlo[1] << " " << b->xlo[2] << " ] , [ "
			      << b->xhi[0] << " " << b->xhi[1] << " "
			      << b->xhi[2] << " ].\n";
		}else if( line == "ITEM: NUMBER OF ATOMS" ){
			std::getline( in, line );
			N = std::stoi( line );
			d_out << "N = " << N << ".\n";
			
			if( block.N == N ){
				// Nice, no need to alloc stuff!
				d_out << "No need for allocating...\n";
				need_alloc = false;
				b = &block;
			}else{
				d_out << "Allocating...\n;";
				b = &bb;
				b->init_empty( N );
			}
			b->tstep = tstep;
			b->N = N;
			
		}else if( starts_with( line, "ITEM: ATOMS" ) ){
			d_out << "Encountered atom data...\n";
			if( need_alloc ){
				b->set_column_indices( line );
				d_out << "Allocating data...";
				d_out.flush();
				b->alloc( b->N );
				d_out << " Done!\n";
			}
			
			
			for( int i = 0; i < b->N; ++i ){
				std::getline( in, line );
				std::stringstream ss(line);
				std::string word;
				int idx = 0;
				while( ss >> word ){
					b->word_to_array( word, idx, i );
					++idx;
				}
			}
			d_out << "Grabbed atom data for " << b->N << " atoms.\n";
			success = true;
			if( need_alloc ){
				// Do NOT do this if an alloc was not needed,
				// as then you worked directly on block.
				block = bb;
			}
			return success;
		}else{
			std::cerr << "Got this line: \"" << line << "\"... :/ \n";
			if( line.empty() ){
				// Probably at EOF:
				d_out << "At EOF!\n";
				return false;
			}else{
				d_out << "Failed to parse line! \n\""
				          << line << "\"\n";
				std::terminate();
			}
		}
	}

	return success;
}


bool last_block_from_istream( block_data &block, std::istream &in )
{
	bool success = true;
	while( success ){
		success = next_block_from_istream( block, in );
	}
	return success;
}


bool skip_block_from_istream( std::istream &in )
{
	std::string line;
	int natoms = 0;
	bool success = false;
	std::size_t t = 0;
	
	while( std::getline( in, line ) ){
		if( line == "ITEM: TIMESTEP" ){
			std::getline( in, line );
			t = std::stoi(line);
		}else if( starts_with( line, "ITEM: BOX BOUNDS" ) ){
			std::getline( in, line );
			std::getline( in, line );
			std::getline( in, line );

		}else if( line == "ITEM: NUMBER OF ATOMS" ){
			std::getline( in, line );
			natoms = std::stoi( line );
			d_out << "Skipping " << natoms << " atoms...\n";
		}else if( starts_with( line, "ITEM: ATOMS" ) ){
			for( int i = 0; i < natoms; ++i ){
				std::getline( in, line );
			}
			return true;
		}else{
			if( line.empty() ){
				// Probably at EOF:
				d_out << "At EOF!\n";
				return false;
			}else{
				d_out << "Failed to parse line! \n\""
				          << line << "\"\n";
				return false;
			}
		}
	}
	return success;
}


void block_data_from_foreign( void *xx, py_int N, py_int *iids,
                              py_int *ttypes, py_int periodic,
                              py_float *xlo, py_float *xhi, py_int dims,
                              py_int tstep, const char *box_line,
                              block_data &b )
{
	// std::cerr << "Allocking for " << N << " atoms.\n";
	b.init_empty(N);
	b.N = N;

	arr3f x( xx, N );
	arr1i ids( iids, N );
	arr1i types( ttypes, N );

	for( py_int i = 0; i < N; ++i ){
		b.x[i][0]  = x[i][0];
		b.x[i][1]  = x[i][1];
		b.x[i][2]  = x[i][2];
		b.ids[i]   = ids[i];
		b.types[i] = types[i];
	}

	b.tstep = tstep;
	std::copy( xlo, xlo + 3, b.xlo );
	std::copy( xhi, xhi + 3, b.xhi );
	b.boxline = box_line;
}


block_data block_data_from_data_file( const char *fname, int &status )
{
	block_data b;
	std::ifstream in(fname);
	std::string line;
	// Skip first two lines:
	std::getline( in, line );
	std::getline( in, line );

	int lc = 3;
	auto my_get_line = [&in, &line, &lc] ( bool silent = true ) -> std::istream& 
	{
		std::getline( in, line );
		++lc;
		if( !silent ){
			d_out << "Line " << lc << ": '" << line << "'\n";
		}
		return in;
		
	};
	int natoms;
	int natom_types;
	my_get_line();
	
	std::stringstream ss(line);

	ss >> natoms;
	ss.str("");
	ss.clear();

	my_get_line();
	ss.str(line);
	ss >> natom_types;
	d_out << "Got " << natom_types << " atom types...\n";
	ss.str("");
	ss.clear();

	b.N = natoms;

	my_get_line();
	my_get_line();

	// Box dimensions:
	ss.str("");
	ss.clear();
	ss.str(line);
	ss >> b.xlo[0] >> b.xhi[0];
	my_get_line();
	ss.str("");
	ss.clear();
	ss.str(line);
	ss >> b.xlo[1] >> b.xhi[1];
	my_get_line();
	ss.str("");
	ss.clear();
	ss.str(line);
	ss >> b.xlo[2] >> b.xhi[2];
	my_get_line();
	// my_get_line();
	
	// Now past meta, set up arrays and read the rest...
	b.init_empty(natoms);
	
	while( my_get_line() ){
		std::string header;
		line = rstrip( line, '#' );
		if( line == "Masses" ){
			std::cerr <<"Ignoring masses...\n";
			my_get_line(  );
			for( int i = 0; i < natom_types; ++i ){
				my_get_line(  );
			}
			my_get_line(  );
		}else if( line == "Pair Coeffs" ){
			std::cerr <<"Ignoring pair coeffs...\n";
			my_get_line(  );
			for( int i = 0; i < natom_types; ++i ){
				my_get_line(  );
			}
			my_get_line(  );
		}else if( line == "Atoms" ){
			std::cerr <<"Grabbing atoms...\n";
			my_get_line();
			for( int i = 0; i < natoms; ++i ){
				my_get_line( );
				ss.str("");
				ss.clear();
				ss.str( line );
				ss >> b.ids[i] >> b.types[i] >> b.x[i][0]
				   >> b.x[i][1] >> b.x[i][2];
			}
			my_get_line();

		}else if( line == "Velocities" ){
			std::cerr <<"Ignoring velocities...\n";
			my_get_line();
			for( int i = 0; i < natoms; ++i ){
				my_get_line();
			}
			my_get_line();
		}else{
			std::cerr << "Unrecognized header '" << line << "'.\n";
			status = -1;
			return b;
		}
	}
	
	status = 0;
	return b;
}

#ifndef DUMP_READER_H
#define DUMP_READER_H

/*!
  @file dump_reader.h
  @brief An interface for reading various dump files.

  \ingroup cpp_lib
*/



#include "util.h"
#include "block_data.h"

#include <iosfwd>

class reader_core
{
public:
	reader_core(){}
	virtual ~reader_core(){}
	
	virtual bool getline( std::string &line );
	virtual void rewind();
};


class dump_interpreter
{
public:
	dump_interpreter(){}
	virtual ~dump_interpreter(){}

	virtual bool next_block( reader_core *r, block_data &block );
	virtual bool last_block( reader_core *r, block_data &block );

	///< Calling next_block_meta only leaves the block in a half-
	///< half-initialised state! You can read out the meta but not the rest.
	virtual bool next_block_meta( reader_core *r, block_data &block );
	virtual bool next_block_body( reader_core *r, block_data &block );
	

	struct block_meta {
		py_int N, tstep;
		py_float xlo[3], xhi[3];
		std::string boxline;
	};

protected:
	std::string last_line;
	block_meta  last_meta;
};


class dump_reader
{
public:
	enum FILE_FORMATS {
		PLAIN   = 0,
		GZIP    = 1,
		ISTREAM = 2,
		BIN     = 3
	};
	
	enum DUMP_FORMATS {
		LAMMPS  = 0,
		GSD     = 1,
		DCD     = 2
	};

	static const char *fformat_to_str( int file_format )
	{
		switch( file_format ){
			default:
				return "UNKNOWN!";
			case PLAIN:
				return "PLAIN TEXT";
			case GZIP:
				return "GZIPPED TEXT";
			case ISTREAM:
				return "INPUT STREAM";
			case BIN:
				return "BINARY";
		}
	}
	static const char *dformat_to_str( int dformat )
	{
		switch( dformat ){
			default:
				return "UNKOWN!";
			case LAMMPS:
				return "LAMMPS";
			case GSD:
				return "GSD";
			case DCD:
				return "DCD";
		}
	}

	dump_reader( std::istream &stream, int dformat = LAMMPS );
	dump_reader( const std::string &fname );
	

	void setup_reader( std::istream &stream );
	void setup_reader( const std::string &fname, int format );

	void setup_interpreter( int dump_format );
	void setup_interpreter_gsd( const std::string &fname );

	/// Tries to read in the next block from file. If successful,
	/// returns true and block_data contains the next block.
	/// Else returns false and block is unchanged.
	virtual bool next_block( block_data &block );

	/// To return read the last block in the file. If successful,
	/// returns true and block_data contains the next block.
	/// Else returns false and block_data is unchanged.
	virtual bool last_block( block_data &block );
	
	virtual ~dump_reader();

	operator bool() const
	{
		return !at_eof;
	}

	/// Fast-forward a number of blocks.
	bool skip_blocks( int Nblocks );

	/// Fast-forward one block.
	bool skip_block ( );

	/// Returns the number of blocks in the file.
	std::size_t block_count();

	/// Rewinds the file to the beginning.
	
private:

	bool at_eof;
	int file_format;
	int dump_format;

	reader_core      *reader;
	dump_interpreter *interp;
};



// Use this to interface with the C++ dump reader.
extern "C" {

struct dump_reader_handle {
	dump_reader *reader;
	block_data  *last_block;
};
	
dump_reader_handle *get_dump_reader_handle( const char *dname );
void release_dump_reader_handle( dump_reader_handle *dh );
bool dump_reader_next_block( dump_reader_handle *dh );
void dump_reader_get_block_meta( dump_reader_handle *dh,
                                 py_int *tstep, py_int *N,
                                 py_float *xlo, py_float *xhi,
                                 py_int *periodic, char *boxline,
                                 py_int *atom_style );
void dump_reader_get_block_data( dump_reader_handle *dh,
                                 py_int N, py_float *x, py_int *ids,
                                 py_int *types, py_int *mol );

bool dump_reader_fast_forward( dump_reader_handle *dh,
                               py_int N_blocks );

	

} // extern "C"






#endif // DUMP_READER_H

#ifndef DUMP_READER_H
#define DUMP_READER_H

/*!
  @file dump_reader.h
  @brief An interface for reading various dump files.

  \ingroup cpp_lib
*/



#include "block_data.h"
#include "dump_interpreter.h"
#include "util.h"


#include <iosfwd>


/// A generic class for reading in dump files.
class dump_reader
{
public:
	/// Specifies various file formats
	enum FILE_FORMATS {
		PLAIN   = 0,
		GZIP    = 1,
		ISTREAM = 2,
		BIN     = 3
	};

	/// Specifies various dump formats
	enum DUMP_FORMATS {
		LAMMPS  = 0,
		GSD     = 1,
		NAMD    = 2
	};

	/**
	   Pretty-prints given file format.

	   \param   file_format The file format to pretty-print.

	   \returns A pretty-printed string of the file format.
	*/
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

	/**
	   Pretty-prints given dump format.

	   \param   dump_format The file format to pretty-print.

	   \returns A pretty-printed string of the dump format.
	*/
	static const char *dformat_to_str( int dformat )
	{
		switch( dformat ){
			default:
				return "UNKOWN!";
			case LAMMPS:
				return "LAMMPS";
			case GSD:
				return "GSD";
			case NAMD:
				return "NAMD";
		}
	}

	
	/// Construct dump reader from file.
	dump_reader( const std::string &fname );
	dump_reader( const std::string &fname, int dformat, int fformat );
	
	/// Cleanup:
	virtual ~dump_reader();

	/// Prepare dump interpreter:
	void setup_interpreter( const std::string &fname );

	/// Tries to read in the next block from file. If successful,
	/// returns true and block_data contains the next block.
	/// Else returns false and block is unchanged.
	int next_block( block_data &block );

	/// To return read the last block in the file. If successful,
	/// returns true and block_data contains the next block.
	/// Else returns false and block_data is unchanged.
        int last_block( block_data &block );

	/// Checks if the internal file is good or not.
	bool eof()  const { return interp->eof(); }
	bool good() const { return interp->good(); }
	

	/// Fast-forward a number of blocks.
        int skip_blocks( int Nblocks );

	/// Fast-forward one block.
	int skip_block();

	/// Rewinds to the beginning.
	void rewind();
	
	/// Returns the number of blocks in the file.
	std::size_t block_count();

	
private:
	int dump_format;  //!< Stores the dump format
	int file_format;  //!< Stores the file format

	dump_interpreter *interp; //!< Points to an internal dump_interpreter

	/// Guesses file type from name
	void guess_file_type ( const std::string &fname );
	/// Guesses dump type from name
	void guess_dump_type ( const std::string &fname );	
};




extern "C" {

/// This struct is used for interfacing between Python and the C++ lib.
struct dump_reader_handle {
	dump_reader_handle() : reader(nullptr), last_block(nullptr){}
	dump_reader *reader;
	block_data  *last_block;
};

/**
   Constructs a dump_reader_handle and passes it back.

   \param dname    The name of the dump file.
   \param dformat  The dump format, should be one of the DUMP_FORMATS
   \param fformat  The file format, should be one of the FILE_FORMATS
   
   Make sure the dump_reader_handle is freed using the matching
   release_dump_reader_handle.
*/
dump_reader_handle *get_dump_reader_handle( const char *dname,
                                            py_int dformat, py_int fformat );

/**
   Releases a dump_reader_handle ensuring proper clean-up.

   \param dh   The dump_reader_handle to free.
*/
void release_dump_reader_handle( dump_reader_handle *dh );

/**
   Makes the dump_reader_handle \p dh read in a new block and store
   it in its matching \p last_block field.

   \param dh  Ptr to the dump_reader_handle that is to read the file.
   
   \returns   An exist signal. 0 on success, negative on an unknown
              failure, positive if the end of file was encountered.
*/
int  dump_reader_next_block( dump_reader_handle *dh );

/**
   Grabs the meta info from the \p last_block field of the dump_reader_handle 
   and stores those in the passed fields.

   \param dh          Ptr to the dump_reader_handle that is to read the file.
   \param tstep       Time step of the block_data.
   \param N           Number of atoms
   \param xlo         Lower bounds of simulation box
   \param xhi         Upper bounds of simulation box
   \param periodic    Periodic bits
   \param boxline     Line describing the box
   \param atom_style  The atom stle of the block_data
*/
void dump_reader_get_block_meta( dump_reader_handle *dh,
                                 py_int *tstep, py_int *N,
                                 py_float *xlo, py_float *xhi,
                                 py_int *periodic, char *boxline,
                                 py_int *atom_style );
void dump_reader_get_block_data( dump_reader_handle *dh,
                                 py_int N, py_float *x, py_int *ids,
                                 py_int *types, py_int *mol );

int dump_reader_fast_forward( dump_reader_handle *dh,
                              py_int N_blocks );

	

} // extern "C"






#endif // DUMP_READER_H

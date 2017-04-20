#ifndef BLOCK_DATA_WRITERS_H
#define BLOCK_DATA_WRITERS_H

#include "block_data.h"
#include <iosfwd>

#include "gsd/gsd.h"

/**
   Writes block data in LAMMPS dump format to output stream.
   
   \param b The block_data to write.
   \param o The output stream to write \p b to.
*/
void write_block_lammps_dump( const block_data &b, std::ostream &o );


/**
   Writes block data in LAMMPS dump format to file
   
   \param b     The block_data to write.
   \param fname The file name to write to.
*/
void write_block_lammps_dump( const block_data &b, const std::string &fname );


/**
   Writes block data in LAMMPS data format to output stream.
   
   \param b     The block_data to write.
   \param fname The file name to write to.
*/
void write_block_lammps_data( const block_data &b, const std::string &fname );


/**
   Writes block data in LAMMPS data format to output stream.
   
   \param b   The block_data to write.
   \param o   The stream to write to.
*/
void write_block_lammps_data( const block_data &b, std::ostream &o );


/**
   Writes block data in HOOMD schema to file named fname.
   \param b     The block_data to write.
   \param fname The name of the GSD file to write \p b to.
*/
void write_block_hoomd_gsd( const block_data &b, const std::string &fname );


/**
   Writes block data in HOOMD schema to given gsd file.
   
   \param b     The block_data to write.
   \param gh    Handle to an open, writeable gsd file.
*/
void write_block_hoomd_gsd( const block_data &b, gsd_handle *gh );


/**
   Reads data file in LAMMPS format to a block_data object.
   
   \param fname   Name of the LAMMPS data file
   \returns       A new block_data object.
*/
block_data read_block_lammps_data( const std::string &fname );




/** C interface to write to files from Python: **/
extern "C" {
void write_block_to_file( const block_data *bh, const char *fname,
                          const char *fformat, const char *dformat );

py_int read_block_from_data_file( block_data *bh, const char *fname,
                                  const char *dformat );
	
	
} // extern "C";


#endif // BLOCK_DATA_WRITERS_H

#ifndef BLOCK_DATA_WRITERS_H
#define BLOCK_DATA_WRITERS_H

#include "block_data.h"
#include <iosfwd>

void write_block_lammps_dump( const block_data &b, std::ostream &o );

void write_block_hoomd_gsd( const block_data &b, const std::string &fname );
void write_block_hoomd_gsd( const block_data &b, gsd_handle *h );


#endif // BLOCK_DATA_WRITERS_H

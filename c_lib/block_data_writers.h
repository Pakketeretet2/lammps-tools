#ifndef BLOCK_DATA_WRITERS_H
#define BLOCK_DATA_WRITERS_H

#include "block_data.h"
#include <iosfwd>

void write_block_lammps_dump( const block_data &b, std::ostream &o );


#endif // BLOCK_DATA_WRITERS_H

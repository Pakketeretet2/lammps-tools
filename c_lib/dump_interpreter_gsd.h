#ifndef DUMP_INTERPRETER_GSD_H
#define DUMP_INTERPRETER_GSD_H

#include "block_data.h"
#include "dump_reader.h"
#include "reader_core_binary.h"

#include <string>
#include <iosfwd>

#include "gsd/gsd.h"


// A dump interpreter for the HOOMD-blue GSD schema.
class dump_interpreter_gsd : public dump_interpreter
{
public:
	dump_interpreter_gsd( const std::string &fname );
	virtual ~dump_interpreter_gsd();
	virtual bool next_block( reader_core *r, block_data &b );

	int get_chunk_data( const std::string &name, void *dest );
	
private:
	int status;
	gsd_handle *gh;

	uint64_t current_frame;
};



#endif // DUMP_INTERPRETER_GSD_H

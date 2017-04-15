#ifndef DUMP_INTERPRETER_GSD_H
#define DUMP_INTERPRETER_GSD_H

#include "block_data.h"
#include "dump_reader.h"

#include <string>
#include <iosfwd>

#include "gsd/gsd.h"

template <typename data_type>
union gsd_utf8_char_templ
{
	data_type s;
	char c[sizeof(data_type)];
};

typedef gsd_utf8_char_templ<char> gsd_utf8_char;


// A dump interpreter for the HOOMD-blue GSD schema.
class dump_interpreter_gsd : public dump_interpreter
{
public:
	enum { TYPE_BUFFER_SIZE = 64 };
	
	dump_interpreter_gsd( const std::string &fname );
	virtual ~dump_interpreter_gsd();
	virtual int next_block( block_data &b );

	int get_chunk_data( const std::string &name, void *dest );
	
private:
	int status;
	gsd_handle *gh;

	uint64_t current_frame;
};



#endif // DUMP_INTERPRETER_GSD_H

#ifndef DUMP_INTERPRETER_H
#define DUMP_INTERPRETER_H

#include <iosfwd>
#include <string>

#include "block_data.h"


class dump_interpreter
{
public:
	dump_interpreter( const std::string &dname );
	dump_interpreter( std::istream &in );
	virtual ~dump_interpreter(){}

	virtual int next_block( block_data &block );
	virtual int last_block( block_data &block );

	///< Calling next_block_meta only leaves the block in a half-
	///< half-initialised state! You can read out the meta but not the rest.
	virtual int next_block_meta( block_data &block );
	virtual int next_block_body( block_data &block );

	virtual bool eof()  const = 0;
	virtual bool good() const = 0;
	
	
	virtual const std::string &dname() const { return dname_; }
	
private:
	std::string dname_;
};


#endif // DUMP_INTERPRETER_H

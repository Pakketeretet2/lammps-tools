#ifndef READER_CORE_PLAIN_H
#define READER_CORE_PLAIN_H

#include "block_data.h"
#include "dump_reader.h"

#include <iosfwd>

class reader_core_plain : public reader_core
{
public:
	reader_core_plain( const std::string &fname );
	reader_core_plain( std::istream &istream );
	
	virtual ~reader_core_plain();

	virtual bool getline( std::string &line );
	virtual void rewind();

private:
	std::istream *in;
	bool got_file;

	int lc;
};





#endif // READER_CORE_PLAIN_H

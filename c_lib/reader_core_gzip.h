#ifndef READER_CORE_GZIP_H
#define READER_CORE_GZIP_H

#include "block_data.h"
#include "dump_reader.h"
#include <fstream>

#ifdef HAVE_BOOST_GZIP
#  include <boost/iostreams/filtering_stream.hpp>
#endif



class reader_core_gzip : public reader_core
{
public:
	reader_core_gzip( const std::string &fname );
	virtual ~reader_core_gzip();

	virtual bool getline( std::string &line );
	virtual void rewind();

	virtual void set_debug(bool);
	virtual operator bool() const
	{
		return static_cast<bool>(infile);
	}

	virtual int peek()
	{
		return in.peek();
	}
	
	
private:
	std::ifstream infile;
#ifdef HAVE_BOOST_GZIP
	boost::iostreams::filtering_istream in;
#else
	std::istream &in;
#endif
	bool debug;
};


#endif // READER_CORE_GZIP_H

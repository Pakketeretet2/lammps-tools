#ifndef TEXT_READERS_H
#define TEXT_READERS_H

#include <iosfwd>
#include <string>

#ifdef HAVE_BOOST_GZIP
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#endif

/** 
    A general interface to read in text files.
*/
class text_reader
{
public:
	text_reader( const std::string &fname ){}
	virtual ~text_reader(){}
	virtual bool getline( std::string &line ){ return 0; }
	virtual int  peek(){ return 0; }

	virtual bool eof()  const = 0;
	virtual bool good() const = 0;
};

/** 
    Reads in plain text files:
*/
class text_reader_plain : public text_reader
{
public:
	text_reader_plain( const std::string &fname );
	virtual ~text_reader_plain();
	virtual bool getline( std::string &line );
	virtual int peek();

	virtual bool eof()  const;
	virtual bool good() const;
	
private:
	std::ifstream *in;
};


/** 
    Reads in gzipped text files:
*/
class text_reader_gzip : public text_reader
{
public:
	text_reader_gzip( const std::string &fname );
	virtual ~text_reader_gzip();

	virtual bool getline( std::string &line );
	virtual int peek();

	virtual bool eof()  const;
	virtual bool good() const;

private:
	std::ifstream *in_file;
#ifdef HAVE_BOOST_GZIP
	boost::iostreams::filtering_istream in;
#endif // HAVE_BOOST_GZIP
};





#endif // TEXT_READERS_H

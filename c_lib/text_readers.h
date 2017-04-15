#ifndef TEXT_READERS_H
#define TEXT_READERS_H

#include <iosfwd>
#include <string>

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

#endif // TEXT_READERS_H

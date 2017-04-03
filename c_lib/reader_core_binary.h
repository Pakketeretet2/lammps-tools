#ifndef READER_CORE_BINARY
#define READER_CORE_BINARY

#include "block_data.h"
#include "dump_reader.h"

#include <iosfwd>
#include <typeinfo>

template <typename T>
struct read_buffer {

	read_buffer( int nbytes ) : buffer(nullptr)
	{
		raw_buffer = new char[nbytes];
		buffer = reinterpret_cast<T*>( raw_buffer );
	};

	~read_buffer()
	{
		if( buffer ) delete [] buffer;
	}

	T    *buffer;
	char *raw_buffer;
};

class reader_core_bin : public reader_core
{
public:

	
	
	reader_core_bin( const std::string &fname );
	reader_core_bin( std::istream &istream );
	
	virtual ~reader_core_bin();

	template <typename T>
	void get_data( T *dest, uint n )
	{
		read_buffer<T> buffer( n*sizeof(T) );
		in->read( buffer.raw_buffer, n*sizeof(T) );
		current_byte += n*sizeof(T);
		for( int i = 0; i < n; ++i ){
			dest[i] = buffer.buffer[i];
		}
	}

	template <typename T>
	T get_data()
	{
		read_buffer<T> buffer( sizeof(T) );
		in->read( buffer.raw_buffer, sizeof(T) );
		current_byte += sizeof(T);
		T t;
		t = buffer.buffer[0];
		return t;
	}
	
	void get_char_data( char *dest, uint n );
	unsigned long int get_current_byte();
	
private:
	std::istream *in;
	bool got_external_stream;
	unsigned long int current_byte;
	
};





#endif // READER_CORE_BINARY

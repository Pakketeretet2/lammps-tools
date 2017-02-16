#ifndef MY_OUTPUT_HPP
#define MY_OUTPUT_HPP

#include <iostream>


class my_ostream {
public:
	my_ostream( std::ostream &o ) : o_(o){}

	template <typename T>
	my_ostream &operator<<( const T &t )
	{
#ifdef VERBOSE_LIB
		o_ << t;
#endif
		return *this;
	}

	
private:	
	std::ostream &o_;
};



#endif // MY_OUTPUT_HPP

#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <algorithm>

#include "list.h"

inline bool starts_with( const std::string &s, const std::string &begin )
{
	return s.compare(0, begin.length(), begin ) == 0;
}

inline bool ends_with( const std::string &s, const std::string &end )
{
	if( end.length() > s.length() ){
		return false;
	}else{
		return s.compare(s.length() - end.length() - 1, s.length()-1, end ) == 0;
	}
}

inline std::string rstrip( std::string s, char delim )
{
	std::size_t idx = s.find( delim );
	if( idx == s.length() ) return s;

	s = s.substr( 0, idx );
	std::size_t top = s.find_last_not_of( " \t" ) + 1;
	return s.substr( 0, top );
}


inline int word_count( const std::string &s )
{
	std::stringstream ss(s);
	std::string w;
	int count = 0;
	while( ss >> w ) ++count;
	return count;
}


template <typename list_type, typename val_type>
inline bool list_has( const list_type &l, val_type val )
{
	return std::find( l.begin(), l.end(), val ) != l.end();
}


template <typename T1, typename T2> inline
T1 max( T1 a, T2 b )
{
	return a > b ? a : b;
}

template <typename T1, typename T2> inline
T1 min( T1 a, T2 b )
{
	return a > b ? b : a;
}


template <typename container_type> inline
double mean( const container_type &c )
{
	double s = std::accumulate( c.begin(), c.end(), 0.0 );
	return s / static_cast<double>( c.size() );
}

template <typename container_type> inline
double var( const container_type &c )
{
	double m = mean(c);
	double s = std::inner_product( c.begin(), c.end(), c.begin(), 0.0 );
	return s / static_cast<double>( c.size() -1 ) - m*m;
}



#endif // UTIL_H

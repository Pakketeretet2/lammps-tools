#include "id_map.h"

id_map::id_map( const arr1i &ids )
{
	for( uint i = 0; i < ids.size(); ++i ){
		m[ids[i]] = i;
	}
}

py_int id_map::operator[]( py_int id ) const
{
	std::map<py_int,py_int>::const_iterator i = m.find(id);
	if( i == m.end() ){
		return -1;
	}else{
		return i->second;
	}
}


py_int id_map::id_to_index( py_int id ) const
{
	return (*this)[id];
}

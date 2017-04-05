#ifndef ID_MAP_H
#define ID_MAP_H

/*!
  \file id_map.h
  @brief Some tools for indexing based on ids.
  
  \ingroup cpp_lib
*/

#include "types.h"
#include <map>
#include <vector>


/*! @brief This class contains some tools for mapping ids onto array indices.

  \ingroup cpp_lib
*/
class id_map {
public:
	/// Empty constructor
	id_map( ) {}
	/// Constructor based on given atom ids
	id_map( const arr1i &ids );
	id_map( const py_int *ids, py_int N );
	
	template <typename int_type>
	id_map( const std::vector<int_type> &ids )
	{
		for( uint i = 0; i < ids.size(); ++i ){
			m[ids[i]] = i;
		}

	}

	/// Returns the index in ids (and other arrays) of given id
	py_int operator[]( py_int id ) const;

	/// Larger, "named" function names for if you are not sure:
	py_int id_to_index( py_int id ) const;
	
private:
	std::map<py_int, py_int> m; ///< std::map the mapping is stored in
};


#endif /* ID_MAP_H */

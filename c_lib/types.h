#ifndef TYPES_H
#define TYPES_H

#include "c_lib.h"

#include <cstdint>
#include <cstddef>
#include <iostream>

/*!
  @file  types.h
  @brief This file contains definitions and other important things.
  
  @detailed
  This header contains some important stuff:
   - Definitions of the types this lib expects Python to pass
   - An interface for dealing with nd.arrays passed from Python as void-ptrs
   - Some useful math constants and typedefs.

   \ingroup cpp_lib
*/


typedef int64_t       bigint;    ///< signed long of guaranteed size (64 bits)
typedef int64_t       py_int;    ///< int as passed from Python
typedef uint64_t      py_uint;   ///< uint as passed from Python
typedef uint64_t      biguint;   ///< unsigned long of guaranteed size (64 bits)
typedef double        py_float;  ///< float as passed from python
typedef unsigned int  uint;      ///< "small" unsigned int
typedef wchar_t       py_char;   ///< widechar type to deal with Python strings


/*!
  @brief To neatly interface with numpy-arrays from Python, use the following
  type trait struct. It defines an implementation for the typical
  access-operators based on the dimensionality of the array.

  @detailed The unspecialized variant handles Nxdim-dimensional arrays,
  making operator[] act like for double arrays (i.e., the first operator[](N)
  returns p + dim*N so that access via [M][N] corresponds to what you would
  expect (namely p + dim*M + N).

  \ingroup cpp_lib
*/
template <typename T, py_int dim, bool read_only = false>
struct access_traits
{
	/// Constructor from void-ptr to data
	access_traits(void *p)       : data( static_cast<T*>(p) ){}

	/// Destructor should be empty, Python handles memory management
	~access_traits(){}

	/// Access operator[] like for a double array, read-only:
	const T* operator[]( uint i ) const
	{ return data + dim*i; }

	/// Access operator[] like for double array, read-write:
	T* operator[]( uint i )
	{ return data + dim*i; }

	/// Access the raw data as void-ptr:
	void *raw_data()
	{ return data; }

	/// Access the raw data as typed ptr:
	T *t_data()
	{ return data; }

	/// Get raw ptr to data, read-only
	const void *raw_data() const
	{ return data; }

	/// Get typed ptr to data, read-only
	const T *t_data() const
	{ return data; }

	/// \private Changes what data points to. Only for experimental stuff!
	void set_data( T* data ) { this->data = data; }
	
private:
	/// Typed ptr to underlying Python data
	T *data;
};



/*! \relates access_traits<T, dim, read_only>
  @detailed Read-only access to pointer data.
*/
template <typename T, py_int dim>
struct access_traits<T, dim, true>
{
	/// Constructor from void-ptr to data
	access_traits(const void *p) : data( static_cast<T*>(p) ){}


	/// Access operator[] like for a double array, read-only:
	const T* operator[]( uint i ) const
	{ return data + dim*i; }

	/// Access the raw data as void-ptr:
	const void *raw_data()
	{ return data; }

	/// Access the raw data as typed ptr:
	const T *t_data()
	{ return data; }

	/// Get raw ptr to data, read-only
	const void *raw_data() const
	{ return data; }

	/// Get typed ptr to data, read-only
	const T *t_data() const
	{ return data; }

	/// \private Changes what data points to. Only for experimental stuff!
	void set_data( const T* data ) { this->data = data; }
	
private:
	/// Typed ptr to underlying Python data
	const T *data;
};



/*! \relates acces_traits<T, dim>

  @detailed Specialized variant of access_traits<T,dim> that handles
  one-dimensional arrays, making operator[] directly return the value at
  the given position.

  \relates access_traits

  \ingroup cpp_lib
*/
template<typename T>
struct access_traits<T, 1>
{
	/// Constructor taking a void-ptr to data
	access_traits(void *p)       : data( static_cast<T*>(p) ){}
	/// Empty destructor, Python does memory management.
	~access_traits(){}

	/// Return value at i, read-only (or rather, copy)
	T operator[]( uint i ) const
	{ return data[i]; }

	/// Return ref to value at i, read/write
	T& operator[]( uint i )
	{ return data[i]; }

	/// Get raw ptr to data, read-only
	void *raw_data()
	{ return data; }

	/// Get typed ptr to data, read-only
	T *t_data()
	{ return data; }

	/// Get raw ptr to data, read-only
	const void *raw_data() const
	{ return data; }

	/// Get typed ptr to data, read-only
	const T *t_data() const
	{ return data; }
	
	/// \private Change address of data, for experimental stuff only!
	void set_data( T* data ) { this->data = data; }

private:
	/// Typed ptr to underlying Python data.
	T *data;
};







/*!
  @brief Class that wraps around a numpy ndarray passed from Python.
  \relates access_traits

  \ingroup cpp_lib
*/
template <typename T, py_int stride, bool read_only = false>
class nd_array : public access_traits<T, stride, read_only>
{
public:
	nd_array( void *data, uint size )
		: access_traits<T, stride,read_only>(data), s(size) {}
	
	~nd_array(){}
	uint size() const { return s; }


	/* Standard C++-like iterators: */
	T *begin()
	{
		return access_traits<T,stride,read_only>::t_data();
	}

	T *end()
	{
		return access_traits<T,stride,read_only>::t_data() + size();
	}


	const T *begin() const
	{
		return access_traits<T,stride,read_only>::t_data();
	}

	const T *end() const
	{
		return access_traits<T,stride,read_only>::t_data() + size();
	}
	
private:
	const uint s;
};



typedef nd_array<py_float,3> arr3f; ///< 3D floating point array
typedef nd_array<py_float,1> arr1f; ///< 1D floating point array
typedef nd_array<py_int,3>   arr3i; ///< 3D signed integer array
typedef nd_array<py_int,1>   arr1i; ///< 1D signed integer array

typedef nd_array<py_float,4> arr4f; ///< 4D floating point array
typedef nd_array<py_int,4>   arr4i; ///< 4D signed integer array


// Same but with read-only access:
typedef nd_array<py_float,3,true> arr3fc; ///< 3D floating point array
typedef nd_array<py_float,1,true> arr1fc; ///< 1D floating point array
typedef nd_array<py_int,3,true>   arr3ic; ///< 3D signed integer array
typedef nd_array<py_int,1,true>   arr1ic; ///< 1D signed integer array

typedef nd_array<py_float,4,true> arr4fc; ///< 4D floating point array
typedef nd_array<py_int,4,true>   arr4ic; ///< 4D signed integer array


namespace math_const {
/// The mathematical constant pi
constexpr const double pi = 3.1415926535897932384626433832795028841971693993751;

}

#endif /* TYPES_H */

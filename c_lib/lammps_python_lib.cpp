#include "c_lib.h"
#include "types.h"

/*!
  \ingroup  cpp_lib
*/



py_int test( py_int N, void *arr_ptr, void *double_ptr,
             void *double_pptr, py_int stride )
{
	// Cast arr to the correct type, int:
	arr1i arr(arr_ptr,N);
	arr1f arrd(double_ptr,N);
	arr3f double_arr(double_pptr,N);
	
	
	fprintf(stdout,"Testing...\n");
	for( int i = 0; i < N; ++i ){
		fprintf(stdout, "%d: %ld\n",i, arr[i] );
	}

	for( int i = 0; i < N; ++i ){
		fprintf(stdout, "%d: %f\n",i, arrd[i] );
	}

	for( int i = 0; i < N; ++i ){
		for( int j = 0; j < 3; ++j ){
			fprintf(stdout, "%d, %d: %f  ",i,j, double_arr[i][j]);
		}
		fprintf(stdout,"\n");
	}

	return N;
}


void interface_test_int( py_int *data, int N )
{
	for( int i = 0; i < N; ++i ){
		fprintf( stdout, "%d   %ld\n", i, data[i] );
	}
	
}

void interface_test_double( py_float *data, int N )
{
	fprintf(stdout, "ptr lives at %p\n", data);
	for( int i = 0; i < N; ++i ){
		fprintf( stdout, "%d   %f\n", i, data[i] );
	}
}


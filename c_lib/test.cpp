#include "types.h"
#include "test.h"

#include <iostream>
#include "util.h"


extern "C" {
void test_types( void *px, py_int N, void *pids, void *ptypes,
                 py_int periodic, py_float *xlo )
{
	std::cout << "Got " << N << " atoms.\n";
	std::cout << "xlo: [ " << xlo[0] << ", " << xlo[1] << ", "
	          << xlo[2] << "  ], periodic = " << periodic << "\n";
	std::cout << "First atom:\n";


	arr1i ids(pids,N);
	arr1i types(ptypes,N);
	double *xx = static_cast<double*>(px);
	
	
	std::cout << ids[0] << " " << types[0]
	          << " " << xx[0] << " " << xx[1] << " " << xx[2] << "\n";

	std::cout << "Using arr3f for x:\n";
	std::cout.flush();
	arr3f x(px,N);
	std::cout << ids[0] << " " << types[0]
	          << " " << x[0][0] << " " << x[0][1] << " " << x[0][2] << "\n";
}

void test_modification( void *px, py_int N, py_int stride)
{
	std::cout << "Got pointer to " << N << " py_float[3]s.\n";
	arr3f x(px,N);
	
	std::cout << "Values:\n";
	for( py_int i = 0; i < N; ++i ){
		std::cout << i << ":";
		for( py_int j = 0; j < stride; ++j ){
			std::cout << " " << x[i][j];
		}
		std::cout << "\n";
	}

	py_int i = min( 3, N-1 );
	py_int j = min( 25, stride-1 );
	std::cout << "modifying value at " << i << ", " << j
	          << " by adding 2.0:\n";
	x[i][j] += 2;

	for( py_int i = 0; i < N; ++i ){
		std::cout << i << ":";
		for( py_int j = 0; j < stride; ++j ){
			std::cout << " " << x[i][j];
		}
		std::cout << "\n";
	}

	
}

} // extern "C"

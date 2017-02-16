#include "domain.h"


#include <cmath>
#include <iostream>

void distance_wrap( py_float *r, const py_float *x1, const py_float *x2,
                    const py_float *xlo, const py_float *xhi, py_int periodic )
{
	/*
	  Here is a written out example for checking the code:
	  x1 = (20.4535, 25.0013, 13.9697 )
	  x2 = ( 22.9008, 13.5519, 15.2423 )
	  L = ( 17.4572, 17.4572, 13.9658 ).
	  dx = ( 4.6694   14.7437   -2.1532 )
	*/
	
	py_float dx = x2[0] - x1[0];
	py_float dy = x2[1] - x1[1];
	py_float dz = x2[2] - x1[2];


	if( periodic & PERIODIC_X ){
		py_float Lx = xhi[0] - xlo[0];
		if( dx >= 0.5*Lx ){ // 4.6697 < 0.5*17.4..
			dx = -Lx + dx;
		}else if( dx < -0.5*Lx ){ // 4.6697 < -...
			dx = Lx + dx;
		}
	}

	if( periodic & PERIODIC_Y ){
		py_float Ly = xhi[1] - xlo[1];
		if( dy >= 0.5*Ly ){ // 14.7437 > 0.5*17.4572
			dy = -Ly + dy; // dy = - 17.4572 +  14.7437 = -2.7133
		}else if( dy < -0.5*Ly ){ // -2.7133 > -0.5*17.45...
			dy = Ly + dy;
		}
	}
	if( periodic & PERIODIC_Z ){
		py_float Lz = xhi[2] - xlo[2];
		if( dz >= 0.5*Lz ){ // -2.1532 < 0.5-13....
			dz = -Lz + dz;
		} else if( dz < -0.5*Lz ){ // -2.1532 > -0.5*13.9658...
			dz = Lz + dz;
		}
	}
	r[0] = dx;
	r[1] = dy;
	r[2] = dz;
}

void distance( py_float *r, const py_float *x1, const py_float *x2 )
{
	distance_wrap( r, x1, x2, nullptr, nullptr, PERIODIC_NONE );
}

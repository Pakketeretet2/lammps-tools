#ifndef TRIANGULATE_H
#define TRIANGULATE_H

#include "types.h"

#include <vector>
#include <cmath>



inline void copy_arr( double *dest, const double *src, int N )
{
	for( int i = 0; i < N; ++i ){
		dest[i] = src[i];
	}
}



struct triangle
{
	py_int i1, i2, i3;
	double x1[3], x2[3], x3[3];

	triangle( py_int j1, py_int j2, py_int j3,
	          const double y1[3], const double y2[3], const double y3[3] )
	{
		// Find the right mapping, which is indices sorted
		// from low to high.
		if ( j1 < j2 ){
			if( j1 < j3 ){
				// Order is j1 < (j2,j3)
				i1 = j1;
				copy_arr( x1, y1, 3 );
				
				if( j2 < j3 ){
					i2 = j2;
					i3 = j3;
					copy_arr( x2, y2, 3 );
					copy_arr( x3, y3, 3 );
				}else{
					i2 = j3;
					i3 = j2;
					copy_arr( x2, y3, 3 );
					copy_arr( x3, y2, 3 );
				}
			}else{
				// Order is j3 < (j1,j2):
				i1 = j3;
				copy_arr( x1, y3, 3 );
				if( j1 < j2 ){
					i2 = j1;
					i3 = j2;
					copy_arr( x2, y1, 3 );
					copy_arr( x3, y2, 3 );
				}else{
					i2 = j2;
					i3 = j1;
					copy_arr( x2, y2, 3 );
					copy_arr( x3, y1, 3 );
				}
			}
		}else{
			if( j2 < j3 ){
				// Order is j2 < (j1,j3)
				i1 = j2;
				copy_arr( x1, y2, 3 );
				if( j1 < j3 ){
					i2 = j1;
					i3 = j3;
					copy_arr( x2, y1, 3 );
					copy_arr( x3, y3, 3 );
				}else{
					i2 = j3;
					i3 = j1;
					copy_arr( x2, y3, 3 );
					copy_arr( x3, y1, 3 );
				}
			}else{
				// Order is j3 < (j1,j2)
				i1 = j3;
				copy_arr( x1, y3, 3 );
				if( j1 < j2 ){
					i2 = j1;
					i3 = j2;
					copy_arr( x2, y1, 3 );
					copy_arr( x3, y2, 3 );
				}else{
					i2 = j2;
					i3 = j1;
					copy_arr( x2, y2, 3 );
					copy_arr( x3, y1, 3 );
				}
			}
		}
		// Some safety checks:
		if( i1 > i2 ){
			std::cerr << "WTF?! i1 > i2 ( " << i1
			          << " > " << i2 << " )!\n";
		}else if( i1 > i3 ){
			std::cerr << "WTF?! i1 > i3 ( " << i1
			          << " > " << i3 << " )!\n";			
		}
	}
	

	bool operator==( const triangle &o ) const
	{
		return (i1==o.i1) && (i2==o.i2) && (i3==o.i3);
	}


	double area() const
	{
		// area using Heron's formula:
		double dx1[3];
		double dx2[3];
		double dx3[3];
		dx1[0] = x1[0] - x2[0];
		dx1[1] = x1[1] - x2[1];
		dx1[2] = x1[2] - x2[2];

		dx2[0] = x1[0] - x3[0];
		dx2[1] = x1[1] - x3[1];
		dx2[2] = x1[2] - x3[2];

		dx3[0] = x3[0] - x2[0];
		dx3[1] = x3[1] - x2[1];
		dx3[2] = x3[2] - x2[2];

		double a = std::sqrt(dx1[0]*dx1[0] + dx1[1]*dx1[1] + dx1[2]*dx1[2]);
		double b = std::sqrt(dx2[0]*dx2[0] + dx2[1]*dx2[1] + dx2[2]*dx2[2]);
		double c = std::sqrt(dx3[0]*dx3[0] + dx3[1]*dx3[1] + dx3[2]*dx3[2]);

		double s = 0.5*( a + b + c );
		double sa = s - a;
		double sb = s - b;
		double sc = s - c;
		
		return std::sqrt( s*sa*sb*sc );
	}
};



extern "C"
{

void triangulate( void *x, py_int N, py_int *ids, py_int *types,
                  py_float rc, py_int periodic,
                  const py_float *xlo, const py_float *xhi, py_int dims,
                  py_int method, const char *pname );
}


void triangulate_impl( const arr3f &x, py_int N, py_int *ids, py_int *types,
                       py_float rc, py_int periodic, const py_float *xlo,
                       const py_float *xhi, py_int dims, py_int method,
                       const char *pname );


void triangulate_block( class block_data &b, py_float rc, py_int periodic,
                        py_int dims, py_int method, std::vector<triangle> &triangles );

double triangulation_area( class block_data &b, std::vector<triangle> &triangles );

#endif // TRIANGULATE_H

#include "similarity.h"
#include "types.h"
#include "random_pcg.h"

#include <iostream>
#include <cmath>
#include <armadillo>

#include <memory>

constexpr const static double HALF_PI = 1.57079632679490;
constexpr const static double PI =      3.14159265358979;
extern "C" {


double map_clusters( void *xi, void *xj,  py_int N,
                     py_int *ids, py_int *jds, py_int periodic,
                     void *xlo, void *xhi, py_int dims,
                     void *rotated, py_float *rot_axis )
{
	arr3f xxi(xi,N), xxj(xj,N);
	arr3f rot(rotated,N);

	test_rotations();

	py_float *xxlo = static_cast<py_float*>(xlo);
	py_float *xxhi = static_cast<py_float*>(xhi);

	return 0;	
	/*
	return map_clusters_impl( xxi, xxj, N, ids, jds,
	                          periodic, xxlo, xxhi, dims,
	                          rot, rot_axis );
	*/    
}


py_float rotate_to_template( void *xxi, void *xxt,  py_int N,
                             py_int periodic, void *xxlo, void *xxhi,
                             py_int dims, void *mmaj_ax, void *mmin_ax,
                             void *mmaj_ax_t, void *mmin_ax_t,
                             void *x_rotated )
{
	arr3f xi( xxi, N );
	arr3f xt( xxt, N );

	py_float *xlo      = static_cast<py_float*>( xxlo );
	py_float *xhi      = static_cast<py_float*>( xxhi );
	py_float *maj_ax_t = static_cast<py_float*>( mmaj_ax_t );
	py_float *min_ax_t = static_cast<py_float*>( mmin_ax_t );
	py_float *maj_ax   = static_cast<py_float*>( mmaj_ax );
	py_float *min_ax   = static_cast<py_float*>( mmin_ax );
	py_float *xrot     = nullptr;
	if( x_rotated ){
		// std::cerr << "Storing rotated values in array at " << x_rotated << ".\n";
		xrot = static_cast<py_float*>( x_rotated );
	}else{
		// std::cerr << "Not storing rotated values...\n";
	}
	
	return rotate_to_template_impl( xi, xt, N, periodic, xlo, xhi,
	                                dims, maj_ax, min_ax, maj_ax_t,
	                                min_ax_t, xrot );

}



	
} // extern "C"


template <typename T> inline
T sign( T& a )
{
	return (a > 0) - (a < 0);
}


void quat_prod( const double *q1, const double *q2, double *result )
{
	double a1 = q1[0];
	double b1 = q1[1];
	double c1 = q1[2];
	double d1 = q1[3];

	double a2 = q2[0];
	double b2 = q2[1];
	double c2 = q2[2];
	double d2 = q2[3];

	result[0] = a1*a2 - b1*b2 - c1*c2 - d1*d2;
	result[1] = a1*b2 + b1*a2 + c1*d2 - d1*c2;
	result[2] = a1*c2 - b1*d2 + c1*a2 + d1*b2;
	result[3] = a1*d2 + b1*c2 - c1*b2 + d1*a2;
}


double dot( const py_float *x, const py_float *y )
{
	return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

void cross( const py_float *x, const py_float *y, py_float *res )
{
	res[0] = x[1]*y[2] - x[2]*y[1];
	res[1] = x[2]*y[0] - x[0]*y[2];
	res[2] = x[0]*y[1] - x[1]*y[0];
}


double norm( const py_float *x, int N )
{
	double sum = 0;
	for( int i = 0; i < N; ++i ){
		sum += x[i]*x[i];
	}
	return sqrt(sum);
}

double norm3_2( const py_float *x )
{
	return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

double norm3( const py_float *x )
{
	return sqrt( norm3_2(x) );
}

double dist_2( const py_float *x, const py_float *y, py_int N )
{
	double r2 = 0.0;
	for( py_int i = 0; i < N; ++i ){
		double dx = x[i] - y[i];
		r2 += dx*dx;
	}
	return r2;
}



double dist3_2( const py_float *x, const py_float *y )
{
	py_float z[3];
	z[0] = x[0] - y[0];
	z[1] = x[1] - y[1];
	z[2] = x[2] - y[2];
	return norm3_2(z);
}

double dist3( const  py_float *x, const py_float *y )
{
	return sqrt( dist3_2(x,y) );
}


bool close_enough( const double *a, const double *b, int N, double tol = 1e-8 )
{
	double dr2 = dist3_2( a, b );

	double dx = a[0] - b[0];
	double dy = a[1] - b[1];
	double dz = a[2] - b[2];
	return dr2 < tol;
}



void rotate_vector( double *rotate, const double *rot_axis, double df )
{
	// Quaternion rotation:
	double q[4], p[4], qt[4];
	p[0] = 0;
	p[1] = rotate[0];
	p[2] = rotate[1];
	p[3] = rotate[2];

	double c = cos(0.5*df);
	double s = sin(0.5*df);
	q[0] = c;
	q[1] = rot_axis[0]*s;
	q[2] = rot_axis[1]*s;
	q[3] = rot_axis[2]*s;
	
	qt[0] =  q[0];
	qt[1] = -q[1];
	qt[2] = -q[2];
	qt[3] = -q[3];

	double result[4];
	quat_prod( q, p, result );
	quat_prod( result, qt, result );
	
	rotate[0] = result[1];
	rotate[1] = result[2];
	rotate[2] = result[3];
}

void gramm_schmid( const py_float *x1, const py_float *x2, py_float *res )
{
	double ab = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
	double b2 = norm3_2(x2);

	ab /= b2;
	
	res[0] = x1[0] - ab * x2[0]; // / b2;
	res[1] = x1[1] - ab * x2[1];
	res[2] = x1[2] - ab * x2[2];
}


// Tests the rotation and quaternion product functions.
void test_rotations()
{
	// Test quaternions:
	py_float q1[4], q2[4];
	q1[0] = 0;
	q1[1] = 1;
	q1[2] = 2;
	q1[3] = 0;

	q2[0] = 1;
	q2[1] = 1;
	q2[2] = q2[3] = 0;

	py_float q3[4];
	quat_prod( q1, q2, q3 );
	std::cerr << "q1 = ( " << q1[0] << ", " << q1[1] << ", " << q1[2] << ", " << q1[3] << " ).\n";
	std::cerr << "q2 = ( " << q2[0] << ", " << q2[1] << ", " << q2[2] << ", " << q2[3] << " ).\n";
	std::cerr << "q3 = ( " << q3[0] << ", " << q3[1] << ", " << q3[2] << ", " << q3[3] << " ).\n";
	// q3 should be
	py_float q3act[4];
	q3act[0] = -1;
	q3act[1] = 1;
	q3act[2] = 2;
	q3act[3] = -2;
	
	py_float x[3];
	py_float y[3];
	py_float z[3];
	
	x[0] = 1; x[1] = x[2] = 0;
	y[0] = y[2] = 0; y[1] = 1;
	z[0] = z[1] = 0; z[2] = 1;
	
	std::cerr << "\n\n";
	
	std::cerr << "x is now ( " << x[0] << ", " << x[1] << ", " << x[2] << " ).\n";
	rotate_vector( x, z, HALF_PI );
	std::cerr << "x is now ( " << x[0] << ", " << x[1] << ", " << x[2] << " ).\n";

	

	if( !close_enough( x, y, 3, 1e-10 ) ){
		std::cerr << "Error in rotations! Aborting.\n";
		std::terminate();
	}
}


void test_cross()
{
	py_float x[3], y[3], z[3];
	x[1] = x[2] = y[0] = y[2] = z[0] = z[1] = 0;
	x[0] = y[1] = z[2] = 1;
	py_float r[3];
	
	cross(x,y,r);
	std::cerr << "x cross y = ( " << r[0] << ", " << r[1] << ", " << r[2] << " ).\n";
	std::cerr << "should be   ( " << z[0] << ", " << z[1] << ", " << z[2] << " ).\n";

	cross(x,z,r);
	std::cerr << "x cross z = ( " << r[0] << ", " <<  r[1] << ", " << r[2] << " ).\n";
	std::cerr << "should be   ( " << y[0] << ", " << -y[1] << ", " << y[2] << " ).\n";

	cross(y,z,r);
	std::cerr << "y cross z = ( " << r[0] << ", " <<  r[1] << ", " << r[2] << " ).\n";
	std::cerr << "should be   ( " << x[0] << ", " <<  x[1] << ", " << x[2] << " ).\n";
	std::cerr << "\n\n\n";
		

}

double rotate_onto( const py_float *v, const py_float *target,
                    py_float *axis )
{
	// Finds axis and angle that rotate v __onto__ target.

	cross( v, target, axis );
	double n = 1.0 / norm3(axis);
	axis[0] *= n;
	axis[1] *= n;
	axis[2] *= n;

	double d = dot( v, target );
	d /= (norm3(v)*norm3(target));
	double angle = acos(d);
	
	return angle;
}

void copy_ensemble( const arr3f &x, arr3f &cp, py_int N )
{
	for( py_int i = 0; i < N; ++i ){
		cp[i][0] = x[i][0];
		cp[i][1] = x[i][1];
		cp[i][2] = x[i][2];
		
	}
}


void rotate_ensemble( arr3f &rot_data, const py_float *rot_ax,
                      py_int N, double angle )
{
	for( py_int i = 0; i < N; ++i ){
		rotate_vector( rot_data[i], rot_ax, angle );
	}
}


py_float rotate_to_template_impl( const arr3f &xi, const arr3f &xt, py_int N,
                                  py_int periodic, py_float *xlo, py_float *xhi,
                                  py_int dims, py_float *maj_ax, py_float *min_ax,
                                  py_float *maj_ax_t, py_float *min_ax_t,
                                  py_float *xrot )
{
	double *rot_data_;
	
	if( xrot ) rot_data_ = xrot;
	else       rot_data_ = new double[3*N];
	
	arr3f rot_data(rot_data_, N);

	double ang1, ang2;
	double rot_ax[3];

	// test_cross();


	// std::cerr << "Rotating " << N << " particles to template.\n";
	// std::cerr << "Source has major and minor axes:\n";
	// std::cerr << "     ( " << maj_ax[0] << ", " <<  maj_ax[1] << ", " << maj_ax[2] << " ).\n";
	// std::cerr << "     ( " << min_ax[0] << ", " <<  min_ax[1] << ", " << min_ax[2] << " ).\n";
	// std::cerr << "Template has major and minor axes:\n";
	// std::cerr << "     ( " << maj_ax_t[0] << ", " <<  maj_ax_t[1] << ", " << maj_ax_t[2] << " ).\n";
	// std::cerr << "     ( " << min_ax_t[0] << ", " <<  min_ax_t[1] << ", " << min_ax_t[2] << " ).\n";

	copy_ensemble( xi, rot_data, N );

	
	// Check if major axes are too much alike:


	double major_dot = dot(maj_ax_t, maj_ax);
	if( major_dot > 1 - 1e-4 ){
		// std::cerr << "Axes almost aligned...\n";
		
	}else if( major_dot < -1 + 1e-4 ){
		// std::cerr << "Axes almost anti-aligned...\n";

		// Idea: Just flip everything upside-down about the minor axis

		rotate_ensemble( rot_data, min_ax, N, PI );
		rotate_vector( maj_ax, min_ax, PI );
	}else{
		// std::cerr << "Axes not aligned...\n";
		

		ang1 = rotate_onto( maj_ax, maj_ax_t, rot_ax );
		rotate_ensemble( rot_data, rot_ax, N, ang1 );

		// Rotate axes as well:
		rotate_vector( maj_ax, rot_ax, ang1 );
		rotate_vector( min_ax, rot_ax, ang1 );
	}

	
	// std::cerr << "After rotation about major axes, this is the state:\n";
	// std::cerr << "Source has major and minor axes:\n";
	// std::cerr << "     ( " << maj_ax[0] << ", " <<  maj_ax[1] << ", " << maj_ax[2] << " ).\n";
	// std::cerr << "     ( " << min_ax[0] << ", " <<  min_ax[1] << ", " << min_ax[2] << " ).\n";
	// std::cerr << "Template has major and minor axes:\n";
	// std::cerr << "     ( " << maj_ax_t[0] << ", " <<  maj_ax_t[1] << ", " << maj_ax_t[2] << " ).\n";
	// std::cerr << "     ( " << min_ax_t[0] << ", " <<  min_ax_t[1] << ", " << min_ax_t[2] << " ).\n";

	double minor_dot = dot(min_ax_t, min_ax);
	if( minor_dot > 1 - 1e-4 ){
		// std::cerr << "Axes almost aligned...\n";
		
	}else if( minor_dot < -1 + 1e-4 ){
		// std::cerr << "Axes almost anti-aligned...\n";

		// Idea: Just flip everything upside-down about the major axis

		rotate_ensemble( rot_data, maj_ax, N, PI );
		rotate_vector( min_ax, maj_ax, PI );
	}else{
		// std::cerr << "Axes not aligned...\n";
		

		ang2 = rotate_onto( min_ax, min_ax_t, rot_ax );
		rotate_ensemble( rot_data, rot_ax, N, ang2 );

		// Rotate axes as well:
		rotate_vector( maj_ax, rot_ax, ang2 );
		rotate_vector( min_ax, rot_ax, ang2 );
	}


	// std::cerr << "After rotation about minor axes, this is the state:\n";
	// std::cerr << "Source has major and minor axes:\n";
	// std::cerr << "     ( " << maj_ax[0] << ", " <<  maj_ax[1] << ", " << maj_ax[2] << " ).\n";
	// std::cerr << "     ( " << min_ax[0] << ", " <<  min_ax[1] << ", " << min_ax[2] << " ).\n";
	// std::cerr << "Template has major and minor axes:\n";
	// std::cerr << "     ( " << maj_ax_t[0] << ", " <<  maj_ax_t[1] << ", " << maj_ax_t[2] << " ).\n";
	// std::cerr << "     ( " << min_ax_t[0] << ", " <<  min_ax_t[1] << ", " << min_ax_t[2] << " ).\n";

	// std::cerr << "\n\n";
	// std::cerr << "Dot of major axes: " << dot(maj_ax,maj_ax_t) << ".\n";
	// std::cerr << "Dot of minor axes: " << dot(min_ax,min_ax_t) << ".\n";

	// Now determine the rmsd:
	py_float rmsd = 0.0;
	for( int i = 0; i < N; ++i ){
		py_float dx, dy, dz;
		
		dx = rot_data[i][0] - xt[0][0];
		dy = rot_data[i][1] - xt[0][1];
		dz = rot_data[i][2] - xt[0][2];
		
		double r2 = dx*dx + dy*dy + dz*dz;
		double mr2 = r2;
		for( int j = 1; j < N; ++j ){
			dx = rot_data[i][0] - xt[j][0];
			dy = rot_data[i][1] - xt[j][1];
			dz = rot_data[i][2] - xt[j][2];
		
			r2 = dx*dx + dy*dy + dz*dz;
			if ( r2 < mr2 ){
				mr2 = r2;
			}
		}
		rmsd += mr2;
	}

	if ( !xrot ){
		delete [] rot_data_;
	}

	return sqrt(rmsd);
}


void test_rotate_template()
{
	
}

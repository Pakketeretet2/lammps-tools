#ifndef DOMAIN_H
#define DOMAIN_H

#include "types.h"

/*!
  \file  domain.h
  @brief Contains functions for dealing with (periodic) domain boundaries.

  \ingroup cpp_lib
*/


/*!
  @brief Encoding of periodicity, so that an int of PERIODIC_X + PERIODIC_Z
  (1 + 4) indicates periodicity in X and Z but not in Y.
*/
enum PERIODICITIES {	
	PERIODIC_NONE = 0,
	PERIODIC_X    = 1,
	PERIODIC_Y    = 2,
	PERIODIC_Z    = 4,
	PERIODIC_FULL = PERIODIC_X + PERIODIC_Y + PERIODIC_Z
};


/*!
  @brief Computes distance between two points in a periodic box

  @param r        Array to store x1 - x2 in
  @param x1       First point
  @param x2       Second point
  @param xlo      Lower bounds of box
  @param xhi      Upper bounds of box
  @param periodic Int that encodes which boundaries are periodic
*/
void distance_wrap( py_float *r, const py_float *x1, const py_float *x2,
                    const py_float *xlo, const py_float *xhi,
                    py_int periodic = 7);

/*!
  @brief Computes distance between two points in a non-periodic box
  
  @param r  Array to store x1 - x2 in
  @param x1 First point
  @param x2 Second point
*/
void distance( py_float *r, const py_float *x1, const py_float *x2 );


#endif /* DOMAIN_H */

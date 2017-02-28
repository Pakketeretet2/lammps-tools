"""!
This module contains some code for finding hexagonal order with respect to some
axis and for stuff on a sphere.

\ingroup lammpstools
"""

import numpy as np
import scipy
import sys, os
import lammpstools
import math


def sign(a):
    return (a > 0) - (a < 0)

def get_theta( x, axis1, axis2 ):
    """ !Computes the rotation theta of the vector x about the given axis. """
    # Determine part along and perpendicular to axis:
    xa = np.dot(x,axis1) * axis1 /np.dot(axis1,axis1)
    xb = np.dot(x,axis2) * axis2 /np.dot(axis2,axis2)
    rp = x - xa - xb

    # print "x, axis1, axis2 = ", x, axis1, axis2

    # Determine sign of angle from cross product:
    cr = np.cross( axis1, x )
    d  = np.dot(x, axis1 ) / (np.linalg.norm(x) * np.linalg.norm(axis1) )
    if   d >  1: d = 1
    elif d < -1: d = -1
    
    ac = math.acos(d)
    
    if np.dot( cr, axis2 ) >= 0:
        theta = ac
    else:
        theta = 2*math.pi - ac

    # print theta
    return theta


def compute_psi_n( b, axis, rc, order, axis2 = np.array( [0.0, 0.0, 1.0] ),
                   dims = 3, itype = 0, jtype = 0, method = 0,
                   neighs = None, quiet = True ):
    """ !Computes the Psi_n order parameter for each atom, and the average.
        Constructs a nearest neighbour list if none given. """

    if neighs == None:
        neighs = lammpstools.neighborize( b, rc, dims, method, itype, jtype )

    im = lammpstools.make_id_map(b.ids)

    # Normalise axis:
    axlen = math.sqrt( np.dot(axis, axis) )
    axis *= 1.0 / axlen
    N = len(b.ids)

    psis = np.zeros( [N, 3] ) # Store real, imag and abs

    # For each atom, loop over its neighbours.
    for n in neighs:
        aid = n[0]
        idx = im[aid]

        xi = b.ids[idx]
        p  = [0.0, 0.0] # Complex number as list

        if len(n) > 1:
            nn = 1.0 / float(len(n) - 1)
        else:
            nn = 1.0
        
        for ajd in n[1:]:
            jdx = im[ajd]
            
            dr = b.x[idx] - b.x[jdx]
            theta = get_theta( dr, axis, axis2 )

            px = math.cos( order * theta )
            py = math.sin( order * theta )
            p[0] += px * nn
            p[1] += py * nn
            if not quiet: print("p += (%f, %f)" % (px, py), file = sys.stderr)
            

        psis[idx][0] = p[0]
        psis[idx][1] = p[1]
        psis[idx][2] = math.sqrt(p[0]*p[0] + p[1]*p[1])
        if not quiet:
            print("Neighs of %d = %d, 1 / nn = %f" % (idx, len(n)-1, nn), file = sys.stderr)

    psi_avg = np.zeros(3)
    n = 1.0 / float(len(psis[:,0]))
    for i in b.ids:
        idx = im[i]
        if not quiet:
            print("psi_6 of %d = (%f, %f), |psi_6| = %f" % \
                  (i, psis[idx][0], psis[idx][1], psis[idx][2]), file = sys.stderr)
        # print psis[idx,0], psis[idx,1], psis[idx,2]
        
        psi_avg[0] += n * psis[idx,0]
        psi_avg[1] += n * psis[idx,1]
        psi_avg[2] += n * psis[idx,2]
        

    
    return psi_avg, psis

"""!
This module contains some functions to calculate the centre of mass

\ingroup lammpstools
"""

import dumpreader
import lammpstools
import numpy as np
import scipy
import math
import sys
from typecasts import *
from ctypes import *


"""
  Computes the centre of mass for periodic box
"""
def com_periodic   ( x, x0, x1, types, masses, filt_indices ):
    xcom = 0.0
    M    = 0.0

    xi_avg = zeta_avg = 0.0
    Lx = x1 - x0
    xavg = 0.0

    xk = zk = 0.0
    
    print >> sys.stderr, "[x0, x1] = [%g, %g]" % (x0,x1)
    
    for i in filt_indices:

        t  = types[i]
        xx = (x[i] - x0) / Lx
        if xx < 0:
            xx += 1.0
        elif xx > 1:
            xx -= 1.0
        m  = masses[t-1]
        M  += m;

        
        theta = xx * 2.0 * math.pi
        if xx < 0 or xx > 1:
            print >> sys.stderr, "theta out of bounds! xx = %g" % xx
        
        xk += m*math.cos( theta )
        zk += m*math.sin( theta )
        
    xk /= M
    zk /= M
    
    theta_avg = math.atan2( -zk, -xk ) + math.pi;
    avg = theta_avg * Lx / (2.0*math.pi)
    xavg = avg + x0


    # Check to make sure it's in bounds:
    if xavg < x0 or xavg > x1:
        print >> sys.stderr, "<x> is out of bounds!"
    
    return xavg

"""
  Computes the centre of mass for non-periodic box

"""
def com_nonperiodic( x, types, masses, filt_indices ):
    xcom = 0.0
    M    = 0.0
    
    for i in filt_indices:
        t  = types[i]
        xx = x[i]
        m  = masses[t-1]
        xcom += xx*m;
        M    += m;
    return xcom / M

"""
  Computes the centre of mass of the block_data b.
  Can be filtered by passing a container with the indices to calculate the com of.

  @param b             Block data
  @param masses        Masses per atom type (indexed by type - 1)
  @param filt_indices  Indices to compute com for (all used if not given)
"""
def compute_com( b, masses, filt_indices = None ):
    if filt_indices == None:
        filt_indices = range(0, b.meta.N)

    x_com = y_com = z_com = 0.0

    x = b.x[:,0]
    y = b.x[:,1]
    z = b.x[:,2]

    types  = b.types
    x0 = b.meta.domain.xlo[0]
    x1 = b.meta.domain.xhi[0]
    y0 = b.meta.domain.xlo[1]
    y1 = b.meta.domain.xhi[1]
    z0 = b.meta.domain.xlo[2]
    z1 = b.meta.domain.xhi[2]
    
    
    if b.meta.domain.periodic & 1:
        x_com = com_periodic( x, x0, x1, types, masses, filt_indices );
    else:
        x_com = com_nonperiodic( x, types, masses, filt_indices );
    
    if b.meta.domain.periodic & 2:
        y_com = com_periodic( y, y0, y1, types, masses, filt_indices );
    else:
        y_com = com_nonperiodic( y, types, masses, filt_indices );        

    if b.meta.domain.periodic & 4:
        z_com = com_periodic( z, z0, z1, types, masses, filt_indices );
    else:
        z_com = com_nonperiodic( z, types, masses, filt_indices );

    return np.array( [x_com, y_com, z_com] )

    

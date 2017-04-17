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
    
    # print("[x0, x1] = [%g, %g]" % (x0,x1), file = sys.stderr)
    
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
            print("theta out of bounds! xx = %g" % xx, file = sys.stderr)
        
        xk += m*math.cos( theta )
        zk += m*math.sin( theta )
        
    xk /= M
    zk /= M
    
    theta_avg = math.atan2( -zk, -xk ) + math.pi;
    avg = theta_avg * Lx / (2.0*math.pi)
    xavg = avg + x0


    # Check to make sure it's in bounds:
    if xavg < x0 or xavg > x1:
        print("<x> is out of bounds!", file = sys.stderr)
    
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

  @param b       Block data
  @param masses  Masses per atom type (indexed by type - 1)
  @param groups  Indicates which atoms are in same group. Group 0 means ignore.
  @param dims    Number of dimensions of the system (2 or 3)
"""
def compute_com( b, masses, groups, dims  ):
    return compute_com_cpp( b, masses, groups, dims )


"""
  Computes the centre of mass for groups of atoms in the block_data b.

  @param b         Block data
  @param masses    Masses per atom type (indexed by type - 1)
  @param groups    Indicates which atoms should be grouped together.
  @param dims      Number of dimensions of the system (2 or 3)
"""    
def compute_com_cpp( b, masses, groups, dims ):

    Ngroups = max(groups)

    com = np.zeros( [Ngroups + 1, dims], dtype = float )
    
    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
    lammpstools.center_of_mass( void_ptr(b.x), c_longlong(b.meta.N),
                                void_ptr(b.ids), void_ptr(b.types),
                                c_longlong(Ngroups), void_ptr(groups),
                                void_ptr(masses), void_ptr(b.meta.domain.xlo),
                                void_ptr(b.meta.domain.xhi),
                                c_longlong(b.meta.domain.periodic),
                                c_longlong(dims), void_ptr( com ) )
    
    return com

"""!
\defgroup lammpstools
This module contains some functions that try to find the minimum RMSD between two sets of points.

@package lammpstools

\ingroup lammpstools
"""

import lammpstools
import dumpreader
import numpy as np
import math
import sys

from ctypes import *
from typecasts import *


def icp( b_a, b_b, dims = 3, rot_mat_guess = np.eye( 3, dtype = float ) ):
    """! Finds the smallest rmsd between block a and b using iterative closest point. """
    N = b_a.meta.N
    if N != b_b.meta.N:
        print >> sys.stderr, "Block sizes must match!"
        return None

    # Leverage everything to libICP.

    # Compute initial rmsd yourself:
    rmsd0 = 0.0
    for i in range( 0, b_a.meta.N ):
        xi = b_a.x[i]
        xj = b_b.x[i]
    
    rmsd, rotated = call_icp_lib( b_a.x, b_b.x, N, b_a.meta.domain, dims )

    return rmsd, rotated


def call_icp_lib( xi, xj, N, domain, dims = 3, rot_mat = np.eye( 3, dtype = float ) ):
    """! Forwards call to the icp-wrapper in c++-lib. """

    rotated = np.zeros( [ N, 3 ], dtype = float )

    
    
    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
    cpp_icp_rot = lammpstools.icp_rot
    lammpstools.icp_rot.restype = ctypes.c_double
    
    rmsd = cpp_icp_rot( void_ptr(xi), void_ptr(xj), c_longlong(N),
                        c_longlong(domain.periodic),
                        void_ptr(domain.xlo), void_ptr(domain.xhi),
                        c_longlong(dims), void_ptr(rot_mat),
                        void_ptr(rotated) )
    return rmsd, rotated

def minimum_rmsd_rotate( xi, xj, N, ids, jds, periodic, xlo, xhi, dims ):
    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")

    rotated = np.zeros( [N, 3], dtype = float )
    axis    = np.zeros( 3, dtype = float )
    dist2 = lammpstools.map_clusters( void_ptr(xi), void_ptr(xj), c_longlong(N),
                                      void_ptr(ids), void_ptr(jds), c_longlong(periodic),
                                      void_ptr( xlo ), void_ptr( xhi ), c_longlong(dims),
                                      void_ptr( rotated ), void_ptr( axis ) )
    
def rotate_to_template(b, template, ax1, ax2, ref_axis1, ref_axis2,
                       rotated = None ):

    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
    lammpstools.rotate_to_template.restype = ctypes.c_double
    if (rotated is None): rot_ptr = void_ptr(0)
    else:                 rot_ptr = void_ptr( rotated.x )

    rmsd =  lammpstools.rotate_to_template( void_ptr(b.x), void_ptr(template.x), b.meta.N,
                                            b.meta.domain.periodic, void_ptr(b.meta.domain.xlo),
                                            void_ptr(b.meta.domain.xhi), 3, void_ptr(ax1),
                                            void_ptr(ax2), void_ptr(ref_axis1), void_ptr(ref_axis2),
                                            rot_ptr )
    
    return rmsd

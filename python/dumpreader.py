"""! @package dumpreader
This package contains functions for extracting dump files.

\ingroup lammpstools
"""

import numpy as np
import gzip
import os
import sys
import copy

import dumpreader_lammps as lammps



def read_raw_pts(fname):
    """! Reads a raw text file, assuming coords are x y z in first three columns. """
    
    xyz = np.loadtxt(fname)
    N = len(xyz[:,0])
    t = 0

    xmin = np.min(xyz[:,0])
    xmax = np.max(xyz[:,0])
    ymin = np.min(xyz[:,1])
    ymax = np.max(xyz[:,1])
    zmin = np.min(xyz[:,2])
    zmax = np.max(xyz[:,2])

    ids = np.zeros(N)
    for i in range(0,N):
        ids[i] = i+1
    types = np.ones(N)
    
    xlo = np.array([xmin-2,ymin-2,zmin-2],dtype=np.float64)
    xhi = np.array([xmax+2,ymax+2,zmax+2],dtype=np.float64)
    dom = domain_data(xlo, xhi, 0, "")

    b = block_data(t=0,N=N,domain=dom,ids=ids,types=types,x=xyz)
    return b







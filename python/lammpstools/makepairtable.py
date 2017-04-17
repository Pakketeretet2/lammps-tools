"""!
This file contains some pair potentials.

\ingroup lammpstools
"""


import lammpstools
import dumpreader
import numpy as np
import math
import sys

def make_pair_table( fname, name, pair_pot, N, mode = "R", lo = 1.0, hi = 10.0  ):
    "Dumps a LAMMPS-style pair table to given file."
    
    if mode == "R":
        dr = (hi - lo)/(N-1)
    elif mode == "RSQ":
        print >> sys.stderr, "Mode RSQ not supported!"
        return -1
    elif mode == "BITMAP":
        print >> sys.stderr, "Mode BITMAP not supported!"
        return -1
    else:
        print >> sys.stderr, "Mode ", mode, " not recognized!"
        return -1

    # First test the given potential:
    if lammpstools.test_potential( pair_pot, lo, hi, 1e-4, 1e-8 ):
        print >> sys.stderr, "Potential not consistent!"
        # return -2

    # Fill table:
    fp = open(fname,"w")
    use_fprime = False
    if hasattr( pair_pot, "force_prime" ):
        "use_fprime = True"

    print >> fp, name
    if use_fprime:
        print >> fp, "N %d %s %f %f" % (N, mode, lo, hi)
    else:
        fplo = pair_pot.force_prime(lo)
        fphi = pair_pot.force_prime(hi)
        print >> fp, "N %d %s %e %e FPRIME %f %f" % (N, mode, lo, hi, fplo, fphi)
    
    print >> fp, ""

    for i in range(0,N):
        r = lo + i*dr
        E  = pair_pot.energy(r)
        f  = pair_pot.force(r)
            
        print >> fp, "%d %e %e %e" % (i,r,E,f)
            

    return 0

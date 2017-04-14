"""! 
Contains some tools for bond analysis

\ingroup lammpstools

"""


import numpy as np
import scipy
import sys, os
import lammpstools
import math

## Computes the \p psi_n order parameter for each atom, and the average.
#
#  \param b       block_data to determine \p psi_n for
#  \param axis    Axis with respect to which \p psi_n should be calculated
#  \param order   Order of \f$\Psi_n\f$, i.e., the \f$n\f$ in \f$\Psi_n\f$
#  \param normal  Out of plane axis, i.e., each bond vector is
#                 projected onto the plane with this normal.
#  \param rc      Cut-off for constructing neighbour list.
#  \param itype   Consider only bonds between types \p i...
#  \param jtype   ... and \p j.
#  \param method  Method used for neighbourising.
#  \param neighs  Pre-calculated neighbour list. This saves computation time.
#  \param quiet   Be quiet?
#
#  \returns The average bond order parameter \f$\left<|\Psi_n|\right>\f$
#           and a numpy array of the individual \f$\Psi_n\f$ values.
#           The entries are in principle complex, the array contains the
#           real and imaginary parts and the absolute value, in that order.
#
# Calculates the \f$\Psi_n\f$ bond order parameter for all particles.
def compute_psi_n( b, axis, order, normal, rc, 
                   itype = 0, jtype = 0, method = 0,
                   neighs = None, quiet = True ):
    """ Calculates the \Psi_n order parameter. """
    if neighs is None:
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
            theta = get_theta( dr, axis, normal )

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



## Computes the rotation \p theta of the vector x about
#  axis1 in plane defined by normal.
#
#  \param x       Vector whose angle to calculate w.r.t. \p axis
#  \param axis    Axis that defines angle of x.
#  \param normal  Defines the plane in which to calculate the angle.
#
def get_theta( x, axis, normal ):
    """ Determine part along and perpendicular to axis. """
    xa = np.dot(x,axis) * axis /np.dot(axis,axis)
    xb = np.dot(x,normal) * normal /np.dot(normal,normal)
    rp = x - xa - xb

    # print "x, axis1, axis2 = ", x, axis1, axis2

    # Determine sign of angle from cross product:
    cr = np.cross( axis, x )
    d  = np.dot(x, axis ) / (np.linalg.norm(x) * np.linalg.norm(axis) )
    if   d >  1: d = 1
    elif d < -1: d = -1
    
    ac = math.acos(d)
    
    if np.dot( cr, normal ) >= 0:
        theta = ac
    else:
        theta = 2*math.pi - ac

    # print theta
    return theta



## Calculates info on bond statistics.
#
#  \param b       block_data to determine statistics for
#  \param rc      Cut-off for constructing the neighbour list, if any
#  \param itype   Consider only bonds between types \p itype...
#  \param jtype   ... and \p jtype (0 for any).
#  \param method  Method used for neighbourising.
#  \param neighs  Pre-calculated neighbour list. This saves computation time.
#
#  This function checks bond_data for bonds formed between pairs of types
#  \p itype and \p jtype. It returns a list of bonds and their distance.
# 
def get_bond_stats( b, rc, itype = 0, jtype = 0, method = 0, neighs = None ):
    """ Calculates statistics on the bonding of atoms of type i and j. """
    if neighs is None:
        neighs = lammpstools.neighborize( b, rc, 3, method, itype, jtype )

    im    = lammpstools.make_id_map( b.ids )
    bonds = []
    for ni in neighs:
        idi = ni[0]
        idx = im[idi]
        if itype and (b.types[idx] != itype):
            continue

        xi = b.x[idx]
        for idj in ni[1:]:
            jdx = im[idj]
            xj = b.x[jdx]
            rij = lammpstools.distance( b, xi, xj )

            bonds.append( idi, idj, math.sqrt( np.dot(rij, rij) ) )
    return bonds
        


## Determines the number of topological defects for each particle and
#  returns a dump_col containing them per particle in block.
# 
#  \param b    Block to determine the topological charges of
#  \param n    Neighbor list corresponding to given block.
#  \param nn0  The expected number of neighbours, i.e., particles with
#              number of neighbours == nn0 will have charge 0.
# 
def topological_defects( b, n, nn0 = 6 ):
    Qs = np.zeros(b.meta.N, dtype=int)
    i = 0
    Qtot = 0
    for ni in n:
        nns = len(ni) - 1
        q   = nn0 - nns
        Qs[i] = q
        i += 1
        Qtot += q
    
    dc = dump_col("q", Qs)

    return dc, Qtot

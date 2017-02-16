"""!
This file contains some analysis stuff related to melting.

\ingroup lammpstools
"""

import lammpstools
import dumpreader
import numpy as np
import math
import sys

import matplotlib.pyplot as plt

def get_avg_shortest_dist( b ):
    """ !Finds the shortest distance for each particle to another. Expensive! """
    N = b.meta.N
    dists = np.ones(N) * 100000
    for i in range(0,N):
        for j in range(0,N):
            if i == j: continue
            
            rij, rvec = lammpstools.distance( b, i, j )
            if rij < dists[i]:
                dists[i] = rij
    a0 = np.average(dists)
    
    dists = []
    for i in range(0,N):
        for j in range(0,N):
            if i == j: continue
            
            rij, rvec = lammpstools.distance( b, i, j )
            if rij < 1.2*a0:
                dists.append(rij)

    a = np.average(dists)
    fp = open("distances.dat","w")
    for r in dists:
        print >> fp, r
    return a



def lindemann( d, rc = None, a = None, method = None, pairs = None ):
    """ !Performs Lindemann-like analysis on given dump file. """
    
    b = d.getblock()
    if rc is None:
        if a is None:
            # Extract a sane a and rc
            a = get_avg_shortest_dist( b )

        rc = 1.25*a

    else:
        if a is None:
            # Find reasonable a from first snapshot
            n = lammpstools.neighborize( b, rc, 3, method = method )
            id_map = lammpstools.make_id_map( b.ids )
            dists = []
            for ni in n:
                idxi = id_map[ ni[0] ]
                for idj in ni[1:]:
                    idxj = id_map[ idj ]
                    rij, rvec = lammpstools.distance( b, idxi, idxj )
                    dists.append(rij)

            fp = open("distances.dat","w")
            for r in dists:
                print >> fp, r
            a = np.average(dists)    
    # Now rc and a have appropriate values.
    print >> sys.stderr,  "Using a = %f and rc = %f" % (a, rc)
    a2 = a*a

    if pairs is None:
        # Store the initial neighbour lists.
        n0     = lammpstools.neighborize( b, rc*1.12, 3, method = method )

        # Convert the neighbour list n0 into a list of pairs.
    
        pairs = []
        for ni in n0:
            i = ni[0]
            for j in ni[1:]:
                if j < i:
                    j, i = i, j
                if [ i, j ] not in pairs:
                    id_map = lammpstools.make_id_map( b.ids )
                    idx = id_map[ i ]
                    jdx = id_map[ j ]
                    rij, rvec = lammpstools.distance( b, idx, jdx )
                    if rij < rc:
                        pairs.append( [ i, j ] )

    Np = len(pairs)
    print >> sys.stderr, "Found %d pairs." % Np
    for p in pairs:
        idi = p[0]
        idj = p[1]
        id_map = lammpstools.make_id_map( b.ids )
        idxi = id_map[ idi ]
        idxj = id_map[ idj ]
        rij, rvec = lammpstools.distance( b, idxi, idxj )
    
    distances = []
    times     = []
    frames = 0
    
    while not d.at_eof:
        dists  = np.zeros( Np )
        gammaL = np.zeros( Np )
        id_map = lammpstools.make_id_map( b.ids )
        times.append( b.meta.t )
        for i, pair in zip( range(0, Np), pairs ):
            if not ( (pair[0] in b.ids ) and (pair[1] in b.ids) ):
                print >> sys.stderr, "Lost pair (%d, %d)!" % (pair[0], pair[1])

            
            j = id_map[ pair[0] ]
            k = id_map[ pair[1] ]

            # You really just want the varaince of rij for each pair.
            rij, rvec = lammpstools.distance( b, j, k )
            dists[i] = rij

        distances.append( dists )
        frames += 1

        if frames % 100 == 0:
            print >> sys.stderr, "At block %d" % frames

        b = d.getblock()

    print >> sys.stderr, "Got %d blocks to compute Lindemann index for." % frames

    distances_in_time = np.zeros( [ len(distances), Np ] )
    for d, nt in zip(distances,range(0,len(distances))):
        for i in range(0,Np):
            distances_in_time[nt, i] = d[i]


    gamma = np.zeros( [ len(distances), Np ] )
    gamma_global = np.zeros( len(distances) )
    # Compute lindemann index:
    c = 1.0 / (2.0*a*a)
    for t, nt in zip( times, range(0,len(times)) ):
        for pp in range(0,Np):
            dunow = distances_in_time[nt,pp]
            du0   = distances_in_time[0,pp]
            gamma[nt,pp] = dunow*dunow + du0*du0 - 2*(dunow - du0)
            gamma[nt,pp] *= c
        gamma_global[nt] = np.average( gamma[nt,:] )

    return gamma_global, gamma, times


"""!
This module contains some normal mode analysis tools

\ingroup lammpstools
"""

import dumpreader
import lammpstools
import numpy as np
import scipy
import sys, os
import math
from typecasts import *
from ctypes import *


from multiprocessing import Process

def normal_mode_analysis(d):
    """ !Performs normal mode analysis on given dump file """

    b = d.getblock()
    print("Initially, t = ", b.meta.t)
    B = []

    while not d.at_eof:
        if b is None:
            break
        B.append(b)
        b = d.getblock()
    
    natoms = b.meta.N
    N = 3*len(b.x)
    Xavg = np.zeros(N)
    im = lammpstools.make_id_map(b.ids)
    nframes = 1
    for b in B:
        X = np.zeros(N)
        for i in range(0,b.meta.N):
            idi = b.ids[i]
            j = idi - 1
            X[3*j]   = b.x[i,0]
            X[3*j+1] = b.x[i,1]
            X[3*j+2] = b.x[i,2]

        Xavg += X
        nframes += 1
    Xavg /= float(nframes)

    print("Average X: ", X)

    # Set up the covariance matrix:
    C = np.zeros([N,N])
    

    for b in B:
        for i in range(0,b.meta.N):
            idi = b.ids[i]
            j = idi - 1
            dxi = b.x[i,0] - Xavg[3*j]
            dyi = b.x[i,1] - Xavg[3*j+1]
            dzi = b.x[i,2] - Xavg[3*j+2]
                        
            for k in range(0,b.meta.N):
                idk = b.ids[k]
                l = idk - 1
                
                dxk = b.x[k,0] - Xavg[3*l]
                dyk = b.x[k,1] - Xavg[3*l+1]
                dzk = b.x[k,2] - Xavg[3*l+2]

                C[3*j  ,3*l  ] += dxi*dxk
                C[3*j  ,3*l+1] += dxi*dyk
                C[3*j  ,3*l+2] += dxi*dzk
                
                C[3*j+1,3*l  ] += dyi*dxk
                C[3*j+1,3*l+1] += dyi*dyk
                C[3*j+1,3*l+2] += dyi*dzk
                
                C[3*j+2,3*l  ] += dzi*dxk
                C[3*j+2,3*l+1] += dzi*dyk
                C[3*j+2,3*l+2] += dzi*dzk
        
    C /= nframes

    w, v = np.linalg.eig(C)

    
    permutation = sorted(range(len(w)), key = lambda x: -w[x])
    
    # Print the largest nlarge:

    nlarge = b.meta.N*3
    print("The largest %d eigenvalues and eigenvectors are:", nlarge)
    for i in range(0,nlarge):
        print("%d: %f, " % (i, w[i]), end = "")
        print(v[:,i])

    fp = open( "eigenvalues_old.dat", "w" )
    for i in range(0,len(w)):
        print("%f" % w[i], file = fp )
    
    # Sort in order of highest to lowest mode:
    w     = w[permutation]
    tmp_v = v
    for k in range(0,len(permutation)):
        v[:,k] = tmp_v[:,permutation[k]]


    dump_out = "lowest_mode.dump"
    for i in range(0, b.meta.N):
        b.ids[i] = i+1
        
    for k, mode in zip(range(0,len(w)), w):
        
        # Output some mode of motion:
        Xmod = np.zeros([b.meta.N,3])
        M    = 50
        times = np.linspace(0,4*math.pi,M)
        bm = dumpreader.block_data( b.meta, b.ids, b.types, np.zeros( [N,3] ) )
        v0 = v[:,k]
        if mode < 0:
            print("WTF, mode %g is negative?" % mode, file = sys.stderr)
        if math.fabs(mode) <= 1e-12:
            mode = 0.0
        A  = math.sqrt(mode)
        print("Eigenvalue %d is %f" % (k, mode))

        for t, n in zip(times, range(0,len(times))):

            bm.meta.t = n
            st = math.sin(t)
            for i in range(0, b.meta.N):
                bm.x[i,0] = Xavg[3*i]   + A*st*v0[3*i]
                bm.x[i,1] = Xavg[3*i+1] + A*st*v0[3*i+1]
                bm.x[i,2] = Xavg[3*i+2] + A*st*v0[3*i+2]

            dumpreader.block_to_dump( bm, dump_out, "a" )

def normal_mode_analysis_clib( d ):
    pass_blocks_to_clib( d )
            
def pass_blocks_to_clib( d ):
    """! Writes block data from a dump file to a pipe. """

    b = d.getblock()
    if b is None:
        print("Block was None, probably dump file %s does not exist?" % dump, file = sys.stderr)
        return

    
    lammpstools = cdll.LoadLibrary("liblammpstools.so")
    pname = '/tmp/lammpstools_neighborize_pipe_' + str(os.getpid())
    
    pipe_status = lammpstools.get_pipe( pname )
    if pipe_status != 0:
        print("Pipe status was not 0 b.meta.t ", pipe_status, file = sys.stderr)
    else:
        print("Pipe named %s succesfully opened." % pname, file = sys.stderr)


    # Make storage for eigenvectors and eigenvalues
    N = 3*b.meta.N
    eigenvalues  = np.zeros(N)
    eigenvectors = np.zeros([N, N])

    
    def start_normal_mode_analysis():
        lammpstools.normal_mode_analysis( pname, void_ptr(eigenvalues),
                                          void_ptr(eigenvectors), c_longlong(N) )
    p = Process( target = start_normal_mode_analysis )
    p.start()
    fp  = open( pname, "w" )
    fp2 = open( "pipe_test", "w" )

    while not d.at_eof:

        if b.meta.N <= 0:
            break
        
        print( b.meta.N, file = fp )
        print( b.meta.N, file = fp2 )
        
        X_sorted = np.zeros( [b.meta.N, 3] )
        
        for i in range(0,b.meta.N):
            idi = b.ids[i]
            j = idi - 1
            X_sorted[j][0] = b.x[i][0]
            X_sorted[j][1] = b.x[i][1]
            X_sorted[j][2] = b.x[i][2]

        for i in range(0,b.meta.N):
            print( "%1.16f %1.16f %1.16f" % (X_sorted[i][0], X_sorted[i][1], X_sorted[i][2]), file = fp )
            print( "%1.16f %1.16f %1.16f" % (X_sorted[i][0], X_sorted[i][1], X_sorted[i][2]), file = fp2 )
            
        
        b = d.getblock()
    
    fp.close()
    pipe_status = lammpstools.close_pipe( pname )

    print("Eigenvalues:")
    for v in eigenvalues:
        print("%f" % v)
    
    

    
    if pipe_status != 0:
        print( "Pipe status was not 0 b.meta.t ", pipe_status, file = sys.stderr )
    else:
        print( "Pipe named %s succesfully closed." % pname, file = sys.stderr )
    

"""!
\file util.py

Contains some utility routines.
\ingroup lammpstools
"""

import numpy as np, sys
from ctypes import *
from lammpstools.typecasts import *
from lammpstools import block_data, block_meta


## Makes an id map in an associative array.
#  The map contains at id_map[ iid ] the index of
#  atom data belonging to  atom with id = iid.
#
#  \param b  The array of ids to make the map out of.
#
def make_id_map( ids, largest_id = None ):
    if largest_id is None:
        Nids = np.max(ids)
        largest_id = Nids

    id_map = -np.ones( largest_id + 1, dtype = int )
    for i in range(0, len(ids)):
        iid = ids[i]
        id_map[ iid ] = i
    return id_map





## Generates an ensemble of particles on given manifold.
# 
# @param N              Number of particles to create
# @param type           Type of the particles to create
# @param domain         Domain to generate particles in
# @param manifold_name  Name of the manifold to generate particles on
# @param manifold_args  List of arguments to the manifold to create
# @param rc             Minimum distance b.meta.tween created particles. 0 by default
# @param seed           Random seed to use. 1 by default
# 
def generate_ensemble_on_manifold( N, typ, domain, manifold_name, manifold_args, rc = 0.0, seed = 1 ):
    x       = np.zeros( [N, 3], dtype = np.float64 )
    ids_arr = np.zeros( N,      dtype = int )
    types   = np.zeros( N,      dtype = int )

    manifold_arg_str = ""
    for s in manifold_args:
        manifold_arg_str += s + " "

    print("Manifold arg string = %s" % manifold_arg_str, file = sys.stderr)

    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
    # Pass the strings as string buffers now.
    #manifold_name_buff    = 
    #manifold_arg_str_buff = 
    lammpstools.ensemble_generate_manifold( void_ptr(x), N, void_ptr(ids_arr),
                                            void_ptr(types), c_longlong(domain.periodic),
                                            void_ptr(domain.xlo), void_ptr(domain.xhi),
                                            c_char_p(manifold_name.encode('ascii')),
                                            c_char_p(manifold_arg_str.encode('ascii')),
                                            c_double(rc), c_longlong(seed) )
    for i in range(0,N):
        types[i] = typ

    
    # Convert the raw to a block-data struct:
    meta = block_meta( 0, N, domain )
    return block_data( meta, ids_arr, types, x )
    


## Computes the shortest distance between atoms indices i and j in block.
#
#  \param b block_data to calculate distance in.
#  \param i Index of particle i
#  \param j Index of particle j
#  
def distance( b, i, j ):
    """ ! """
    xi = b.x[i]
    xj = b.x[j]

    return domain_distance( b.meta.domain, xi, xj )

## Computes the shortest distance between atoms indices i and j
#  in different blocks
#
#  \param b block_data to calculate distance in.
#  \param i Index of particle i
#  \param j Index of particle j
#  
def distance_diff_blocks( b1, b2, i, j ):
    """ !Computes the shortest distance between atoms i in block b1 and j in block b2. """
    xi = b1.x[i]
    xj = b2.x[j]
    return domain_distance( b1.meta.domain, xi, xj )



def domain_distance( domain, xi, xj ):
    """! Calculates the distance between xi and xj according to domain. """
    r  = xi - xj
    Lx = domain.xhi[0] - domain.xlo[0]
    Ly = domain.xhi[1] - domain.xlo[1]
    Lz = domain.xhi[2] - domain.xlo[2]
    
    if domain.periodic & 1:
        if r[0] >  0.5*Lx:
            r[0] -= Lx
        elif r[0] <= -0.5*Lx:
            r[0] += Lx
    
    if domain.periodic & 2:
        if r[1] >  0.5*Ly:
            r[1] -= Ly
        elif r[1] <= -0.5*Ly:
            r[1] += Ly
    
    if domain.periodic & 4:
        if r[2] >  0.5*Lz:
            r[2] -= Lz
        elif r[2] <= -0.5*Lz:
            r[2] += Lz
    dist = np.linalg.norm(r)
    
    return dist, r


def block_filter( block, indices ):
    """ !Returns a block containing only the listed indices. """
    N = len(indices)
    x     = np.zeros( [N, 3], dtype = float )
    ids   = np.zeros(N, dtype=int)
    types = np.zeros(N, dtype = int)

    Ncols = len(block.other_cols)

    other_cols = []
    for i in range(0,Ncols):
        other_cols.append( dumpreader.dump_col( block.other_cols[i].header, []) )
    
    j = 0
    for i in indices:
        x[j] = block.x[i]
        ids[j] = block.ids[i]
        types[j] = block.types[i]

        # copy other cols:
        for k in range(0,Ncols):
            other_cols[k].data.append( block.other_cols[k].data[i] )
        
        j += 1
 
    bmod = block_data( block.meta, ids, types, x )
    bmod.meta.N = len(indices)
    bmod.other_cols = other_cols

    for k in range(0,Ncols):
        bmod.other_cols[k].N = len(bmod.other_cols[k].data)
    return bmod



def recenter( b, rescale_box = False ):
    """ !Recenters particles in block so that total COM is at (0,0,0). """
    com = np.zeros(3,dtype=float)
    for i in range(0,b.meta.N):
        com += b.x[i]

    c = 1.0 / b.meta.N
    for i in range(0,b.meta.N):
        b.x[i][0] -= com[0]*c
        b.x[i][1] -= com[1]*c
        b.x[i][2] -= com[2]*c

    print(b.meta.domain.xlo)
    print(b.meta.domain.xhi)
    
    b.meta.domain.xlo[0] -= com[0]*c
    b.meta.domain.xlo[1] -= com[1]*c
    b.meta.domain.xlo[2] -= com[2]*c
    b.meta.domain.xhi[0] -= com[0]*c
    b.meta.domain.xhi[1] -= com[1]*c
    b.meta.domain.xhi[2] -= com[2]*c

    if rescale_box:
        xmax = np.max( b.x[:,0] )
        xmin = np.min( b.x[:,0] )
        ymax = np.max( b.x[:,1] )
        ymin = np.min( b.x[:,1] )
        zmax = np.max( b.x[:,2] )
        zmin = np.min( b.x[:,2] )

        Lx = 1.2*(xmax - xmin)
        Ly = 1.2*(ymax - ymin)
        Lz = 1.2*(zmax - zmin)

        b.meta.domain.xlo[0] = -0.5*Lx
        b.meta.domain.xhi[0] =  0.5*Lx
        b.meta.domain.xlo[1] = -0.5*Ly
        b.meta.domain.xhi[1] =  0.5*Ly
        b.meta.domain.xlo[2] = -0.5*Lz
        b.meta.domain.xhi[2] =  0.5*Lz
        
        return b
    else:
        return b




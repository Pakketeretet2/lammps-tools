"""!
\defgroup lammpstools
This module contains all Python-side tools for analysis of LAMMPS data.

@package lammpstools
This package defines the python interface to the C++ lib and some
(higher level) analysis tools written in Python itself.

\ingroup lammpstools
"""

import numpy as np, struct, time, math, sys, os, multiprocessing

from ctypes    import *
from typecasts import *


# Import all submodules:
import histogram, normal_mode_analysis, melting_analysis, minimum_rmsd, \
    fit_einstein_crystal, multiprocessing, ribbon_analysis, compute_com, \
    bond_analysis, potentials, makepairtable, block_data


## Computes the RDF of atoms of types itype and jtype from block data b
#
#  \param b          Block of data to compute RDF for
#  \param r0         Inner cutoff radius where RDF is computed
#  \param r1         Outer cutoff radius where RDF is computed
#  \param nbins      Number of bins to use, so resolution \f$ dr =
#                    (r1 - r0)/(nbins-1) \f$
#  \param itype      Type of atoms 1 to consider (0 for all)
#  \param jtype      Type of atoms 2 to consider (0 for all)
#  \param dims       Dimension of simulation box (used in normalisation)
#  \param method     Neighborisation method to use.
#
def compute_rdf( b, r0, r1, nbins, itype, jtype, dims, method = None ):
    """ Computes the RDF of atoms of types itype and jtype for block_data b. """

    rdf    = np.zeros(nbins, dtype=np.float64);
    coords = np.zeros(nbins, dtype=np.float64);
    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
    if method is None:
        # Guess a good method based on b.meta.N:
        if b.meta.N < 200: method = 0
        else:         method = 1
    
    lammpstools.compute_rdf( void_ptr(b.x), c_longlong(b.meta.N), void_ptr(b.ids),
                             void_ptr(b.types), c_double(r0),
                             c_double(r1), c_longlong(nbins),
                             c_longlong(itype), c_longlong(jtype),
                             void_ptr(b.meta.domain.xlo), void_ptr(b.meta.domain.xhi),
                             c_longlong(b.meta.domain.periodic), c_longlong(dims),
                             c_longlong(method), void_ptr(rdf), void_ptr(coords) )

    pts = np.zeros(nbins,dtype=np.float64);
    for i in range(0,nbins):
        pts[i] = r0 + (r1-r0)*i/float(nbins-1)
    return pts, rdf, coords


## Computes the ADF of atoms of types itype and jtype from block data b
#  for particles on a sphere of radius R.
# 
# @param b          Block of data to compute ADF for
# @param R          Radius of template
# @param nbins      Number of bins to use, so resolution dr = (r1 - r0)/(nbins-1)
# @param itype      Type of atoms 1 to consider (0 for all)
# @param jtype      Type of atoms 2 to consider (0 for all)
# @param method     Neighborization method to use. Defaults to distance (either
#                   binned or not, depending on the size of b.meta.N)
#
def compute_adf( b, R, nbins, itype, jtype, method = None ):
    """ Computes ADF of atoms of types itype and jtype from block_data. """
    adf    = np.zeros(nbins, dtype=np.float64);
    coords = np.zeros(nbins, dtype=np.float64);
    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
    if method is None:
        print("Please set method explicitly.", file = sys.stderr)
        return None
    
    lammpstools.compute_adf( void_ptr(b.x), c_longlong(b.meta.N), void_ptr(b.ids),
                             void_ptr(b.types), c_longlong(nbins),
                             c_longlong(itype), c_longlong(jtype),
                             c_double(R), c_longlong(method), void_ptr(adf),
                             void_ptr(coords) )
    
    pts = np.zeros(nbins,dtype=np.float64);
    for i in range(0,nbins):
        pts[i] = math.pi * i / float(nbins-1)
    return pts, adf, coords


# Makes a neighbor list of all particles in block.
# 
# @param b         Block of data to neighborize
# @param rc        Cut-off for distance criterion (ignored if not needed)
# @param dims      DImensions of simulation box
# @param method    Neighborization method to use. 
# 
def neighborize( b, rc, dims, method = None, itype = 0, jtype = 0,
                 quiet = True ):
    if method is None:
        # Guess a good method based on b.meta.N:
        if b.meta.N < 200: method = 0
        else:         method = 1
    if rc is None:
        if method < 2:
            print("Methods 0 and 1 DO need rc!", file = sys.stderr)
            return
        else:
            rc = 0.0 # Set to dummy value


    pname_base = '/tmp/lammpstools_neighborize_pipe_'
    pname = pname_base + str(os.getpid())

    os.mkfifo(pname)
    neighs = []
    
    try:
        lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
        pname_buffer = create_string_buffer( pname.encode('ascii') )

        def start_neigh():
            # Make the child neighborize while the main thread reads from the pipe
            lammpstools.neighborize(void_ptr(b.x), c_longlong(b.meta.N),
                                    void_ptr(b.ids), void_ptr(b.types),
                                    c_double(rc), c_longlong(b.meta.domain.periodic),
                                    void_ptr(b.meta.domain.xlo),
                                    void_ptr(b.meta.domain.xhi), c_longlong(dims),
                                    c_longlong(method), pname_buffer,
                                    c_longlong(itype), c_longlong(jtype))

        p = multiprocessing.Process(target = start_neigh)
        p.start()

        # Read in the file and store in neighs
        current_l = []
        int_size = 8
        int_fmt  = 'q'
        ints_read = 0
        i = 0

        with open(pname,"rb") as fifo:
            byte = b'        '
            while True:
                byte = fifo.read(int_size)
                
                if byte == '' or byte == b'':
                    # Out of data apparently.
                    break;

                x = struct.unpack(int_fmt,byte)[0]
                ints_read += 1
            
                if x > 0:
                    current_l.append(x)
                else:
                    neighs.append(current_l)
                    current_l = []
        p.join()

    finally:
        if not quiet:
            print("Received %d ints (%d b.meta.tes) through named pipe %s" \
                  % (ints_read, int_size*ints_read, pname), file = sys.stderr)

        if os.path.exists( pname ):
            os.unlink( pname )
            # time.sleep(0.1)
    return neighs



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



## Attempts to identify clusters, based on a threshold criterion.
# 
#  \param neighs  Neighbor list to identify clusters in
#  \param ids     List of ids present in neighs
#  \param thresh  Cluster size threshold, i.e., only
#                 clusters of size >= thresh are returned.
def find_clusters( neighs, ids, thresh = 1 ):

    max_id = max(ids)

    ids_out = np.zeros(max_id+1,dtype=int)
    
    clusters = []
    imap = make_id_map(ids)

    if( len(neighs) > 400 ):
        print("Finding clusters, might take a while...", file = sys.stderr, end = "")
    for i in ids:
        if ids_out[i]: continue

        clust_i = [i]
        
        for j in clust_i:
            for n in neighs:
                if j in n:
                    for k in n:
                        if not ids_out[k]:
                            clust_i.append( k )
                            ids_out[k] = 1

        # Filter on thresh:
        if( len(clust_i) >= thresh ):
            clust_i.sort()
            clusters.append( clust_i )


    # Filter uniques:
    if( len(neighs) > 400 ):
        print(" Done!", file = sys.stderr)
    clusters_f = []
    
    for c in clusters:
        # Filter each:
        c_f = []
        for idi in c:
            if not idi in c_f:
                c_f.append(idi)
        clusters_f.append(c_f)

    return clusters_f


## Creates a block_data object for only those particles in the original
#  block_data \p block that are also in cluster \p cluster.
# 
#  \param block   Original block to select particles from
#  \param cluster Cluster that ocntains only the particles that should
#                 be in the new \p block_data
#
def cluster_to_block_data( block, cluster ):
    n_parts = len(cluster)

    ids   = np.empty(n_parts,dtype = int)
    types = np.empty(n_parts,dtype = int)
    x     = np.empty( [n_parts, 3], dtype = np.float64 )

    im = make_id_map(block.ids)

    j = 0
    for idi in cluster:
        idx = im[idi]
        x[j][0]  = block.x[idx][0]
        x[j][1]  = block.x[idx][1]
        x[j][2]  = block.x[idx][2]
        ids[j]   = block.ids[idx]
        types[j] = block.types[idx]
        
        j += 1
    
    
    b_filt = block_data( block.meta, ids, types, x)
    b_filt.meta.N = len(ids)
    # TODO: Filter other_cols later

    return b_filt


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
    lammpstools.ensemble_generate_manifold( void_ptr(x), N, void_ptr(ids_arr),
                                            void_ptr(types), c_longlong(domain.periodic),
                                            void_ptr(domain.xlo), void_ptr(domain.xhi),
                                            c_char_p(manifold_name), c_char_p(manifold_arg_str),
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



def triangulate( b, rc, dims, method = None ):
    """ ! Triangulates given block. """

    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")

    pname_base = '/tmp/lammpstools_triangulate_pipe_'
    pname = pname_base + str(os.getpid())

    os.mkfifo(pname)
    triangles = []
    
    try:

        def start_triangulate():
            # Make the child neighborize while the main thread reads from the pipe
            lammpstools.triangulate(void_ptr(b.x), c_longlong(b.meta.N), void_ptr(b.ids),
                                    void_ptr(b.types), c_double(rc),
                                    c_longlong(b.meta.domain.periodic),
                                    void_ptr(b.meta.domain.xlo), void_ptr(b.meta.domain.xhi),
                                    c_longlong(dims), c_longlong(method), c_char_p(pname))
            
        p = Process(target = start_triangulate)
        p.start()

        # Read in the file and store in neighs
        current_l = []
        int_size = 8
        int_fmt  = 'q'
        ints_read = 0
        i = 0
        lc = 1
        with open(pname,"rb") as fifo:
            while True:
                line = fifo.readline()

                if line == '':
                    # Out of data apparently.
                    break;
                words = line.split()
                t1 = int(words[0])
                t2 = int(words[1])
                t3 = int(words[2])
                
                triangles.append( [ t1, t2, t3 ] )
                lc += 1
        p.join()

    finally:
        if os.path.exists( pname ):
            os.unlink( pname )
            # time.sleep(0.1)

    return triangles


def triangulation_area( b, triangles ):
    """ Determines the area of given triangulation. """
    im = make_id_map( b.ids )
    A  = 0.0
    it = 1
    for t in triangles:
        # Heron's formula:
        x1 = b.x[t[0]]
        x2 = b.x[t[1]]
        x3 = b.x[t[2]]
        aa = np.linalg.norm( x1 - x2 )
        bb = np.linalg.norm( x1 - x3 )
        cc = np.linalg.norm( x2 - x3 )
        s = (aa + bb + cc) / 2.0
        aaa = math.sqrt( s*(s-aa)*(s-bb)*(s-cc) )
        A  += aaa
        it += 1
    return A

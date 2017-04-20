"""!
\file neighborize.py
\module lammpstools.py

Contains some routines relating to nearest neighbours and the like.
\inpackage lammpstools
"""

import os, multiprocessing, struct
from ctypes import *

from lammpstools.typecasts import *
from lammpstools import util

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

    rdf    = np.zeros(nbins, dtype=np.float64)
    coords = np.zeros(nbins, dtype=np.float64)
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

    pts = np.zeros(nbins,dtype=np.float64)
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
    adf    = np.zeros(nbins, dtype=np.float64)
    coords = np.zeros(nbins, dtype=np.float64)
    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
    if method is None:
        print("Please set method explicitly.", file = sys.stderr)
        return None
    
    lammpstools.compute_adf( void_ptr(b.x), c_longlong(b.meta.N), void_ptr(b.ids),
                             void_ptr(b.types), c_longlong(nbins),
                             c_longlong(itype), c_longlong(jtype),
                             c_double(R), c_longlong(method), void_ptr(adf),
                             void_ptr(coords) )
    
    pts = np.zeros(nbins,dtype=np.float64)
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
                    break

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

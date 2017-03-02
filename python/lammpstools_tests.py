"""! @package lammpstools_tests
This package contains some functions that test the lammpstools interface


\ingroup lammpstools
"""

import lammpstools
from dumpreader import *
from typecasts  import *
from b.meta.tterprint import printf
from ctypes import *
import matplotlib.pyplot as plt
import sys
import random
import math


def test_cluster():
    """! Tests the cluster finder on a fake system and a real dump file. """
    
    print "Testing cluster finder on fake system:"
    
    ids = [1,2,3,6,7,9,4,5,8,10]
    neighs = [[1, 2],
              [2, 3],
              [3, 4],
              [4 ,5],
              [5, 6],
              [6, 7],
              [7, 8],
              [8, 9],
              [9, 10] ]
    print "Neigh lists:"
    for n in neighs:
        print "[%s]" % ", ".join(map(str,n))
    
    lammpstools.find_clusters(neighs,ids,1)

    printf("\n\n")
    print "Testing cluster finder on real data:"
    # Test on real data:
    fname = "dump.cluster.gz"
    r = dumpreader(fname);
    b = r.getblock()
    if b:
        neighs = neighborize( b, 1.45, 7, 3 )
        clusters = find_clusters( b, neighs, 5 )


def test_rdf():
    """! Tests the RDF computations on real data """

    fname = "dump.rdf.gz"
    r = dumpreader(fname)
    for i in range(0,10):
        b = r.getblock()
    pts, rdf, coord  = lammpstools.compute_rdf( b, 0.0, 8.0, 101, 0, 0, 3 )
    b.meta.domain.periodic = 0
    pts, rdf2, coord = lammpstools.compute_rdf( b, 0.0, 8.0, 101, 0, 0, 3 )
    pid = os.fork()
    if pid == 0:
        plt.plot( pts, rdf, pts, rdf2 )
        plt.show()
        sys.exit(0)
    

def test_neighborize():
    """ !Tests the neighborization on real data """

    
    fname = "dump.gz"
    r = dumpreader(fname)
    for i in range(0,10):
        b = r.getblock()
    neighs = lammpstools.neighborize( b, 1.25, 3, 2 )
    print "Got neigh list of size %d" % len(neighs)
    
    #print "Neighs of atoms idx/id:"
    #for n in neighs:
    #    print "[%s]" % ", ".join(map(str,n))

    #printf("Last two lines should be:\n");
    #printf("7998/4153 ( or 4153 ): [3921, 3638, 7773]\n")
    #printf("7999/3962 ( or 3962 ): [3638, 7999, 7741, 4011]\n")

    #print "Neighs of atom with id 352 should be [27, 34, 7826, 393]"

def test_neighborize_2d_delaunay():
    """ !Tests 2D delaunay triangulation, b.meta.th periodic and not. """
    N = 8
    xx = [ [  1.1,  1.0,  0.0 ],
           [  4.5,  2.7,  0.0 ],
           [  4.6,  5.0,  0.0 ],
           [  2.4,  6.2,  0.0 ],
           [  3.4,  8.9,  0.0 ],
           [  8.0,  7.5,  0.0 ],
           [  7.0,  4.0,  0.0 ],
           [  8.8,  1.1,  0.0 ] ]
    x = np.array( xx, dtype=float )
        
    
    ids   = np.array( [1,2,3,4,5,6,7,8,9,10] )
    types = np.array( [1,1,1,1,1,1,1,1,1,1] )

    xlo = np.array( [ 0.0,  0.0,  0.0] )
    xhi = np.array( [10.0, 10.0, 10.0] )
    periodic = 0
    domain = domain_data(xlo, xhi, periodic)
    b = block_data(0, N, domain, ids, types, x)
    print "xlo: [%s]" % ", ".join(map(str, b.meta.domain.xlo) )
    print "xhi: [%s]" % ", ".join(map(str, b.meta.domain.xhi) )

    lammpstools.neighborize( b, None, 2, 2 )
    b.meta.domain.periodic = 7
    lammpstools.neighborize( b, None, 2, 2 )

def test_convex_hull():
    """ !Tests the convex hull on a set of points on a sphere """
    d = dumpreader("icosahedral.dump")
    b = d.getblock()
    neighs = lammpstools.neighborize( b, None, 3, 3 )

    print "Neigh lists:"
    for n in neighs:
        print "%d: [%s]" % (n[0], ", ".join(map(str,n[1:])))
    nns = np.zeros(b.meta.N)
    i = 0
    for n in neighs:
        # On the Python side the first id is of the particle itself, so -1:
        nns[i] = len(n)-1
        i += 1
    nns_col = dump_col("nns", nns)
    append_col(b, nns_col)
    block_to_dump(b, "icosahedral_nns.dump", overwrite=True)
    block_to_dump(b, "icosahedral_nns.gz",   overwrite=True)

    


def test_interface():
    """ !Tests if the types are passed correctly """

    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")

    fname = "dump.rdf.gz"
    r = dumpreader(fname)
    for i in range(0,10):
        b = r.getblock()
        lammpstools.test_types( void_ptr(b.x), c_longlong(b.meta.N),
                                void_ptr(b.ids), void_ptr(b.types),
                                c_longlong( b.meta.domain.periodic ), void_ptr(b.meta.domain.xlo) )

    lammpstools.test_modification( void_ptr(b.x), c_longlong(b.meta.N), 3 )
    
    nbins = 25
    rdf    = np.zeros(nbins, dtype=float);
    coords = np.zeros(nbins, dtype=float);
    
    lammpstools.test_modification( void_ptr(rdf), c_longlong(nbins), 1 )

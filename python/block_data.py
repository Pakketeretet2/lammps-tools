"""! @package block_data
This package contains definitions for block_data, which is atom info per time step

\ingroup lammpstools
"""

import numpy as np
import gzip
import os
import sys
import copy
from ctypes import *


class domain_data:
    """! A struct containing info about the domain of the dump file"""
    def __init__(self,xlo,xhi,periodic,box_line):
        self.xlo = xlo
        self.xhi = xhi
        self.periodic = periodic
        self.box_line = box_line

class block_meta:
    """! Contains the block "meta"-data, like number of atoms, domain, etc. """
    def __init__(self,t,N,domain):
        self.t = t
        self.N = N
        self.domain = domain
        self.atom_style = "atomic"

    
class block_data:
    """! Just a 'struct' that contains some basic data."""
    
    def __init__(self,meta,ids,types,x, mol = None):
        """! Initialises the block data from given block meta, 
        ids, types and positions."""
        
        self.meta = meta
        self.ids = ids
        self.types = types
        self.x = x
        if mol is None:
            self.mol = None # np.zeros( meta.N, dtype = int )
        else:
            self.mol = mol
            self.meta.atom_style = "molecular"
            
        self.other_cols = []

    
    def add_particle( self, xi, idi, typei, moli = None, idx = None ):
        """! Adds given particle to the block_data and
        updates all members accordingly.
        Do not use this to append a lot of particles as that is slow.
        Use \f merge for that.
        """
        if self.meta.atom_style == "molecular" and (moli is None and not self.mol is None):
            print("For atom_styles that have mol, pass a mol id ",
                  "as well (can be 0)!")
        # Some checks:
        if idi in self.ids:
            print("id ", idi, " is already present in ids!")
            return None

        
        self.meta.N += 1
        
        # Resize all arrays:
        grow_arrays_( self, self.meta.N )
        

        if idx is None:
            idx = self.meta.N-1
            self.ids[idx]   = idi
            self.types[idx] = typei
            self.x[idx]     = xi
            if self.mol is not None:
                self.mol[idx] = moli
        else:
            # Gotta swap...
            jdx = self.meta.N-1
            self.ids[jdx] = self.ids[idx]
            self.types[jdx] = self.types[idx]
            self.x[jdx] = self.x[idx]

            self.ids[idx]   = idi
            self.types[idx] = typei
            self.x[idx]     = xi

            if self.mol is not None:
                self.mol[jdx] = self.mol[idx]
                self.mol[idx] = moli


    def merge( self, other_block, error_on_double_id = True ):
        """! Merges \p other_block into this one. """
        if error_on_double_id:
            self.merge_keep_ids_( other_block )
        else:
            self.merge_fresh_ids_( other_block )

        # After merging all particle data, merge the metadata.
        newxlo = np.zeros( 3, dtype = float )
        newxhi = np.zeros( 3, dtype = float )

        my_dom = self.meta.domain
        o_dom  = other_block.meta.domain

        for i in range(0,2):
            newxlo[i] = min( my_dom.xlo[i], o_dom.xlo[i] )
            newxhi[i] = max( my_dom.xhi[i], o_dom.xhi[i] )

        self.meta.domain.xlo = newxlo
        self.meta.domain.xhi = newxhi

        if (my_dom.periodic != o_dom.periodic or
            my_dom.box_line != o_dom.box_line):
           print("Some domain settings were not consistent between blocks. ",
                 "Throwing away info of other block!", file = sys.stderr )

    
    def merge_fresh_ids_( self, other_block ):
        """! Merges \p other_block but sets the ids of its particles
        to 'fresh' values starting from the largest id in self + 1. """
        old_max = self.meta.N
        new_id = np.max( self.ids ) + 1
        new_N = self.meta.N + other_block.meta.N
        self.grow_arrays_( new_N )
        
        for j in range( old_max, new_N ):
            i = j-old_max
            self.ids[j]   = new_id
            self.x[j]     = other_block.x[i]
            self.types[j] = other_block.types[i]
            if self.meta.atom_style == "molecular":
                self.mol[j] = other_block.mol[i]
            
            new_id += 1
            

    def merge_keep_ids_( self, other_block ):
        """! Merges \p other_block and keeps its ids. Errors if double 
        id is encountered. """
        old_max = self.meta.N
        new_N = self.meta.N + other_block.meta.N
        
        self.grow_arrays_( new_N )
        for j in range( old_max, new_N ):
            i = j-old_max
            if other_block.ids[i] in self.ids[0:old_max]:
                raise RuntimeError("Double ID encountered!")
            
            self.ids[j]   = other_block.ids[i]
            self.x[j]     = other_block.x[i]
            self.types[j] = other_block.types[i]
            if self.meta.atom_style == "molecular":
                self.mol[j] = other_block.mol[i]
            
            new_id += 1
        
            
    def grow_arrays_( self, new_N ):
        """! Resizes the arrays and particle number. """
        self.ids   = np.resize(self.ids,   new_N)
        self.types = np.resize(self.types, new_N)
        self.x     = np.resize(self.x,     [new_N,3])
        
        if self.meta.atom_style == "molecular":
            self.mol = np.resize(self.mol, self.meta.N)
        self.meta.N = new_N


class dump_col:
    """! For if you have/want additional columns in your dump file."""
    def __init__(self,header,data):
        self.header = header
        self.data   = data
        self.N      = len(data)

    


def append_col(block, col):
    """! Appends column to other_cols of block"""
    if( block.meta.N != col.N):
        print("Size mismatch between column and block info!")
        return
    block.other_cols.append(col)
    return block


def read_data(fname, atom_style = None):
    """! Reads a data file, converts it to block_data and returns that. """
    with open(fname) as fin:
        
        N = 0
        t = 0
        xlo = np.array([0,0,0],dtype=np.float64)
        xhi = np.array([0,0,0],dtype=np.float64)
        dom = domain_data(xlo, xhi, 0, "")
        mol = None

        skip_two = 0
        if atom_style is None:
            print("Assuming atom_style is atomic.", file=sys.stderr)
            atom_style = "atomic"
        
        if atom_style == "atomic":
            nwords_expect = 5

        elif atom_style == "molecular":
            nwords_expect = 6

        else:
            print("Atom style ", atom_style, " not supported.")
            return None
            
        while True:
            line = fin.readline()
            if skip_two < 2:
                skip_two += 1
                continue
            
            if len(line) == 0:
                # At EOF.
                break
        
            line = line.rstrip()
            if not len(line) == 0:
                # Split line in words:
                words = line.split()
                if( len(words) > 1 ):
                    if words[1] == "atoms":
                        N = int(words[0])

                        ids = np.zeros(N,dtype=int)
                        types = np.zeros(N,dtype=int)
                        x     = np.zeros([N,3], dtype=float)

                        if atom_style == "molecular":
                            mol = np.zeros(N, dtype = int)
        
                if( len(words) > 3 ):
                    if( words[2] == "xlo" and words[3] == "xhi" ):
                        dom.xlo[0] = float(words[0])
                        dom.xhi[0] = float(words[1])
                    if( words[2] == "ylo" and words[3] == "yhi" ):
                        dom.xlo[1] = float(words[0])
                        dom.xhi[1] = float(words[1])
                    if( words[2] == "zlo" and words[3] == "zhi" ):
                        dom.xlo[2] = float(words[0])
                        dom.xhi[2] = float(words[1])
                if( words[0] == "Atoms" ):
                    # Read the entire Atoms block.
                    line = fin.readline()
                    jjj = 0
                    
                    while jjj < N:
                        line = fin.readline()
                        line = line.rstrip()

                        words = line.split()
                        if len(words) != nwords_expect:
                            print("Incorrect formatting of data file for ",
                                  "atom_style ", atom_style)
                            return None
                        if atom_style == "atomic":
                            idi = int(words[0])
                            ity = int(words[1])
                            xi  = float(words[2])
                            yi  = float(words[3])
                            zi  = float(words[4])
                        elif atom_style == "molecular":
                            idi = int(words[0])
                            im  = int(words[1])
                            ity = int(words[2])
                            xi  = float(words[3])
                            yi  = float(words[4])
                            zi  = float(words[5])

                        ids[jjj] = idi;
                        types[jjj] = ity;
                        x[jjj,0] = xi
                        x[jjj,1] = yi
                        x[jjj,2] = zi
                        if mol is not None:
                            mol[jjj] = im;
                        
                        
                        jjj += 1

    bm = block_meta(t,N,dom)
    b = block_data(bm,ids=ids,types=types,x=x, mol = mol)
    return b



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




def block_to_dump_write(b, fname, write_mode = "w", file_mode = "plain"):
    if file_mode == "plain":
        dump_file = open(fname, write_mode)
    elif file_mode == "gz":
        dump_file = gzip.GzipFile(fname,write_mode)
    else:
        print("File mode %s does not exist!" % file_mode, file = sys.stderr)
        sys.exit(-1)

    print("ITEM: TIMESTEP", file = dump_file)
    print(b.meta.t, file = dump_file)
    print("ITEM: NUMBER OF ATOMS", file = dump_file)
    print(b.meta.N, file = dump_file)
    print("ITEM: BOX BOUNDS", b.meta.domain.box_line, file = dump_file)
    print("%f  %f" % ( b.meta.domain.xlo[0], b.meta.domain.xhi[0] ), file = dump_file)
    print("%f  %f" % ( b.meta.domain.xlo[1], b.meta.domain.xhi[1] ), file = dump_file)
    print("%f  %f" % ( b.meta.domain.xlo[2], b.meta.domain.xhi[2] ), file = dump_file)

    if b.meta.atom_style == "atomic":
        atom_string = "ITEM: ATOMS id type x y z"
    elif b.meta.atom_style == "molecular":
        atom_string = "ITEM: ATOMS id mol type x y z"
    
    for c in b.other_cols:
        if c.N != b.meta.N:
            print("Error! Column size does not match rest of block data!")
        atom_string += " " + c.header
    print(atom_string, file = dump_file)
    for i in range(0,b.meta.N):
        #atom_line = "%d %d %f %f %f" % (b.ids[i], b.types[i],
        #                                b.x[i][0], b.x[i][1], b.x[i][2])
        if b.meta.atom_style == "atomic":
            print(b.ids[i], b.types[i], b.x[i][0], b.x[i][1], b.x[i][2], end = "", file = dump_file)
        elif b.meta.atom_style == "molecular":
            print(b.ids[i], b.mol[i], b.types[i], b.x[i][0], b.x[i][1], b.x[i][2], end = "", file = dump_file)


        for c in b.other_cols:
            print(c.data[i], file = dump_file, end = "")
        print("", file = dump_file)
    

    
def block_to_dump(b,fname, write_mode):
    """ !Writes given block data to a dump file. """
    
    if fname[len(fname)-2:] == "gz":
        file_mode = "gz"
    else:
        file_mode = "plain"

    block_to_dump_write(b, fname, write_mode, file_mode)

    
    
def block_to_data(b,fname, overwrite = False):
    """ !Writes given block data to a dump file. """

    write_mode = "w"    

    
    # Determine some important params for the data file:
    ntypes = np.max( b.types )
    natoms = b.meta.N
    
    f = open(fname, write_mode)

    print("LAMMPS data file from lammpstools.py.", file = f )
    print("", file = f )
    print("%d atoms"      % natoms, file = f )
    print("%d atom types" % ntypes, file = f )
    print("", file = f )
    print("%f %f xlo xhi" % (b.meta.domain.xlo[0], b.meta.domain.xhi[0]), file = f )
    print("%f %f ylo yhi" % (b.meta.domain.xlo[1], b.meta.domain.xhi[1]), file = f )
    print("%f %f zlo zhi" % (b.meta.domain.xlo[2], b.meta.domain.xhi[2]), file = f )
    print("", file = f )
    print("Atoms", file = f )
    print("", file = f )

    if( b.meta.atom_style == "molecular" ):
        for i in range(0,b.meta.N):
            print("%d %d %d %f %f %f" % (b.ids[i], b.mol[i], b.types[i], b.x[i][0], b.x[i][1], b.x[i][2]), file = f )
    else:
        for i in range(0,b.meta.N):

            print("%d %d %f %f %f" % (b.ids[i], b.types[i], b.x[i,0], b.x[i,1], b.x[i,2]), file = f )
    


def block_to_xyz_write(b, fname, write_mode = "w"):
    """! Dumps block to xyz file. """

    fp = open(fname,write_mode)

    print(b.meta.N, file = fp)
    print("Atoms. Timestep ", b.meta.t, file = fp )
    
    # xyz needs to be sorted along ids...

    permutation = np.argsort( b.ids )

    for i in permutation:
        print(b.types[i], b.x[i][0], b.x[i][1], b.x[i][2], file = fp, end = "")

        for c in b.other_cols:
            print(c.data[i], file = fp, end = "")
        print("", file = fp)


def copy_meta_new_block( b_old, Xnew, id_new = None, type_new = None, mol_new = None ):
    """! Makes new block_data object for new positions with old meta. """

    bm   = copy.deepcopy(b_old.meta)
    bm.N = Xnew.shape[0]
    X    = np.zeros( [ bm.N, 3 ], dtype = float )

    ids   = np.ones( bm.N, dtype = int )
    types = np.ones( bm.N, dtype = int )
    mol   = np.ones( bm.N, dtype = int )

    
    for i in range(0,bm.N):
        X[i,:]   = Xnew[i,:]
        
        if id_new is None:
            ids[i] = i+1
        else:
            ids[i] = id_new[i]

        if type_new is None:
            types[i] = 1
        else:
            types[i] = type_new[i]

        if mol_new is None:
            mol[i]   = ids[i]
        else:
            mol[i] = mol_new[i]
    
    bb = block_data( bm, ids, types, X, mol )
    return bb


def block_data_from_foreign( X, ids, types, mol, periodic, xlo, xhi,
                             dims, tstep, boxline ):
    """! Constructs a block_data object from given arrays. """

    N    = len(X)
    dom  = block_data.domain_data( xlo, xhi, periodic, boxline )
    meta = block_data.block_meta( 0, N, dom )
    b    = block_data.block_data( meta, ids, types, X, mol )
    return b

def blocks_from_xyz( fname, pad = None, xxlo = None, xxhi = None ):
    """! Constructs a list of block_datas from a given xyz file. """

    blocks = []
    if pad is None:
        pad = 2

    tstep = 0
    with open(fname,"r") as f:
        for l in f:
            l = l.rstrip()
            N = int(l)
            X     = np.zeros( [N,3] )
            ids   = np.array( [ i for i in range(1,N+1) ], dtype = int )
            types = np.ones( N, dtype = int )
            l = f.readline()
            xlo = np.zeros(3)
            xhi = np.zeros(3)
        
            for i in range(0,N):
                l = f.readline()
                l = l.rstrip()
                w = l.split()

                X[i][0] = float(w[1])
                X[i][1] = float(w[2])
                X[i][2] = float(w[3])

                for k in range(0,3):
                    if X[i][k] < xlo[k]:
                        xlo[k] = X[i][k]
                    elif X[i][k] > xhi[k]:
                        xhi[k] = X[i][k]

            xlo -= pad
            xhi += pad
            if not xxlo is None:
                xlo = xxlo
            if not xxhi is None:
                xhi = xxhi
            dom  = domain_data( xlo, xhi, 0, "ITEM: BOX BOUNDS ff ff ff" )
            meta = block_meta( tstep, N, dom )
            b    = block_data( meta, ids, types, X, None )
            blocks.append(b)
            tstep += 1
    return blocks


def new_block_data_cpp( b ):
    """! Creates a C++-style block_data struct and returns the handle to it. """
    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
    bh = lammpstools.new_block_data()

    if b.mol is None:
        molref = None
    else:
        molref = b.mol.ctypes.data_as(POINTER(c_longlong))
        
    boxline_buff = create_string_buffer( b.meta.domain.box_line.encode('ascii') )

    lammpstools.set_block_data( bh, b.meta.N, b.meta.t,
                                b.x.ctypes.data_as(POINTER(c_double)),
                                b.ids.ctypes.data_as(POINTER(c_longlong)),
                                b.types.ctypes.data_as(POINTER(c_longlong)),
                                molref,
                                b.meta.domain.xlo.ctypes.data_as(POINTER(c_double)),
                                b.meta.domain.xhi.ctypes.data_as(POINTER(c_double)),
                                b.meta.domain.periodic, boxline_buff )
    
    
    return bh

        
def free_block_data_cpp( bh ):
    """! Deletes a C++-style block_data struct. """
    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
    lammpstools.free_block_data( bh )
    

def write_block_data( b, file_name, file_format, data_format ):
    """! Writes given block_data to specified file and format.
    
    This is just a wrapper around write_block_data_cpp, see that
    documentation for the parameters. """
    bh = new_block_data_cpp( b )
    try:
        write_block_data_cpp( bh, file_name, file_format, data_format )
    finally:
        free_block_data_cpp( bh )
        
                      
def write_block_data_cpp( bh, file_name, file_format, data_format ):
    """! Writes given block_data to specified file and format.
    @param b            Block data to write.
    @param fname        File name to write to.
    @param file_format  File format to write. Currently supported:
                        'PLAIN', 'GZIP', 'BIN'
    @param data_format  Data format to write. Currently supported:
                        'LAMMPS', 'HOOMD'
    """
    
    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
    fformat = file_format.encode( 'ascii' )
    dformat = data_format.encode( 'ascii' )
    fname   = file_name.encode( 'ascii' )
    
    lammpstools.write_block_to_file( bh, fname, fformat, dformat )

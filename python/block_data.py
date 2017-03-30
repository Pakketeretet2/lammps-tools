"""! @package block_data
This package contains definitions for block_data, which is atom info per time step

\ingroup lammpstools
"""

import numpy as np
import gzip
import os
import sys
import copy


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
        self.meta = meta
        
        self.ids = ids
        self.types = types
        self.x = x
        if mol is None:
            self.mol = None # np.zeros( meta.N, dtype = int )
        else:
            self.mol = mol
            
        self.other_cols = []


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


def read_data(fname):
    """! Reads a data file, converts it to block_data and returns that. """
    with open(fname) as fin:
        
        N = 0
        t = 0
        xlo = np.array([0,0,0],dtype=np.float64)
        xhi = np.array([0,0,0],dtype=np.float64)
        dom = domain_data(xlo, xhi, 0, "")

        skip_two = 0
        
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
                        idi = int(words[0])
                        ity = int(words[1])
                        xi  = float(words[2])
                        yi  = float(words[3])
                        zi  = float(words[4])

                        ids[jjj] = idi;
                        types[jjj] = ity;
                        x[jjj,0] = xi
                        x[jjj,1] = yi
                        x[jjj,2] = zi
                        
                        jjj += 1

    bm = block_meta(t,N,dom)
    b = block_data(bm,ids=ids,types=types,x=x)
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
    elif b.meta.atom_style == "molecule":
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
        elif b.meta.atom_style == "molecule":
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
            print("%d %d %f %f %f" % (b.ids[i], b.types[i], b.x[i][0], b.x[i][1], b.x[i][2]), file = f )
    


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


"""! \file block_data.py

This file contains definitions for block_data, which is atom info per time step

\ingroup lammpstools
"""

import numpy as np
import gzip
import os
import sys
import copy
from ctypes import *

## Contains some basic info about a simulation domain.
#
class domain_data:
    ## Constructor
    #  
    #  \param xlo        Lower bounds of the box
    #  \param xhi        Upper bounds of the box
    #  \param periodic   Int that contains the periodicity flags. Its value is
    #                    0 + (x periodic?)*1 + (y periodic?)*2 + (z periodic?)*4
    #  \param box_line   A string that describes the box in the LAMMPS format.
    def __init__(self,xlo,xhi,periodic,box_line):
        self.xlo = xlo             ## Lower bounds of box
        self.xhi = xhi             ## Upper bounds of box
        self.periodic = periodic   ## Int that contains the periodicity flags
        self.box_line = box_line   ## LAMMPS-style string describing the box

    ## Calculates the distance for this geometry.
    #  
    #  \param xi    Position of particle i
    #  \param xj    Position of particle j
    #
    #  \returns     A tuple with the Euclidian distance and the distance vector
    def distance(self, xi, xj):
        Lx = self.xhi[0] - self.xlo[0]
        Ly = self.xhi[1] - self.xlo[1]
        Lz = self.xhi[2] - self.xlo[2]
        r  = xi - xj
        if self.periodic & 1:
            if r[0] >  0.5*Lx:
                r[0] -= Lx
            elif r[0] <= -0.5*Lx:
                r[0] += Lx
    
        if self.periodic & 2:
            if r[1] >  0.5*Ly:
                r[1] -= Ly
            elif r[1] <= -0.5*Ly:
                r[1] += Ly
    
        if self.periodic & 4:
            if r[2] >  0.5*Lz:
                r[2] -= Lz
            elif r[2] <= -0.5*Lz:
                r[2] += Lz
                
        dist = np.linalg.norm(r)
        
        return dist, r

## Contains some metadata about the time step.
#
#  Members:
#  \param  t           Time step of the block_data
#  \param  N           Number of particles
#  \param  domain      Domain data for the block_data
#  \param  atom_style  Atom style
# 
class block_meta:
    ## Constructor
    #
    #  \param t        The time step
    #  \param N        The number of particles
    #  \param domain   The domain data
    def __init__(self,t,N,domain):
        self.t = t                  ## 
        self.N = N
        self.domain = domain
        self.atom_style = "atomic"

## Block data struct:
#
#  Members:
#  \param meta   Block metadata    
#  \param ids    Particle ids      
#  \param types  Particle types    
#  \param x      Particle positions
#  \param mol    Molecule ids (None if meta.atom_style doesn't support mols)
#
class block_data:
    ## Constructor that takes arrays and converts them into a POD struct
    #  \param meta   Block metadata    
    #  \param ids    Particle ids      
    #  \param types  Particle types    
    #  \param x      Particle positions
    #  \param mol    Molecule ids
    def __init__(self,meta,ids,types, x, mol = None):

        # Check lengths:
        if not all(meta.N == length for length in
                   [ len(ids), len(types), len(x) ] ):
            raise RuntimeError("Length mismatch in arrays passed to block_data!")
        if (not mol is None) and (len(mol) != meta.N):
            raise RuntimeError("Length mismatch in arrays passed to block_data!")
        
        self.meta  = meta
        self.ids   = ids
        self.types = types
        self.x     = x
        if mol is None:
            self.mol = None
        else:
            self.mol = mol
            self.meta.atom_style = "molecular"
        
        self.other_cols = []


## A class for additional dump columns.
#
#  Members:
#  \param header  String describing the column
#  \param data    Data of the column
#  \param N       Number of atoms/data entries
class dump_col:
    ## Constructor.
    #
    #  \param header  String describing the column
    #  \param data    Data of the column
    def __init__(self,header,data):
        self.header = header
        self.data   = data
        self.N      = len(data)

        
## Adds given particle to \p block and returns a new block_data.
#  Do not use this to append a lot of particles as that is slow.
#  Use \f merge for that.
#
#  \param block  The block to append the particle to.
#  \param xi     The particle position
#  \param idi    The particle ID
#  \param typei  The particle type
#  \param moli   The molecule the particle belongs to
#    
def add_particle( block, xi, idi, typei, moli = None ):
    if ( block.meta.atom_style == "molecular" and
         (moli is None and not block.mol is None) ):
        print("For atom_styles that have mol, pass a mol id ",
              "as well (can be 0)!")
        raise RuntimeError("Inconsistent mol id passed!")
    # Some checks:
    if idi in block.ids:
        raise RuntimeError("id ", idi, " is already present in ids!")

    have_mol = (not block.mol is None)
    N     = block.meta.N
    x     = np.zeros( [N+1,3], dtype = float )
    ids   = np.zeros( N+1, dtype = int )
    types = np.zeros( N+1, dtype = int )
    if have_mol:
        mol = np.zeros( N+1, dtype = int )
    
    for i in range(0,N):
        x[i] = block.x[i]
        ids[i] = block.ids[i]
        types[i] = block.types[i]
        if have_mol:
            mol[i] = block.mol[i]

    x[N] = xi
    ids[N] = idi
    types[N] = typei
    if have_mol: 
        mol[N] = moli

    meta = copy.deepcopy(block.meta) # I hate Python.
    #meta = block.meta
    meta.N = N + 1

    if meta.N == block.meta.N:
        raise RuntimeError("Object mutated! I hate Python!")
    if have_mol:
        return block_data( meta, ids, types, x, mol )
    else:
        return block_data( meta, ids, types, x, None )
    

## Merges two block_data objects into a fresh one.
#
#  \param block         The first block
#  \param other_block   The second block
#  
def merge( block, other_block ):
    # After merging all particle data, merge the metadata.
    newxlo = np.zeros( 3, dtype = float )
    newxhi = np.zeros( 3, dtype = float )

    my_dom = self.meta.domain
    o_dom  = other_block.meta.domain

    for i in range(0,2):
        newxlo[i] = min( my_dom.xlo[i], o_dom.xlo[i] )
        newxhi[i] = max( my_dom.xhi[i], o_dom.xhi[i] )

    if (my_dom.periodic != o_dom.periodic or
        my_dom.box_line != o_dom.box_line):
        raise RuntimeError("Some domain settings were not consistent between blocks!")

    M = block.meta.N + other_block.meta.N

    
    meta   = copy.deepcopy( block.meta )
    meta.N = M
    ids    = np.zeros( M, dtype = int )
    types  = np.zeros( M, dtype = int )
    x      = np.zeros( [M,3], dtype = float )
    
    if block.meta.atom_style != other_block.meta.atom_style:
        raise RuntimeError("Atom styles inconsistent!")
    if block.meta.atom_style == "atomic":
        mols = None
    else:
        mols   = np.zeros( M, dtype = int )

    for i in range(0,block.meta.N):
        ids[i]   = block.ids[i]
        types[i] = block.types[i]
        x[i]     = block.x[i]
        if not mols is None:
            mols[i] = block.mols[i]

    offset = block.meta.N
    for i in range(0,other_block.meta.N):
        j = i + offset
        ids[j]   = other_block.ids[i]
        types[j] = other_block.types[i]
        x[j]     = other_block.x[i]
        if not mols is None:
            mols[j] = other_block.mols[i]
    

    return block_data( meta, ids, types, x, mols )
    
    



def append_col(block, col):
    """! Appends column to other_cols of block"""
    if( block.meta.N != col.N):
        raise RuntimeError("Size mismatch between column and block info!")
    block.other_cols.append(col)
    return block



def block_data_from_foreign( X, ids, types, mol, periodic, xlo, xhi,
                             dims, tstep, boxline ):
    """! Constructs a block_data object from given arrays. """

    N    = len(X)
    dom  = domain_data( xlo, xhi, periodic, boxline )
    meta = block_meta( 0, N, dom )
    b    = block_data( meta, ids, types, X, mol )
    return b


## Writes given block_data to file.
#
#  \param b      block to write.
#  \param fname  File name to write to.
def block_to_data( b, fname ):
    bh = new_block_data_cpp( b )
    write_block_data( b, fname, "PLAIN", "LAMMPS_DATA" )
    
    free_block_data_cpp( bh )


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



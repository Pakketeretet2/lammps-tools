"""! \package lammpstools
\file dumpreader_cpp.py
This package contains functions for extracting dump files by delegating
it to the C++ lib. This leads to better performance and the code is cleaner. XD

\ingroup lammpstools
"""

import numpy as np
import gzip
import os
import sys
import copy

from ctypes import *
import ctypes

from lammpstools.block_data import *



class dumpreader_cpp:
    """
    Interface to the C++ dump reader. This is the preferred method of
    reading dump files from Python due to its enhanced performance and better
    flexibility.
    """
    
    def __init__(self, fname, dformat = None, fformat = None):
        ## This is the constructor for the dump reader.
        #  @arg fname   Dump file name
        #  @arg dformat Dump format (0 = LAMMPS, 1 = GSD, 2 = DCD)
        #  @arg fformat File format (0 = plain text, 1 = GZipped, 3 = binary)
        #
        #  It opens a handle to an instance of a C++ dump_reader which does
        #  all of the heavy lifting. The handle is released in the destructor
        ## 
        lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")

        if fformat == 2:
            raise LogicError("A file format of 2 (Input stream) is not supported!")

        if dformat is None:
            dformat = -1
        if fformat is None:
            fformat = -1

        self.handle = lammpstools.get_dump_reader_handle( fname.encode(),
                                                          dformat, fformat )
            
        print("Opened dump reader handle @ ", hex(self.handle),
              file = sys.stderr)
        self.at_eof = False
        

    def __del__(self):
        ## This is the deconstructor.
        ## It releases handle to the C++ dump_reader.
        lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
        print("Releasing dump reader handle @ ", hex(self.handle),
              file = sys.stderr)
        lammpstools.release_dump_reader_handle( self.handle )


    def __iter__(self):
        ## Makes the dumpreader iterable, i.e., you can do for b in d:
        return self

    def __next__(self):
        ## Iterates to next block.
        if not self.at_eof:
            tmp_b = self.getblock()
            if tmp_b is None:
                raise StopIteration
            else:
                return tmp_b
        else:
            raise StopIteration

        

    def getblock(self):
        ## Returns a block_data containing the info of the next block,
        ## or None if the next block could somehow not be read.
        # print("Gonna grab new block now...", file=sys.stderr)
        lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
        domain = domain_data( np.array( [0, 0, 0] ), np.array( [0, 0, 0] ),
                              0, "ITEM: BOX BOUNDS pp pp pp" )

        N = ctypes.c_longlong(0)
        tstep = ctypes.c_longlong(0)
        xlo = np.zeros( 3, dtype = float )
        xhi = np.zeros( 3, dtype = float )
        atom_style = ctypes.c_longlong(0)
        
        periodic = ctypes.c_longlong(0)

        # Should fit "ITEM: BOX BOUNDS pp pp pp\0"
        box_line_buff = ctypes.create_string_buffer(26)

        if lammpstools.dump_reader_next_block( self.handle ):
            # print("An error happened reading the next block!",file=sys.stderr)
            self.at_eof = True
            return None

        lammpstools.dump_reader_get_block_meta( self.handle, byref(tstep), byref(N),
                                                ctypes.c_void_p(xlo.ctypes.data),
                                                ctypes.c_void_p(xhi.ctypes.data),
                                                byref(periodic),
                                                box_line_buff,
                                                byref(atom_style) )
        
        box_line = str( box_line_buff, 'ascii' )
        x     = np.empty( [N.value, 3], dtype = float )
        ids   = np.empty( N.value, dtype = int )
        types = np.empty( N.value, dtype = int )
        mol   = np.zeros( N.value, dtype = int )


        if atom_style.value == 0:
            atom_style_named = "atomic"
        elif atom_style.value == 1:
            atom_style_named = "molecular"
        else:
            raise RuntimeError("Unknown atom style ", atom_style.value,
                               " encountered!")
        
        lammpstools.dump_reader_get_block_data(
            self.handle, N, x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ids.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong)),
            types.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong)),
            mol.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong)) )
        
        
        dom  = domain_data( xlo, xhi, periodic.value, box_line )
        meta = block_meta( tstep.value, N.value, dom )
        meta.atom_style = atom_style_named

        b = block_data( meta, ids, types, x, mol )
        return b

    
    def fast_forward(self,Nblocks=1):
        ## Fast-forwards a number of blocks.
        #  @arg Nblocks  The number of blocks to skip.
        ##
        if self.at_eof:
            print("At EOF already, cannot get more blocks!", file = sys.stderr)
            return None

        lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")

        success = lammpstools.dump_reader_fast_forward( self.handle, Nblocks )
        if success:
            pass
        else:
            print("Fast-forward failed, probably encountered EOF.", file = sys.stderr)
            self.at_eof = True

    def get_all_blocks(self):
        ## Returns Returns all blocks in the data file in a list of blocks.
        ## This can be very memory-heavy!
        all_data = []
        for b in self:
            all_data.append(b)
        return all_data

        

"""! @package dumpreader_cpp
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
from block_data import *


class dumpreader_cpp:
    def __init__(self, fname ):
        lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
        self.handle = lammpstools.get_dump_reader_handle( fname.encode() )
        print("Opened dump reader handle @ ", hex(self.handle))
        self.at_eof = False

    def __del__(self):
        lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
        print("Releasing dump reader handle @ ", hex(self.handle))
        lammpstools.release_dump_reader_handle( self.handle )

    def __iter__(self):
        return self

    def __next__(self):
        if not self.at_eof:
            tmp_b = self.getblock()
            if tmp_b is None:
                raise StopIteration
            else:
                return tmp_b
        else:
            raise StopIteration

        

    def getblock(self):
        lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")
        domain = domain_data( np.array( [0, 0, 0] ), np.array( [0, 0, 0] ),
                              0, "pp pp pp" )

        N = ctypes.c_longlong(0)
        tstep = ctypes.c_longlong(0)
        xlo = np.zeros( 3, dtype = float )
        xhi = np.zeros( 3, dtype = float )
        atom_style = ctypes.c_longlong(0)
        
        periodic = ctypes.c_longlong(0)

        # Should fit "ITEM: BOX BOUNDS pp pp pp\0"
        box_line_buff = ctypes.create_string_buffer(32)

        if not lammpstools.dump_reader_next_block( self.handle ):
            # print("An error happened reading the next block!",file=sys.stderr)
            self.at_eof = True
            return None


        lammpstools.dump_reader_get_block_meta( self.handle, byref(tstep), byref(N),
                                                ctypes.c_void_p(xlo.ctypes.data),
                                                ctypes.c_void_p(xhi.ctypes.data),
                                                byref(periodic),
                                                box_line_buff,
                                                byref(atom_style) )
        box_line = str( box_line_buff, 'utf-8' )
        
        x     = np.empty( [N.value, 3], dtype = float )
        ids   = np.empty( N.value, dtype = int )
        types = np.empty( N.value, dtype = int )
        mol   = np.zeros( N.value, dtype = int )

        
        lammpstools.dump_reader_get_block_data(
            self.handle, N, x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ids.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong)),
            types.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong)),
            mol.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong)) )
        
        
        dom  = domain_data( xlo, xhi, periodic.value, box_line )
        meta = block_meta( tstep.value, N.value, dom )
        b = block_data( meta, ids, types, x, mol )
        
        return b

    
    def fast_forward(self,Nblocks=1):
        """! Fast-forwards Nblocks. """
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
        """! Returns all blocks in the data file in a list of blocks.
             This can be very memory-heavy!"""
        all_data = []
        for b in self:
            all_data.append(b)
        return all_data

        

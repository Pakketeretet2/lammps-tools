import lammpstools, block_data, dumpreader
import sys, numpy as np, os, ctypes

if len(sys.argv) < 2:
    print( "Pass a dump file!", file = sys.stderr )
    sys.exit(-1)

d = dumpreader.dumpreader_cpp( sys.argv[1] )
b = d.getblock()

block_data.write_block_data( b, "test.gsd",  "BIN",   "HOOMD" )
block_data.write_block_data( b, "test.dump", "PLAIN", "LAMMPS" )

d = dumpreader.dumpreader_cpp( "test.gsd" )
b = d.getblock()


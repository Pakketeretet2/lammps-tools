import lammpstools, block_data, dumpreader
import sys, numpy as np, os, ctypes

if len(sys.argv) < 2:
    print( "Pass a dump file!", file = sys.stderr )
    sys.exit(-1)

d = dumpreader.dumpreader( sys.argv[1] )
b = d.getblock()

bh = block_data.new_block_data_cpp( b )
block_data.free_block_data_cpp( bh )




import lammpstools as lt
import sys, numpy as np, os, ctypes

def test_dump_read_write( fname ):
    d = lt.dumpreader_cpp( fname )
    b = d.getblock()
    
    lt.write_block_data( b, "test.gsd",  "BIN",   "HOOMD" )
    lt.write_block_data( b, "test.dump", "PLAIN", "LAMMPS" )
    
    d = lt.dumpreader_cpp( "test.gsd" )
    b = d.getblock()

def test_data_read_write( fname ):
    d = lt.dumpreader_cpp(fname)
    b = d.getblock()
    lt.write_block_data( b, "test.data", "PLAIN", "LAMMPS_DATA" )
    b2 = lt.read_data( "test.data" )
    print("b2.x[51] = ", b2.x[51], ", b.x[51] = ", b.x[51])

def test_data_read_two():
    b2 = lt.read_data( "coords_min_lj.data" )
    lt.block_to_data(b2, "-")

if len(sys.argv) < 2:
    print( "Pass a dump file!", file = sys.stderr )
    sys.exit(-1)

#test_dump_read_write( sys.argv[1] )
#print("Test dump read write OK")
#
#test_data_read_write( sys.argv[1] )
#print("Test data read write OK")

test_data_read_two()


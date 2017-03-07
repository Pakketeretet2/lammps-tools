import dumpreader
import sys

dump_name = "melt.dump"
d = dumpreader.dumpreader( dump_name, x_tag = "xs", y_tag = "ys", z_tag = "zs" )
b = d.getblock()

while not d.at_eof:
    print("t = ", b.meta.t)
    b = d.getblock()

d  = dumpreader.dumpreader( dump_name, x_tag = "xs", y_tag = "ys", z_tag = "zs" )
b  = d.getblock()
dumpreader.block_to_dump_write( b, "/dev/stdout", "w" )

d2 = dumpreader.dumpreader_cpp( dump_name )
b = d2.getblock()

dumpreader.block_to_dump_write( b, "/dev/stdout", "w" )

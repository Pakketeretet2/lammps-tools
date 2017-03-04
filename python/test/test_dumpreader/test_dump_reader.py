import dumpreader
import sys

dump_name = "melt.dump"
d = dumpreader.dumpreader( dump_name, x_tag = "xs", y_tag = "ys", z_tag = "zs" )
b = d.getblock()

while not d.at_eof:
    print("t = ", b.meta.t)
    b = d.getblock()


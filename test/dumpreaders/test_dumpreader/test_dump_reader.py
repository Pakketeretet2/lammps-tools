import dumpreader
import sys
import time

dump_files = [ 'melt.dump', 'melt.dump.gz', 'init.gsd' ]

for dname in dump_files:
    d = dumpreader.dumpreader_cpp(dname)
    for b in d:
        print("t = ",b.meta.t,", N = ",b.meta.N)


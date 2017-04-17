import lammpstools
import dumpreader
import numpy as np
import sys

if len(sys.argv) < 2:
    print >> sys.stderr, "Pass a dump file!"
    sys.exit(-1)


fname = sys.argv[1]
d = dumpreader.dumpreader(fname)
b = d.getblock()

# Assume particles on sphere.
R2 = 0.0;
c  = 1.0 / float(b.meta.N)
for i in range(0,b.meta.N):
    R2 += np.dot( b.x[i], b.x[i] ) * c

R = math.sqrt(R2)
print >> sys.stderr, "<R> = %f" % R

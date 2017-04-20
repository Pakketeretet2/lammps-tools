import lammpstools as lt
import numpy as np

d = lt.dumpreader_cpp( "melt.dump" )
b  = d.getblock()
b2 = lt.add_particle( b, np.array( [0.0,0.0,1.0] ), np.max(b.ids) + 1, 1, 0 )


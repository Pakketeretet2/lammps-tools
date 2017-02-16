import lammpstools
import dumpreader
import hexorder
import numpy as np
import math

fname = "traj.dump"
fname2 = "psi6.dump"


d = dumpreader.dumpreader(fname) #, x_tag = "xs", y_tag = "ys", z_tag = "zs")
b = d.getblock()

ax = np.array( [math.sqrt(2), 0*math.sqrt(2), 0.0] )
rc = 1.9
ax /= np.linalg.norm(ax)

i = 0

while not d.at_eof:
    psi_avg, psis = hexorder.compute_psi_n( b, ax, rc, order = 6 )
    N = len(psis[:,1])
    psi6r_col = dumpreader.dump_col( "psi6r", psis[:,0] )
    psi6i_col = dumpreader.dump_col( "psi6i", psis[:,1] )
    psi6a_col = dumpreader.dump_col( "psi6a", psis[:,2] )

    b.other_cols.append( psi6r_col )
    b.other_cols.append( psi6i_col )
    b.other_cols.append( psi6a_col )

    dumpreader.block_to_dump( b, fname2, "a" )

    if (i%2) == 0:
        print "at t = %d, |<psi6>| = %f" % (b.meta.t, psi_avg[2])
    
    i += 1
    b = d.getblock()


# print "<Psi_6> = (%f, %f), |<Psi_6>| = %f" % (psi_avg[0], psi_avg[1], psi_avg[2])

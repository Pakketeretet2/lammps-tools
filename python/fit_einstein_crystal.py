import lammpstools, dumpreader
import sys, os, numpy as np



def get_einstein_approx( d, ids, T, dims, alpha = 0.9 ):
    """ ! Attempts to extract the chemical potential by determining
          effective spring constaints for given ids and using the
          chemical potential of an Einstein crystal. """

    # Prepare initial average:
    b = d.getblock()
    im = lammpstools.make_id_map(b.ids)

    if ids is None:
        ids = []
        n = lammpstools.neighborize( b, 1.3, 3, method = 0 )
        for ni in n:
            if len(ni) > 5:
                ids.append( ni[0] )

    Ncheck = len(ids)
    xavg   = np.zeros( [Ncheck, 3] )
    k_eff  = np.zeros( Ncheck )
    dx2    = np.zeros( [Ncheck, 4] )

    k_vals = []

    for j in range(0,Ncheck):
        i = im[ids[j]]
        xavg[j] = b.x[i]*alpha + (1 - alpha)*xavg[j]

    b = d.getblock()
    # Now determine in time the relative fluctuations with respect to the
    # moving average lattice position:
    while not d.at_eof:
        im = lammpstools.make_id_map(b.ids)
        k_eff  = np.zeros( Ncheck )
        for j in range(0,Ncheck):
            i = im[ids[j]]

            dx = b.x[i][0] - xavg[j][0]
            dy = b.x[i][1] - xavg[j][1]
            dz = b.x[i][2] - xavg[j][2]

            dx2[j][0] = dx*dx
            dx2[j][1] = dy*dy
            dx2[j][2] = dz*dz
            dx2[j][3] = dx2[j][0] + dx2[j][1] + dx2[j][2]

            k_eff[j] = dims*T / dx2[j][3]
            print("%d %d %g %g" % (j, b.meta.t, dx2[j][3], k_eff[j]))
            
        k_vals.append(k_eff)
        for j in range(0,Ncheck):
            xavg[j] = b.x[ im[ids[j]] ]*alpha + (1 - alpha)*xavg[j]
        b = d.getblock()
    return k_vals

import numpy as np


def make_histogram( yy, y0, y1, Nbins ):
    count = 0.0;
    hist = np.zeros(Nbins, dtype=float)
    dx = (y1 - y0) / (Nbins-1.0)
    a = 1.0 / (dx*Nbins)
    bins = np.zeros(Nbins, dtype=float)
    for i in range(0,Nbins):
        bins[i] = y0 + i*dx
    
    misses = 0
    mean   = 0.0
    modal  = 0.0

    for y in yy:
    
        ibin = int( (y - y0) / dx)
        if ibin >= Nbins or ibin < 0:
            ++misses
        else:
            hist[ibin] += a;
            count += 1.0;
            mean  += y

    if count > 0:
        mean /= count;
    modal_bin = np.argmax( hist )
    modal     = bins[modal_bin]
    
    return bins, hist, mean, modal


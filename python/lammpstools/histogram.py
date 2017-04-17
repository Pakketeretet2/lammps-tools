"""!
\ingroup lammpstools

Produces a histogram in a more intuitive way than numpy does.
"""


import numpy as np



def make_histogram( yy, y0, y1, Nbins ):
    """! Produces a histogram from data in yy.
    
    @param yy    Data to histogram
    @param y0    Lower bound of histogram
    @param y1    Upper bound of histogram
    @param Nbins Number of bins.

    The number of bins, y0 and y1 together implicitly define
    the resolution dy = (y1 - y0) / (Nbins-1). Note that data
    outside of the bracket [y0,y1] is discarded.
    """
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


import numpy as np

def make_histogram( data, x0, x1, nbins, normed = False, return_cdf = False ):
    """! Makes a histogram of given data b.meta.tween ranges x0 and x1 """
    dx = (x1-x0)/(nbins-1.0)
    hist = np.zeros(nbins)
    cdf  = np.zeros(nbins)
    bins = np.linspace(x0,x1,nbins)
    
    added = 0
    for x in data:
        idx = int( (x - x0) / dx )
        added += 1.0
        if idx < 0 or idx >= nbins: continue
        
        hist[idx] += 1.0/dx
        
    for i in range(0,len(hist)):
        hist[i] /= float(added)

    I = 0.0;
    counts = 0
    for i in range(1,len(hist)):
        I += 0.5*dx*(hist[i] + hist[i-1])
        cdf[i] = I
        counts += 1
    if normed:
        for i in range(0,len(hist)):
            cdf[i]  /= I
            hist[i] /= (dx*counts)
    
    if return_cdf:
        return bins, hist, cdf
    else:
        return bins, hist

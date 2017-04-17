"""!
This module my own ribbon analysis stuff. 

\ingroup lammpstools
"""

from ctypes import *
from typecasts import *
from dumpreader import *
import numpy as np



def insideness(b,dims = 3):
    """Determines the insideness of each particle. """
    insideness = np.zeros(b.meta.N, dtype=float)
    liblammpstools = cdll.LoadLibrary('liblammpstools.so')
    gi = liblammpstools.get_insideness
    
    gi( void_ptr(b.x), c_longlong(b.meta.N), void_ptr(b.ids),
        void_ptr(b.types), c_longlong(b.meta.domain.periodic),
        void_ptr(b.meta.domain.xlo), void_ptr(b.meta.domain.xhi),
        c_longlong(dims), void_ptr(insideness) )

    return insideness                
    

def edt(b,R,dims=3):
    """Determines the Euclidian distance transform of each particle. """
    edt = np.zeros(b.meta.N, dtype=float)
    lammpstools = cdll.LoadLibrary('liblammpstools.so')
    gedt = lammpstools.get_euclidian_distance_transform
    
    ins = insideness(b,dims)
    
    gedt( void_ptr(b.x), c_longlong(b.meta.N), void_ptr(b.ids),
          void_ptr(b.types), c_longlong(b.meta.domain.periodic),
          void_ptr(b.meta.domain.xlo), void_ptr(b.meta.domain.xhi),
          c_longlong(dims),
          c_double(R), void_ptr(ins), void_ptr(edt) )


    return edt

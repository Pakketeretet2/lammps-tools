"""! @package typecasts

This package contains some functions that help with interfacing C++ and Python.
The interface is based on ctypes.

\ingroup lammpstools

"""


import ctypes
import numpy as np

def ctype_arr_to_np(ptr,Nx,Ny):
    # check if memory addresses are contiguous
    first = pnt[0]
    return np.ctypeslib.as_array(first, (Nx,Ny))


def void_ptr(a):
    "Casts a to a c_type void_ptr"
    return a.ctypes.data_as(ctypes.c_void_p)

"""!
\defgroup lammpstools
This module contains all Python-side tools for analysis of LAMMPS data.

\package lammpstools
This package defines the python interface to the C++ lib and some
(higher level) analysis tools written in Python itself.

\ingroup lammpstools
"""

import math
import struct
import sys
import os
import multiprocessing
from ctypes import c_longlong, c_double

import numpy as np


from typecasts import void_ptr


# Import all submodules:
#import histogram
#import normal_mode_analysis
#import melting_analysis
#import minimum_rmsd
#import fit_einstein_crystal
#import multiprocessing
#import ribbon_analysis
#import compute_com
#import bond_analysis
#import potentials
#import makepairtable
import bond_analysis
import block_data
import lammpstools

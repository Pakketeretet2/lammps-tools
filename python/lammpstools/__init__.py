"""!
\defgroup lammpstools
This module contains all Python-side tools for analysis of LAMMPS data.

\package lammpstools
This package defines the python interface to the C++ lib and some
(higher level) analysis tools written in Python itself.

\ingroup lammpstools
"""

from lammpstools.dumpreader_cpp import *
from lammpstools.bond_analysis  import *
from lammpstools.block_data     import *
from lammpstools.util           import *



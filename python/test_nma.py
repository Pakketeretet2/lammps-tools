# tests the communication stuff.

import lammpstools
from ctypes import *
from typecasts import *
import dumpreader


from normal_mode_analysis import pass_blocks_to_clib
from normal_mode_analysis import normal_mode_analysis


d = dumpreader.dumpreader("traj.dump", x_tag = "f_x[1]", y_tag = "f_x[2]", z_tag = "f_x[3]" )
b = d.getblock()
pass_blocks_to_clib( d )






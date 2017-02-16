#!/bin/zsh
#
#
#

make -j4
make -j4 -f Makefile_no_cgal
make -j4 -f Makefile_static
make -j4 -f Makefile_no_cgal_static
make install
make install -f Makefile_no_cgal 

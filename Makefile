# Makefile for all

.PHONY = all dox clean-dox

all :
	$(MAKE) -C c_lib

install :
	$(MAKE) -C c_lib install

dox :
	doxygen

clean-dox :
	rm -rf html
	rm -rf latex

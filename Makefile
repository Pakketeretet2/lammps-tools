# Makefile for all

.PHONY = all dox clean-dox

all :
	$(MAKE) -C c_lib

install :
	$(MAKE) -C c_lib install

dox :
	cd python; doxygen;
	ln -s --force ./python/html/index.html manual.html

clean-dox :
	rm -rf html
	rm -rf latex

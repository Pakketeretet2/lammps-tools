# Makefile for all

.PHONY = all dox clean-dox

all :
	$(MAKE) -C c_lib $@
	$(MAKE) -C c_lib/test_lib $@

dox :
	doxygen

clean-dox :
	rm -rf html
	rm -rf latex

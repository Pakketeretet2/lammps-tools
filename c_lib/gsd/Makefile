# A simple Makefile for the gsd library.
#

SRC = gsd.c
EXE = libgsd.so


EXE : $(SRC)
	gcc -O3 -g -fPIC -shared -o $(EXE) $(SRC)

install :
	cp $(EXE) /usr/local/lib/$(EXE)


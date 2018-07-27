# Simple Makefile to compile the cc programs.

ifndef CXX
	CC=g++
endif

ifndef OMP_DIRECTIVE
	OMP_DIRECTIVE=-fopenmp
endif

sc_green.cc.x: sc_green.cc
	$(CC) $(OMP_DIRECTIVE) sc_green.cc -o sc_green.cc.x

sc_measure.cc.x: sc_measure.cc
	$(CC) $(OMP_DIRECTIVE) sc_measure.cc -o sc_measure.cc.x

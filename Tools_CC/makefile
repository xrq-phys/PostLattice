# Simple Makefile to compile the cc programs.

ifndef CXX
	CC=g++
endif

ifndef OMP_DIRECTIVE
	OMP_DIRECTIVE=-fopenmp
endif

default: sc_measure.cc.x sc_green.cc.x

sc_green.cc.x: sc_green.cc sc_classes.hh
	$(CXX) $(OMP_DIRECTIVE) sc_green.cc -o sc_green.cc.x

sc_measure.cc.x: sc_measure.cc sc_classes.hh
	$(CXX) $(OMP_DIRECTIVE) sc_measure.cc -o sc_measure.cc.x


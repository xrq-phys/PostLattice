# Simple Makefile to compile the cc programs.

ifndef CXX
	CXX=g++
endif

ifndef OMP_DIRECTIVE
	OMP_DIRECTIVE=-fopenmp
endif

default: sc_measure.cc.x sc_green.cc.x

%.cc.x: %.cc sc_classes.hh
	$(CXX) $(OMP_DIRECTIVE) $< -o $@


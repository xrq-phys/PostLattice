# Simple Makefile to compile the cc programs.

ifndef CXX
	CXX=g++
endif

ifndef OMP_DIRECTIVE
	OMP_DIRECTIVE=-fopenmp
endif

ifndef EXTRA_FLAGS
	EXTRA_FLAGS=-O3
endif

default: sc_measure.cc.x sc_green.cc.x core.x

%.cc.x: %.cc sc_classes.hh
	$(CXX) $(OMP_DIRECTIVE) $(EXTRA_FLAGS) $< -o $@

core.x: lattice.o green_gen.o operators.o trig_lattice.o core.o
	$(CXX) $(OMP_DIRECTIVE) $(EXTRA_FLAGS) -o core.x lattice.o green_gen.o operators.o trig_lattice.o

%.o: %.cc
	$(CXX) $(OMP_DIRECTIVE) $(EXTRA_FLAGS) -c -std=c++17 $<


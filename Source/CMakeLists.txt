ADD_LIBRARY (PostLattice
	SHARED
	lattice.cc operators.cc trig_lattice.cc
	lattice.hh operators.hh trig_lattice.hh)

ADD_EXECUTABLE (Post
	green_gen.hh
	green_gen.cc
	core.cc)

TARGET_LINK_LIBRARIES (Post PostLattice)

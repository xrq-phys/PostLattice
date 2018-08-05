ADD_LIBRARY (PostLattice
	SHARED
	lattice.cc operators.cc trig_lattice.cc
	lattice.hh operators.hh trig_lattice.hh)

ADD_EXECUTABLE (Post
	green_gen.hh measure.hh
	green_gen.cc measure.cc
	core.cc)

TARGET_LINK_LIBRARIES (Post PostLattice inih)

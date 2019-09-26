ADD_LIBRARY (PostLattice
	STATIC
	lattice.cc operators.cc trig_lattice.cc honeycomb.cc lieb.cc
	lattice.hh operators.hh trig_lattice.hh honeycomb.hh lieb.hh)

ADD_EXECUTABLE (Post
	green_gen.hh measure.hh plot.cc
	green_gen.cc measure.cc plot.hh
	core.cc)

TARGET_LINK_LIBRARIES (Post PostLattice inih)
INSTALL (TARGETS Post PostLattice
	RUNTIME DESTINATION bin
	ARCHIVE DESTINATION lib)

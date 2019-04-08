/**
 * @file plot.hh
 * Lattice plotter.
 */
#pragma once
#include "lattice.hh"

/**
 * @brief Plot lattice geometry to target file.
 * 
 * @param mp_fnm Target METAPOST filename.
 * @param system Lattice system object.
 */
void plot_lattice(const char *mp_fnm, lattice::lattice &system, const bool label_on);

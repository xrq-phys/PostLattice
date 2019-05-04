/**
 * @file plot.hh
 * Lattice plotter.
 */
#pragma once
#include "lattice.hh"

double site_func_const(double *val_mat, int idx) { return 4; }

/**
 * @brief Plot lattice geometry to target file.
 * 
 * @param mp_fnm Target METAPOST filename.
 * @param system Lattice system object.
 */
void plot_lattice(const char *mp_fnm, lattice::lattice &system, const bool label_on,
                  double (*size_func)(double *, int) = site_func_const, double *size_mat = nullptr,
                  double (*spin_func)(double *, int) = site_func_const, double *spin_mat = nullptr);

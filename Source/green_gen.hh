/**
 * @file green_gen.hh
 * Green function generator.
 */
#include "lattice.hh"

/**
 * @brief Write 1-body Green function definition used by AF/SC analysis.
 * 
 * @param ca_fnm 1-body Green function definition filename.
 * @param system Lattice system object.
 */
void write_ca(const char *ca_fnm, lattice::lattice &system);

/**
 * @brief Write 2-body Green function definition used by AF/SC analysis.
 * 
 * @param caca_fnm 2-body Green function definition filename.
 * @param system Lattice system object.
 */
void write_caca(const char *caca_fnm, lattice::lattice &system,
        bool sc = true, bool sc_p = false, bool sc_simp = false);

/**
 * @file measure.hh
 * Main routine for performing measurements.
 */
#pragma once
#include "lattice.hh"
#include "operators.hh"

/**
 * @brief Options for operators being evaluated.
 */
struct operator_options {
    const char *g1e_fnm = nullptr;
    const char *g2e_fnm = nullptr;
    bool db = true;
    bool af = true;
    char sc = '-'; ///< Waveform
    int  sc_n = 4;
    int *af_n = nullptr;
    bool verbose = false;
    const char *sc_fnm = nullptr;
    const char *af_fnm = nullptr;
    const char *db_fnm = nullptr;

    ~operator_options()
    { delete[] af_n; }
};

/**
 * @brief Main extrance for measurement.
 *
 * @param physics Lattice system.
 * @param options Options struct.
 */
void measure(lattice::lattice physics, operator_options &options);

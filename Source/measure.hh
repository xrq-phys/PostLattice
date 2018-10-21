/**
 * @file measure.hh
 * Main routine for performing measurements.
 */
#pragma once
#include "lattice.hh"
#include "operators.hh"
#include <string>

/**
 * @brief Options for operators being evaluated.
 */
struct operator_options {
    std::string g1e_fnm = "";
    std::string g2e_fnm = "";
    bool db = true;
    bool af = true;
    char sc = '-'; ///< Waveform
    char sc_stat = 'S'; ///< Statistic mode
    int  sc_n = 4;
    int *af_n = nullptr;
    bool direct = false;
    bool verbose = false;
    std::string sc_fnm = "";
    std::string af_fnm = "";
    std::string db_fnm = "";

    ~operator_options()
    { delete[] af_n; }
};

/**
 * @brief Main extrance for measurement.
 *
 * @param physics Lattice system.
 * @param options Options struct.
 */
void measure(lattice::lattice &physics, operator_options &options);

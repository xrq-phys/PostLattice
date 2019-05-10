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
    std::string g1e_idx_fnm = "";
    std::string g2e_idx_fnm = "";
    std::string sc = "-"; ///< Waveform, or specified manually.
    bool db = true;
    bool af = true;
    char sc_stat = 'S'; ///< Statistic mode
    int  sc_use_p = 0; ///< Parallel SC pairing
    int  sc_simple = 0; ///< Count SC correlation only from site index-0
    int  sc_n = 4;
    int *af_n = nullptr;
    bool direct = false;
    bool verbose = false;
    bool chk_duplicate = true;
    bool ca_verbose = true; ///< Green function output has i|j|k|l|x or x|x|x|... form
    bool ca_legacy = true; ///< Green function indices has r|s|r|s or r|r|s form
    std::string sc_fnm = "";
    std::string af_fnm = "";
    std::string db_fnm = "";
    std::string mp_fnm = "";

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

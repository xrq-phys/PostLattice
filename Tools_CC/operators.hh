/**
 * @file operators.hh
 * General operator interface and some instances.
 */
#pragma once
#include "lattice.hh"

namespace operators
{
    /**
     * @brief Superconductive correlation function.
     */
    struct sc_corr
    {
        lattice::lattice &system;
        char form; ///< Wave form factor of superconductivity: s, d or p.
        int rc_count; ///< Number of R_c's to use.
        double *values;

        sc_corr(lattice::lattice &system_i, const int rc_count_i, const char form_i)
        : system(system_i), rc_count(rc_count_i), form(form_i)
        { values = new double[rc_count]; for (int i = 0; i < rc_count; i++) values[i] = 0; }

        ~sc_corr()
        { delete[] values; }

        /**
         * @brief Do measurement.
         * 
         * @param a  Site index.
         * @param sa Spin index.
         * @param i  Site index.
         * @param si Spin index.
         * @param b  Site index.
         * @param sb Spin index.
         * @param j  Site index.
         * @param sj Spin index.
         * @param x  < c c a a >.
         */
        void measure(int a, int sa, int b, int sb, int i, int si, int j, int sj, double x);

        /**
         * @brief Validation if a term MIGHT have contribution to SC correlation.
         * 
         * @param a  Site index.
         * @param sa Spin index.
         * @param i  Site index.
         * @param si Spin index.
         * @param b  Site index.
         * @param sb Spin index.
         * @param j  Site index.
         * @param sj Spin index.
         * @return If the term exists.
         */
        bool validate(int a, int sa, int b, int sb, int i, int si, int j, int sj)
        { return si != sj && sa != sb && system.nn[a][b] && system.nn[i][j]; }
    };

    struct spin_struct
    {
        lattice::lattice &system;
        int* points; ///< Q-poitns.
        double *values; /// Measured values.
        // Pauli matrices.
        int pauli_x[2][2];
        int pauli_y[2][2];
        int pauli_z[2][2];

        // NB: USE DOUBLE R ROUTINES WHEN EVALUATING
    };
}

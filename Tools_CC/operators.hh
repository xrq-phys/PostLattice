/**
 * @file operators.hh
 * General operator interface and some instances.
 */
#pragma once
#include "lattice.hh"

namespace operators
{
    /**
     * @brief Operator base class.
     */
    struct operators
    {
        double *values; ///< Measured values.

        operators()
        : values(nullptr) { }

        ~operators()
        { if (values != nullptr) delete[] values; }

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
         * @param x  < c c a a > or < c a c a > depending on form of the operator.
         */
        virtual void measure(int a, int sa, int b, int sb, int i, int si, int j, int sj, double x);
    };

    /**
     * @brief Superconductive correlation function.
     */
    struct sc_corr
    : operators
    {
        lattice::lattice &system;
        char form; ///< Wave form factor of superconductivity: s, d or p.
        int rc_count; ///< Number of R_c's to use.

        sc_corr(lattice::lattice &system_i, const int rc_count_i, const char form_i)
        : operators::operators(), system(system_i), rc_count(rc_count_i), form(form_i)
        { values = new double[rc_count]; for (int i = 0; i < rc_count; i++) values[i] = 0; }

        /**
         * @brief Implement measurement.
         * @param x < c c a a >
         */
        virtual void measure(int a, int sa, int b, int sb, int i, int si, int j, int sj, double x) override;

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

    /**
     * @brief Spin structure.
     */
    struct spin_struct
    : operators
    {
        lattice::lattice &system;
        int n_points;
        int **points; ///< Q-points.
        // Pauli matrices.
        int pauli_x[2][2];
        int pauli_y[2][2];
        int pauli_z[2][2];

        /**
         * @brief Construct a new spin struct object.
         * 
         * @param system_i System.
         * @param ndiv Number of q-points in each dimension.
         */
        spin_struct(lattice::lattice system_i, int *ndiv);

        ~spin_struct()
        { for (int i = 0; i < n_points; i++) delete[] points[i]; delete[] points; }

        /**
         * @brief Implement measurement.
         * @param x < c a c a >
         */
        virtual void measure(int a, int sa, int b, int sb, int i, int si, int j, int sj, double x) override;
    };
}

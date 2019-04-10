/**
 * @file operators.hh
 * General operator interface and some instances.
 */
#pragma once
#include "lattice.hh"
#include <complex>

namespace operators
{
    /**
     * @brief Operator base class.
     */
    struct operators
    {
        std::complex<double> *values; ///< Measured values.

        operators()
        : values(nullptr) { }

        virtual ~operators()
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
         * @param x  < c c a a > or < c a c a > depending on form of the operator.
         */
        virtual void measure(int a, int sa, int b, int sb, int i, int si, int j, int sj, std::complex<double> x);

        /**
         * @brief Refresh operator data so that it's ready for output.
         */
        virtual void refresh() { }
    };

    /**
     * @brief Double occupation.
     */
    struct doublon
    : operators
    {
        doublon()
        : operators() { values = new std::complex<double>; }

        /**
         * @brief Implement measurement.
         * @param x < c a c a >
         */
        virtual void measure(int a, int sa, int b, int sb, int i, int si, int j, int sj, std::complex<double> x) override;
    };

    /**
     * @brief Superconductive correlation function.
     */
    struct sc_corr
    : operators
    {
        lattice::lattice &system;
        std::complex<double> *val_mat; ///< Storage for SC correlation at index difference i.
        std::complex<double> *form; ///< Wave form factor of superconductivity.
        // char form; ///< Wave form factor of superconductivity: s, d or p. TODO: Delete it.
        int use_p; ///< Whether to calculate parallel spin paring: 0AP, 1P.
        int rc_count; ///< Number of R_c's to use.

        sc_corr(lattice::lattice &system_i, const int rc_count_i, const int use_p_i, std::string form_s);

        virtual ~sc_corr()
        { delete[] val_mat; delete[] form; }

        /**
         * @brief Implement measurement.
         * @param x < c c a a >
         */
        virtual void measure(int a, int sa, int b, int sb, int i, int si, int j, int sj, std::complex<double> x) override;

        virtual void refresh() override;

        virtual void refresh(char mode);

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
        { return (use_p ? si == sj : si != sj && sa != sb) && system.nn[a][b] && system.nn[i][j]; }
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

        /**
         * @brief Construct a new spin struct object.
         * 
         * @param system_i System.
         * @param ndiv Number of q-points in each dimension.
         */
        spin_struct(lattice::lattice &system_i, const int *ndiv);

        virtual ~spin_struct() override
        { for (int i = 0; i < n_points; i++) delete[] points[i]; delete[] points; }

        /**
         * @brief Implement measurement.
         * @param x < c a c a >
         */
        virtual void measure(int a, int sa, int b, int sb, int i, int si, int j, int sj, std::complex<double> x) override;
    };
}

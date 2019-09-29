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

        virtual ~operators()
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
        : operators() { values = new double; }

        /**
         * @brief Implement measurement.
         * @param x < c a c a >
         */
        virtual void measure(int a, int sa, int b, int sb, int i, int si, int j, int sj, double x) override;
    };

    /**
     * @brief Superconductive correlation function.
     */
    struct sc_corr
    : operators
    {
        lattice::lattice &system;
        double *val_mat; ///< Storage for SC correlation at index difference i.
        char form; ///< Wave form factor of superconductivity: s, d or p.
        int use_p; ///< Whether to calculate parallel spin paring: 0AP, 1P.
        int rc_count; ///< Number of R_c's to use.

        sc_corr(lattice::lattice &system_i, const int rc_count_i, const char form_i, const int use_p_i)
        : operators::operators(), system(system_i), form(form_i), use_p(use_p_i),
          rc_count(rc_count_i < system_i.rc_n ? rc_count_i : system_i.rc_n)
        { values = new double[rc_count];
          for (int i = 0; i < rc_count; i++) values[i] = 0;
          val_mat = new double[system_i.n * system_i.ncell];
          for (int i = 0; i < system_i.n * system_i.ncell; i++) val_mat[i] = 0; }

        virtual ~sc_corr()
        { delete[] val_mat; delete[] values; val_mat = nullptr; values = nullptr; }

        /**
         * @brief Implement measurement.
         * @param x < c c a a >
         */
        virtual void measure(int a, int sa, int b, int sb, int i, int si, int j, int sj, double x) override;

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
     * @brief Base class for 2-site correlation.
     */
    struct site_corr
    : operators
    {
        lattice::lattice &system;
        int **connection; ///< Correlation Site Indices
        int conn_slave;
        double *val_mat;

        /**
         * @brief Construct a new 2-site correlation object.
         * 
         * @param system_i System.
         * @param ndiv Number of q-points in each dimension.
         */
        site_corr(lattice::lattice &system_i);

        /**
         * @brief Construct a new 2-site correlation object from the same system,
         *        copying connection data from previous object.
         *
         * @param opr_i Operator object to copy connection data from.
         */
        site_corr(site_corr &opr_i);

        virtual ~site_corr() override;

        virtual void refresh() override {}

        /**
         * @brief Calculate real-space correlation vs. R.
         * @param mode Average or find maximum.
         * @param rc_count Number of R's to count.
         */
        virtual void refresh(char mode, int rc_count);

        /**
         * @brief Perform Fourier transformation and stores k-space vectors to values.
         * @param ndiv Number of Q-Points per dimension.
         */
        virtual void refresh(int *ndiv);
    };

    /**
     * Spin-structure factor - 2-point spin inner product.
     */
    struct spin_struct
    : site_corr
    {
        spin_struct(lattice::lattice &system_i)
        : site_corr(system_i) { }

        spin_struct(site_corr &opr_i)
        : site_corr(opr_i) { }

        /**
         * @brief Implement measurement.
         * @param x < c a c a >
         */
        virtual void measure(int a, int sa, int b, int sb, int i, int si, int j, int sj, double x) override;
    };

    /**
     * Charge-structure factor - 2-point density correlation.
     */
    struct charge_struct
    : site_corr
    {
        charge_struct(lattice::lattice &system_i)
        : site_corr(system_i) { }

        charge_struct(site_corr &opr_i)
        : site_corr(opr_i) { }

        /**
         * @brief Implement measurement.
         * @param x < c a c a >
         */
        virtual void measure(int a, int sa, int b, int sb, int i, int si, int j, int sj, double x) override;
    };
}

/**
 * @file trig_lattice.hh
 * Triangular Lattice
 */
#pragma once
#include "lattice.hh"

namespace lattice{
    struct trig2d
    : lattice::lattice
    {
        int a0W, ///< Vector a0W with a0L = 0 (otherwise not supported here).
            a1L, a1W; ///< Vector a1.
        int a1Lr, a1Wr; ///< Smallest integral direction.

        trig2d(int a0W_i, int a1L_i, int a1W_i);

        /**
         * @brief Implementation of r. Note that here r in diamond coordinate is returned.
         * 
         * @param r Location will be written here. DIM=2
         * @param i Index.
         * @param r_q Block ID. DIM=2
         */
        virtual void r(int *r, int i, int *r_q = nullptr) override;

        /**
         * @brief Implementation of real r.
         * 
         * @param r Location will be written here. DIM=2
         * @param i Index.
         */
        virtual void r(double *r, int i) override;

        virtual int calc_rmin(int i, int j) override;

        /**
         * @brief Get integral x and y coordinate in diamond x-y system.
         * 
         * @param x x INDEX.
         * @param y y INDEX/COORDINATE.
         * @return int x COORDINATE
         */
        inline int x_diam(int x, int y)
        { return x + ((y + (a1W > 0 ? 1 : 0)) * a1W) / a1L; }

        /**
         * @brief Get x-index of diamond x-y coordinate.
         * 
         * @param x x COORDINATE.
         * @param y y COORDINATE.
         * @return int x INDEX.
         */
        inline int x_idx(int x, int y)
        { return x - ((y + (a1W > 0 ? 1 : 0)) * a1W) / a1L; }
    };
}
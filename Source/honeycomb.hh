/**
 * \file honeycomb.hh
 * Honeycomb Lattice
 */
#pragma once
#include "lattice.hh"

namespace lattice {
    struct honeycomb
    : lattice::lattice
    {
        int a0W, ///< Vector a0W with a0L = 0 (otherwise not supported here).
            a1L, a1W; ///< Vector a1.

        honeycomb(int a0W_i, int a1L_i, int a1W_i);

        /**
         * @brief Implementation of r.
         *        Note that the CELL diamond coordinate instead of lattice position is returned.
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

        /**
         * @brief Implementation of equivalent index finder.
         *        Note that if the indice is negative, one should start from site 1 instead of 0.
         * 
         * @param x Site i.
         * @param y Site j.
         * @return int The equivalent site index.
         */
        virtual int idx_rij(int x, int y) override;

        virtual int calc_rmin(int i, int j) override;

        /**
         * @brief Index to coordinate converter for X (cell index, not lattice).
         *
         * @param x x INDEX.
         * @param y y INDEX/COORDINATE.
         * @return int x COORDINATE
         */
        inline int x_diam(int x, int y)
        { return x + ((y + (a1W > 0 ? 1 : 0)) * a1W) / a1L; }

        /**
         * @brief Coordinate to index converter for X (cell index).
         *
         * @param x x COORDINATE.
         * @param y y COORDINATE.
         * @return int x INDEX.
         */
        inline int x_idx(int x, int y)
        { return x - ((y + (a1W > 0 ? 1 : 0)) * a1W) / a1L; }
    };
}
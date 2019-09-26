/**
 * \file lieb.hh
 * Lieb Lattice (square base only)
 */
#pragma once
#include "lattice.hh"

namespace lattice {
    struct lieb
    : lattice
    {
        int w, l; ///< Number of row/columns.

        lieb(int w_i, int l_i);

        /**
         * \brief Implementation of integral R.
         *
         * @param r Output.
         * @param i Index.
         * @param r_q Block indices.
         */
        virtual void r(int *r, int i, int *r_q = nullptr) override;

        /**
         * \brief Real-valued R.
         *
         * @param r Output.
         * @param i Index.
         */
        virtual void r(double *r, int i) override;

        /**
         * \brief Equivalent index finder.
         *        Here's situation here: I don't git a s**t about vectors starting from
         *        non-fully connected sites. So finally the value should be multiplied
         *        by a factor of 3 as it should be P(SC)/Ncells instead of P(SC)/Nsites.
         *
         * @param i starting point
         * @param j end point
         * @return site index starting from 0. Of course.
         */
        virtual int idx_rij(int i, int j) override;

        /**
         * \brief Minimal R^2.
         *        Note that this only calculates from 0 as in this version only distance
         *        form 0 is used. If no one proposed extension that utilizes it,
         *        Removing index i will be done in future.
         *
         * @param i From index.
         * @param j To index.
         * @return Square of R.
         */
        virtual int calc_rmin(int i, int j) override;
    };
}


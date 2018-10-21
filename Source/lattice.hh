/**
 * \file lattice.hh
 * Basic facilities & some common instances for lattice systems.
 */
#pragma once 

namespace lattice 
{
    struct lattice
    {
        int n;
        int dim;    ///< Dimension.
        int rc_n;   ///< Number of R_c to Measure.
        int *r_c;   ///< Minimal-r List.
        int *r_n;   ///< Number of Lattice Points at Each R_c.
        int **rmin; ///< Minimal r LUT.
        int **nn;   ///< Nearest-Neighhbour LUT.

        lattice(int n_i, int dim_i)
        : n(n_i), dim(dim_i)
        {
            nn = new int*[n];
            rmin = new int*[n];
            for (int i = 0; i < n; i++) {
                nn[i] = new int[n];
                rmin[i] = new int[n];
                for (int j = 0; j < n; j++) {
                    nn[i][j] = 0;
                    rmin[i][j] = -1;
                }
            }
        }

        virtual ~lattice()
        { free_lut(); delete[] r_c; delete[] r_n; }

        /**
         * @brief Query integral location of site index i.
         * 
         * @param r Location will be written here. Make sure it's >= DIM.
         * @param i Index.
         * @param r_q If the region's sliced into blocks, return block id in each direction.
         */
        virtual void r(int *r, int i, int *r_q = nullptr);

        /**
         * @brief Query real fractional location of site index i.
         * 
         * @param r Location will be written here. Make sure it's >= DIM.
         * @param i Index
         */
        virtual void r(double *r, int i);

        /**
         * @brief Returns site index corresponding to the lattice vector
         *        which is the distance between site i and j.
         * @param i Site i.
         * @param j Site j.
         * @return The index.
         */
        virtual int idx_rij(int i, int j);

        /**
         * @brief Calculates minimal (integral) distance squared between two sites.
         * 
         * @param i Site 1.
         * @param j Site 2.
         * @return Distance.
         */
        virtual int calc_rmin(int i, int j);

        /**
         * @brief Frees default look-up-tables.
         */
        inline void free_lut()
        {
            for (int i = 0; i < n; i++) {
                delete[] nn[i];
                delete[] rmin[i];
            }
            delete[] nn;
            delete[] rmin;
        }
    };

    struct square2d
    : lattice::lattice
    {
        int w; ///< Number of columns (x-length).
        int l; ///< Number of rows    (y-length).

        square2d(int w_i, int l_i);

        /**
         * @brief Implementation of integral r.
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
         * @param r_q Block ID. DIM=2
         */
        virtual void r(double *r, int i) override;

        virtual int idx_rij(int x, int y) override;

        /**
         * @brief Implementation of prototype.
         * 
         * @param i Site 1.
         * @param j Site 2.
         * @return Distance.
         */
        virtual int calc_rmin(int i, int j) override;

        /**
         * @brief Returns site index for location.
         * @param x x location.
         * @param y y location.
         */
        int idx_2d(int x, int y) 
        { return y*w + x; }
    };

    /**
     * @brief Initializes The System's R_c For SC Calculation. 
     *        Sets R_c list, number of sites with R=R_c and total number of R_c.
     * @param system The lattice system object.
     */
    void sorted_rc(lattice &system);
}

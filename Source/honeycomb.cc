#include "honeycomb.hh"
#ifdef _DEBUG
#include <iostream>
#endif

// {
// These are Tool Functions.
// TODO: Put it somewhere else.
int inner2_honeycomb(const int *a, const int ai, const int *b, const int bi)
{
    // AI and BI are site indices in the cell parsed by A and B.
    int da1 = (a[0] - b[0]) * 3,
        da2 = (a[1] - b[1]) * 3;
    int l2;
    if (((ai % 2) + (bi % 2)) % 2) {
        if (ai % 2) {
            da1 += 1;
            da2 += 1;
        } else {
            da1 -= 1;
            da2 -= 1;
        }
    }
    // This sum of product must be a multiply of 3.
    l2 = da1 * da1 + da2 * da2 + da1 * da2;
#ifdef _DEBUG
    if (l2 % 3)
        std::cerr << "Err: Undefined honeycomb inner product." << std::endl;
#endif
    return l2 / 3;
}
// }

lattice::honeycomb::honeycomb(int a0W_i, int a1L_i, int a1W_i)
: a0W(a0W_i), a1L(a1L_i), a1W(a1W_i), lattice::lattice(a0W_i * a1L_i * 2, 2, 2)
{
    sorted_rc(*this);

    // Update LUT.
    int ri[2], ynn, xnn; ///< They are all CELL indices and coordinates.
    for (int i = 1; i < n; i+=2) { // Work on odd indices.
        r(ri, i);

        // In-cell connection, two-way binding.
        nn[i][i - 1] = 1;
        nn[i - 1][i] = 1;

        // X + 1 direction, two-way binding.
        xnn = x_idx(ri[0], ri[1]);
        nn[i][(a0W * ri[1] + (xnn == a0W - 1 ? 0 : xnn + 1)) * 2] = 2;
        nn[(a0W * ri[1] + (xnn == a0W - 1 ? 0 : xnn + 1)) * 2][i] = 2;
        // Y + 1, two-way.
        ynn = ri[1] + 1;
        if (ynn == a1L) {
            ynn = 0;
            // x-index is conserved.
            xnn = x_idx(ri[0], ri[1]);
            // Boundary.
            if (a1W % a1L) xnn++;
        } else {
            // x-coordinate conserved.
            xnn = x_idx(ri[0], ynn);
        }
        xnn = ((xnn % a0W) + a0W) % a0W;
        nn[i][(a0W * ynn + xnn) * 2] = 3;
        nn[(a0W * ynn + xnn) * 2][i] = 3;
    }
}

void lattice::honeycomb::r(int *r, int i, int *r_q)
{
    int icl = i / 2; // Cell index.
    int xi = icl % a0W, yi = icl / a0W; // x-y INDEX, not coordinate.

    if (r_q != nullptr) {
        if (xi < a0W / 2)
            r_q[0] = 0;
        else 
            r_q[0] = 1;
        if (yi < a1L / 2)
            r_q[1] = 0;
        else 
            r_q[1] = 1;
    }

    // Set coordinates.
    r[1] = yi;
    r[0] = x_diam(xi, yi);
}

void lattice::honeycomb::r(double *r, int i)
{
    int ri[2];
    this->r(ri, i);
    r[1] = (double(ri[1]) + (i % 2 ? 1. / 3 : 0)) / a1L;
    r[0] = (ri[0] - r[1] * a1W + (i % 2 ? 1. / 3 : 0)) / a0W;
}

int lattice::honeycomb::idx_rij(int i, int j)
{
    int idc, k, shift;
    int ri[2],
        rj[2];
    r(ri, i);
    r(rj, j);
    ri[0] -= rj[0];
    ri[1] -= rj[1];
    // The horizontal PBC.
    if (ri[1] < 0) {
        ri[1] += a1L;
        ri[0] += a1W;
    }
    // The slant-vertical PBC.
    if (ri[0] < x_diam(0, ri[1]))
        ri[0] += a0W;
    idc = ri[1] * a0W + x_idx(ri[0], ri[1]);
    if (j % 2 && (i + 1) % 2)
        return idc * 2 + n;
    else if (i % 2 && (j + 1) % 2)
        return idc * 2 + 1;
    else
        return idc * 2;
}

int lattice::honeycomb::calc_rmin(int i, int j) 
{
    if (rmin[i][j] >= 0)
        return rmin[i][j];
    else {
        double ri_f[2], rj_f[2];
        int ri[2], rj[2];
        int rtmp, rnew;
        r(ri, i);
        r(rj, j);
        r(ri_f, i);
        r(rj_f, j);

        // [0, 0]
        rtmp = inner2_honeycomb(ri, i, rj, j);

        // X Displacement.
        if (ri_f[0] < rj_f[0])
            ri[0] += a0W;
        else
            rj[0] += a0W;
        rnew = inner2_honeycomb(ri, i, rj, j);
        if (rtmp > rnew)
            rtmp = rnew;

        // XY Displacement.
        if (ri_f[1] < rj_f[1]) {
            ri[0] += a1W;
            ri[1] += a1L;
        } else {
            rj[0] += a1W;
            rj[1] += a1L;
        }
        rnew = inner2_honeycomb(ri, i, rj, j);
        if (rtmp > rnew)
            rtmp = rnew;

        // Y Displacement.
        if (ri_f[0] < rj_f[0])
            ri[0] -= a0W;
        else
            rj[0] -= a0W;
        rnew = inner2_honeycomb(ri, i, rj, j);
        if (rtmp > rnew)
            rtmp = rnew;

        // Update LUT.
        rmin[i][j] = rtmp;
        rmin[j][i] = rtmp;

        return rtmp;
    }
}
#include "trig_lattice.hh"
#include <cstring>

// {
// These are Tool Functions.
// TODO: Put it somewhere else.
int inner2_trig(const int *a, const int *b)
{
    return (a[0] - b[0]) * (a[0] - b[0]) +
           (a[1] - b[1]) * (a[1] - b[1]) +
           (a[1] - b[1]) * (a[0] - b[0]);
}
// }

lattice::trig2d::trig2d(int a0W_i, int a1L_i, int a1W_i)
: a0W(a0W_i), a1L(a1L_i), a1W(a1W_i), lattice::lattice(a0W_i * a1L_i, 2, 1, 6)
{
    sorted_rc(*this);

    // Update LUT.
    int ri[2], ynn, xnn;
    for (int i = 0; i < n; i++) {
        r(ri, i);

        // 0 Deg direction.
        xnn = x_idx(ri[0], ri[1]);
        nn[i][a0W * ri[1] + (xnn == 0 ? a0W - 1 : xnn - 1)] =-1;
        nn[i][a0W * ri[1] + (xnn == a0W - 1 ? 0 : xnn + 1)] = 1;
        // Y + 1.
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
        nn[i][a0W * ynn + xnn] = 2;
        nn[i][a0W * ynn + (xnn == 0 ? a0W - 1 : xnn - 1)] = 3;
        // Y - 1.
        // Do the same.
        ynn = ri[1] - 1;
        if (ynn < 0) {
            ynn = a1L - 1;
            xnn = x_idx(ri[0], ri[1]);
            if (a1W % a1L) xnn--;
        } else {
            xnn = x_idx(ri[0], ri[1]);
        }
        xnn = ((xnn % a0W) + a0W) % a0W;
        nn[i][a0W * ynn + xnn] =-2;
        nn[i][a0W * ynn + (xnn == a0W - 1 ? 0 : xnn + 1)] =-3;
    }
}

void lattice::trig2d::r(int *r, int i, int *r_q)
{
    int xi = i % a0W, yi = i / a0W; // x-y INDEX, not coordinate.

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

void lattice::trig2d::r(double *r, int i)
{
    int ri[2];
    this->r(ri, i);
    r[0] = (ri[0] - double(ri[1]) * a1W / a1L) / a0W;
    r[1] = double(ri[1]) / a1L;
}

int lattice::trig2d::idx_rij(int i, int j)
{
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
    return ri[1] * a0W + x_idx(ri[0], ri[1]);
}

int lattice::trig2d::calc_rmin(int i, int j) {
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
        rtmp = inner2_trig(ri, rj);

        // X Displacement.
        if (ri_f[0] < rj_f[0])
            ri[0] += a0W;
        else
            rj[0] += a0W;
        rnew = inner2_trig(ri, rj);
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
        rnew = inner2_trig(ri, rj);
        if (rtmp > rnew)
            rtmp = rnew;

        // Y Displacement.
        if (ri_f[0] < rj_f[0])
            ri[0] -= a0W;
        else
            rj[0] -= a0W;
        rnew = inner2_trig(ri, rj);
        if (rtmp > rnew)
            rtmp = rnew;

        // Update LUT.
        rmin[i][j] = rtmp;
        rmin[j][i] = rtmp;

        return rtmp;
    }
}

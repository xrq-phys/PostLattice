#include "lieb.hh"

// Tool Routines {
int inner2(int *a, int *b, bool cell = false)
{ return ((a[0] - b[0]) * (a[0] - b[0])) + ((a[1] - b[1]) * (a[1] - b[1])); }
// }

lattice::lieb::lieb(int w_i, int l_i)
: w(w_i), l(l_i), lattice::lattice(w_i * l_i * 3, 2, 3, 4)
{
    sorted_rc(*this);

    int icell, jcell;
    for (int x = 0; x < w; x++)
        for (int y = 0; y < l; y++) {
            icell = y * w + x;
            // Inner cell.
            nn[icell * ncell][icell * ncell + 1] = 1;
            nn[icell * ncell][icell * ncell + 2] = 2;
            nn[icell * ncell + 1][icell * ncell] =-1;
            nn[icell * ncell + 2][icell * ncell] =-2;
            // Inter cell X.
            jcell = y * w + ((x + 1) % w);
            nn[icell * ncell + 1][jcell * ncell] = 1;
            nn[jcell * ncell][icell * ncell + 1] =-1;
            // Inter cell Y.
            jcell = ((y + 1) % w) * w + x;
            nn[icell * ncell + 2][jcell * ncell] = 2;
            nn[jcell * ncell][icell * ncell + 2] =-2;
        }
}

void lattice::lieb::r(int *r, int i, int *r_q)
{
    int isite = i % ncell, ///< Site index in cell.
        icell = i / ncell;
    int xcell = icell % w,
        ycell = icell / w;
    r[0] = xcell * 2,
    r[1] = ycell * 2;

    r[0] += isite % 2;
    r[1] += isite / 2;

    if (r_q != nullptr) {
        r_q[0] = (xcell >= (w + 1) / 2);
        r_q[1] = (ycell >= (l + 1) / 2);
    }
}

void lattice::lieb::r(double *r, int i)
{
    int r_i[2];
    this->r(r_i, i);
    r[0] = double(r_i[0]) / (w * 2);
    r[1] = double(r_i[1]) / (l * 2);
}

int lattice::lieb::idx_rij(int i, int j) {
    int icell = i / ncell,
        isite = i % ncell,
        jcell = j / ncell,
        jsite = j % ncell;
    int xicl = icell % w,
        yicl = icell / w,
        xjcl = jcell % w,
        yjcl = jcell / w;
    int xrcl = (xjcl - xicl + w) % w,
        yrcl = (yjcl - yicl + l) % l;

    return (yrcl * w + xrcl) * ncell + jsite + isite * n;
}

int lattice::lieb::calc_rmin(int i, int j) {
    // Currently copying square lattice's code.
    if (rmin[i][j] >= 0)
        return rmin[i][j];
    else {
        int ri[2], rj[2], ri_q[2], rj_q[2];
        int rtmp, rnew;
        r(ri, i, ri_q);
        r(rj, j, rj_q);

        // [0, 0]
        rtmp = inner2(ri, rj);

        // Different X quarter
        if (ri_q[0] != rj_q[0]) {
            if (ri_q[0] == 0)
                ri[0] += w * 2;
            else
                rj[0] += w * 2;
            rnew = inner2(ri, rj);
            if (rtmp > rnew)
                rtmp = rnew;
        }

        // Different Y quarter
        if (ri_q[1] != rj_q[1]) {
            if (ri_q[1] == 0)
                ri[1] += l * 2;
            else
                rj[1] += l * 2;
            rnew = inner2(ri, rj);
            if (rtmp > rnew)
                rtmp = rnew;

            // Correction if both X and Y quarter are different
            if (ri_q[0] != rj_q[0]) {
                if (ri_q[0] == 0)
                    ri[0] -= w * 2;
                else
                    rj[0] -= w * 2;
                rnew = inner2(ri, rj);
                if (rtmp > rnew)
                    rtmp = rnew;
            }
        }

        // Update LUT
        rmin[i][j] = rtmp;
        rmin[j][i] = rtmp;

        return rtmp;
    }
}


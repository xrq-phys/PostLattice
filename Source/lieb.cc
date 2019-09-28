#include "lieb.hh"

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
    if (i % ncell)
        return 0; // Sorry for temporarily using 0-0 correlation as trash can.

    int icell = i / ncell,
        jcell = j / ncell,
        jsite = j % ncell;
    int xicl = icell % w,
        yicl = icell / w,
        xjcl = jcell % w,
        yjcl = jcell / w;
    int xrcl = (xjcl - xicl + w) % w,
        yrcl = (yjcl - yicl + l) % l;

    return (yrcl * w + xrcl) * ncell + jsite;
}

int lattice::lieb::calc_rmin(int i, int j) {
    if (i != 0)
        return -1;
    int r[2], r_q[2];

    this->r(r, j, r_q);
    if (r_q[0]) r[0] -= w * 2;
    if (r_q[1]) r[1] -= l * 2;

    return r[0]*r[0] + r[1]*r[1];
}


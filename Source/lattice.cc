#include "lattice.hh"
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>

// {
// Local Tool Functions.
int inner2(int *a, int*b)
{ return ((a[0] - b[0]) * (a[0] - b[0])) + ((a[1] - b[1]) * (a[1] - b[1])); }
// }

// {
// Virtual methods that MUST be overriden.
void lattice::lattice::r(double *r, int i)
{ abort(); }
void lattice::lattice::r(int *r, int i, int *r_q)
{ abort(); }
int lattice::lattice::calc_rmin(int i, int j)
{ abort(); }
// }

void lattice::sorted_rc(lattice::lattice &system)
{
    std::vector<int> r_v(system.n, 0);
    std::vector<int> n_v(system.n, 1);
    for (int i = 0; i < system.n; i++)
        r_v[i] = system.calc_rmin(0, i);
    std::sort(r_v.begin(), r_v.end());
    for (int i = 1; i != r_v.size(); )
        if (r_v[i] == r_v[i - 1]) {
            r_v.erase(r_v.begin() + i);
            n_v.erase(n_v.begin() + i);
            n_v[i - 1]++;
        } else
            i++;
    system.r_c = new int[r_v.size()];
    system.r_n = new int[n_v.size()];
    system.rc_n = r_v.size();
    for (int i = 0; i < r_v.size(); i++) {
        system.r_c[i] = r_v[i];
        system.r_n[i] = n_v[i];
    }
}

lattice::square2d::square2d(int w_i, int l_i)
: w(w_i), l(l_i), lattice::lattice(w_i * l_i, 2)
{
    sorted_rc(*this);

    for (int xi = 0; xi < w; xi++)
        for (int yi = 0; yi < l; yi++) {
            nn[idx_2d(xi, yi)][idx_2d(xi, (yi + 1)     % l)] = 1;
            nn[idx_2d(xi, yi)][idx_2d(xi, (yi - 1 + l) % l)] =-1;
            nn[idx_2d(xi, yi)][idx_2d((xi + 1)     % w, yi)] = 2;
            nn[idx_2d(xi, yi)][idx_2d((xi - 1 + w) % w, yi)] =-2;
        }
}

void lattice::square2d::r(int *r, int i, int *r_q)
{ 
    // This routine sets not only \vec r, but also the quarter r_i is in
    r[0] = i % w; r[1] = i / w; 
    if (r_q != nullptr) {
        if (r[0] < w / 2)
            r_q[0] = 0;
        else 
            r_q[0] = 1;
        if (r[1] < l / 2)
            r_q[1] = 0;
        else 
            r_q[1] = 1;
    }
}

void lattice::square2d::r(double *rd, int i)
{
    int ri[2];
    r(ri, i);
    rd[0] = double(ri[0]) / w;
    rd[1] = double(ri[1]) / l;
}

int lattice::square2d::calc_rmin(int i, int j)
{
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
                ri[0] += w;
            else
                rj[0] += w;
            rnew = inner2(ri, rj);
            if (rtmp > rnew)
                rtmp = rnew;
        }

        // Different Y quarter
        if (ri_q[1] != rj_q[1]) {
            if (ri_q[1] == 0)
                ri[1] += l;
            else 
                rj[1] += l;
            rnew = inner2(ri, rj);
            if (rtmp > rnew)
                rtmp = rnew;

            // Correction if both X and Y quarter are different
            if (ri_q[0] != rj_q[0]) {
                if (ri_q[0] == 0)
                    ri[0] -= w;
                else 
                    rj[0] -= w;
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
#pragma once 
#include <cmath>
#include <cstdio>

struct lattice
{
    int w;
    int n;
    int **rmin;
    int **nn;

    lattice(int w_i)
    : w(w_i)
    {
        n = w * w;
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

        for (int xi = 0; xi < w; xi++)
            for (int yi = 0; yi < w; yi++) {
                nn[idx_2d(xi, yi)][idx_2d(xi, (yi + 1)     % w)] = 1;
                nn[idx_2d(xi, yi)][idx_2d(xi, (yi - 1 + w) % w)] = 1;
                nn[idx_2d(xi, yi)][idx_2d((xi + 1)     % w, yi)] = 2;
                nn[idx_2d(xi, yi)][idx_2d((xi - 1 + w) % w, yi)] = 2;
            }
    }

    ~lattice()
    {
        for (int i = 0; i < n; i++) {
            delete[] nn[i];
            delete[] rmin[i];
        }
        delete[] nn;
        delete[] rmin;
    }

    void r(int *r, int i, int *r_q = nullptr) 
    { 
        // This routine sets not only \vec r, but also the quarter r_i is in
        r[0] = i % w; r[1] = i / w; 
        if (r_q != nullptr) {
            if (r[0] < w / 2)
                r_q[0] = 0;
            else 
                r_q[0] = 1;
            if (r[1] < w / 2)
                r_q[1] = 0;
            else 
                r_q[1] = 1;
        }
    }

    int idx_2d(int x, int y) 
    { return y*w + x; }

    int inner2(int *a, int*b)
    { return ((a[0] - b[0]) * (a[0] - b[0])) + ((a[1] - b[1]) * (a[1] - b[1])); }

    int calc_rmin(int i, int j)
    {
        if (rmin[i][j] > 0)
            return rmin[i][j];
        else {
            int ri[2], rj[2], ri_q[2], rj_q[2];
            int rtmp, rnew;
            r(ri, i, ri_q);
            r(rj, j, rj_q);

            // hard-code...
            // [todo] do optimization
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
                    ri[1] += w;
                else 
                    rj[1] += w;
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
};

struct sc_corr
{
    lattice &system;
    int run_mode;
    int rc_count;
    int rc_lst[12];
    int rx_lst[12];
    int ry_lst[12];
    double *values;

    sc_corr(lattice &system_i, int rc_count_i, int run_mode_i = false)
    : system(system_i), 
      rc_lst { 0, 1, 1 + 1, 4, 1 + 4, 4 + 4, 9, 9 + 1, 9 + 4, 16, 9 + 9, 16 + 1 },
      rx_lst { 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 3, 4 },
      ry_lst { 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 3, 1 }
    {
        rc_count = rc_count_i;
        run_mode = run_mode_i;
        values = new double[rc_count];
        for (int i = 0; i < rc_count; i++)
            values[i] = 0;
    }

    ~sc_corr()
    { delete[] values; }

    void measure(int a, int sa, int b, int sb, int i, int si, int j, int sj, double x)
    {
        int sign;
        int rmin;
        int rb[2], 
            ri[2];
        
        if (validate(a, sa, b, sb, i, si, j, sj)) {
            sign = 1;

            if (sb != si) 
                sign = -sign;
            if (system.nn[a][b] != system.nn[i][j])
                sign = -sign;

            if (run_mode) {
                system.r(rb, b);
                system.r(ri, i);
                for (int ii = 0; ii < rc_count; ii++)
                    if (ri[0] - rb[0] + system.w % system.w == rx_lst[ii] &&
                        ri[1] - rb[1] + system.w % system.w == ry_lst[ii]) {
#pragma omp critical
                        values[ii] += x * sign / (2. * system.n);
                        printf("%d %s %d %s %d %s %d %s %d %d %f\n",
                               a, sa ? "↓" : "↑", b, sb ? "↓" : "↑",
                               i, si ? "↓" : "↑", j, sj ? "↓" : "↑",
                               rx_lst[ii], ry_lst[ii], x);
                    }
            } else {
                rmin = system.calc_rmin(b, i);
                for (int ii = 0; ii < rc_count; ii++)
                    if (rmin == rc_lst[ii])
#pragma omp critical
                        values[ii] += x * sign / (2. * system.n);
            }
        }
    }

    bool validate(int a, int sa, int b, int sb, int i, int si, int j, int sj)
    { return si != sj && sa != sb && system.nn[a][b] && system.nn[i][j]; }
};

#include "operators.hh"
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <ccomplex>
#include <cmath>

// Pauli matrices.
static const int pauli_x[2][2] = { 0,  1, 1,  0 };
static const int pauli_y[2][2] = { 0, -1, 1,  0 };
static const int pauli_z[2][2] = { 1,  0, 0, -1 };

// {
// These are Tool Functions.
// TODO: Put it somewhere else.
inline double inner_r_qidx(double *r, const int *qidx, const int ndim)
{
    double sum = 0;
    for (int i = 0; i < ndim; i++)
        sum += 2 * M_PI * qidx[i] * r[i];
    return sum;
}

inline int idx_nd(const int *iqpt, const int *nqpt, const int ndim)
{
    int idx = 0, nelem = 1;
    for (int i = ndim - 1; i >=0; i--) {
        idx += iqpt[i] * nelem;
        nelem *= nqpt[i];
    }
    return idx;
}

void allocate_qpoints(int **qpts, const int *nqpt, int *iqpt, const int dim, const int ndim)
{
    if (dim < ndim) {
        for (int i = 0; i < nqpt[dim]; i++) {
            iqpt[dim] = i;
            allocate_qpoints(qpts, nqpt, iqpt, dim + 1, ndim);
        }
    } else {
        int idx = idx_nd(iqpt, nqpt, ndim);
        qpts[idx] = new int[ndim];
        std::memcpy(qpts[idx], iqpt, ndim * sizeof(int));
    }
}
// }

// {
// Virtual methods that MUST be overridden.
void operators::operators::measure(int a, int sa, int b, int sb, 
                                   int i, int si, int j, int sj, double x)
{ abort(); }
// }

void operators::doublon::measure(int i, int si, int j, int sj,
                                 int k, int sk, int l, int sl, double x)
{ values[0] += (i == j && k == l && i == k && si == sj && sk == sl && si != sk) ? x / 2. : 0; }

void operators::sc_corr::measure(int a, int sa, int b, int sb, 
                                 int i, int si, int j, int sj, double x)
{
    int sign;
    int rmin;
    int rb[2], 
        ri[2];
    
    if (validate(a, sa, b, sb, i, si, j, sj)) {
        sign = 1;

        if (sb != si) 
            sign = -sign;

        // Waveform factor only affects this part.
        if (form == 's' /* S */) { /* Do nothing. */ }
        else if (form == 'd' /* D_{XY} */) {
            int alpha_c = std::abs(system.nn[a][b]),
                alpha_a = std::abs(system.nn[i][j]);
            if (alpha_a != alpha_c &&
                alpha_a +  alpha_c == 3)
                sign = -sign;
        } else if (form < '4' && form > '0' /* P in 3 directions */) {
            if (std::abs(system.nn[a][b]) != std::abs(system.nn[i][j]) ||
                std::abs(system.nn[a][b]) != form - '0')
                return; // sign = 0;
            sign *= system.nn[a][b] == system.nn[i][j] ? 1 : -1;
        } else if (form == 'f' /* F_{X^3-3XY^2} */) {
            sign *= system.nn[a][b] > 0 ? (system.nn[a][b] % 2) * 2 - 1
                                        : (system.nn[a][b] % 2) * 2 + 1;
            sign *= system.nn[i][j] > 0 ? (system.nn[i][j] % 2) * 2 - 1
                                        : (system.nn[i][j] % 2) * 2 + 1;
        } else
            return; // Waveform not supported.

        rmin = system.calc_rmin(b, i);
        for (int ii = 0; ii < rc_count; ii++)
            if (rmin == system.r_c[ii])
#pragma omp critical
            values[ii] += x * sign / (2. * system.n);
    }
}

operators::spin_struct::spin_struct(lattice::lattice &system_i, const int *ndiv)
: operators::operators(), system(system_i)
{
    // Total number of q-points.
    n_points = 1;
    for (int i = 0; i < system.dim; i++) 
        n_points *= ndiv[i];

    values = new double[n_points];
    points = new int*[n_points];
    int qpt_buff[system.dim];
    allocate_qpoints(points, ndiv, qpt_buff, 0, system.dim);
}

void operators::spin_struct::measure(int ri, int si, int rj, int sj,
                                     int rk, int sk, int rl, int sl, double x)
{
    if (ri != rj || rk != rl)
        return;
    
    // Si \cdot Sk.
    double spin_part = (pauli_x[si][sj] * pauli_x[sk][sl] 
                      + pauli_z[si][sj] * pauli_z[sk][sl] 
                      - pauli_y[si][sj] * pauli_y[sk][sl]) / 4.;

    // exp(r * q) == cos(r * q) (by symmetry).
    double *dr = new double[system.dim],
           *r1 = new double[system.dim],
           *r2 = new double[system.dim];
    system.r(r1, ri);
    system.r(r2, rk);
    for (int i = 0; i < system.dim; i++)
        dr[i] = r1[i] - r2[i];
    for (int i = 0; i < n_points; i++) {
        double space_part = 1 / (3. * system.n * system.n) *
                            std::cos(inner_r_qidx(dr, points[i], system.dim));
#pragma omp critical
        values[i] += spin_part * space_part * x;
    }
}

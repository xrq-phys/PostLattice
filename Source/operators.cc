#define _USE_MATH_DEFINES
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
    for (int i = ndim - 1; i >= 0; i--) {
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

void mini_fourier(int cdim, int *ndiv, int *nidxspc, int *qidx,
                  lattice::lattice &system, double *val_mat, double *values) {
    if (cdim != system.dim)
        for (int i = 0; i < ndiv[cdim]; i++) {
            qidx[cdim] = i;
            mini_fourier(cdim + 1, ndiv, nidxspc, qidx, system, val_mat, values);
        } else {
            // exp(r * q) == cos(r * q) (by symmetry).
            double *dr = new double[system.dim],
                   *r1 = new double[system.dim],
                   *r2 = new double[system.dim];
            for (int i = 0; i < system.n; i++)
                for (int j = 0; j < system.ncell; j++) {
                    int idx = 0;
                    system.r(r1, i);
                    system.r(r2, j);
                    for (int k = 0; k < system.dim; k++)
                        dr[i] = r1[i] - r2[i];
                    for (int k = 1; k < system.dim; k++)
                        idx += qidx[k] * nidxspc[k];
                    values[idx] += 1 / (3. * system.n * system.n) *
                                   std::cos(inner_r_qidx(dr, qidx, system.dim)) *
                                   val_mat[system.idx_rij(i, j)];
                }
            delete[] dr; delete[] r1; delete[] r2;
        }
}

void stat_rc(double *values, lattice::lattice &system, int rc_count, char mode, double *val_mat)
{
    // Compute according specified modes.
    if (mode == 'M' || mode == 'm')
        for (int i = 0; i < system.n * system.ncell; i++) {
            int rmin = system.calc_rmin(i / system.n /* from nth site in 0th cell */,
                                        i % system.n /* to nth site in whole lattice */);
            for (int j = 0; j < rc_count; j++)
                if (rmin == system.r_c[j])
                    if (std::abs(values[j]) < std::abs(val_mat[i]))
                        values[j] = val_mat[i];
        }
    else
        for (int i = 0; i < system.n * system.ncell; i++) {
            int rmin = system.calc_rmin(i / system.n, i % system.n);
            for (int j = 0; j < rc_count; j++)
                if (rmin == system.r_c[j])
                    switch (mode) {
                    case 'A':
                    case 'a':
                        if (system.ncell > 1)
                            values[j] += val_mat[i] * pow(-1, i / system.n);
                        else
                            if (system.dim == 2) {
                                int r[2]; system.r(r, i);
                                values[j] += val_mat[i] * pow(-1, r[0] + r[1]);
                            }
                        break;
                    default:
                        values[j] += abs(val_mat[i]);
                    }
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

        val_mat[system.idx_rij(i, b)] += x * sign / (2. * system.n);
    }
}

void operators::sc_corr::refresh()
{ refresh('S'); }

void operators::sc_corr::refresh(char mode)
{ stat_rc(values, system, rc_count, mode, val_mat); }

operators::site_corr::site_corr(lattice::lattice &system_i)
: operators::operators(), system(system_i)
{
    int idx;
    // Initialize points
    val_mat = new double[system_i.n * system_i.ncell];
    connection = new int*[system_i.n * system_i.ncell];
    conn_slave = 0;
    for (int i = 0; i < system_i.n * system_i.ncell; i++) {
        val_mat[i] = 0;
        connection[i] = new int[2];
    }
    for (int i = 0; i < system_i.n; i++)
        for (int j = 0; j < system_i.ncell; j++) {
            idx = system_i.idx_rij(i, j);
            connection[idx][0] = i;
            connection[idx][1] = j;
        }
}

operators::site_corr::site_corr(site_corr &opr_i)
: operators::operators(), system(opr_i.system)
{
    lattice::lattice &system_i = opr_i.system;
    val_mat = new double[system_i.n * system_i.ncell];
    connection = opr_i.connection;
    conn_slave = 1;
    for (int i = 0; i < system_i.n * system_i.ncell; i++)
        val_mat[i] = 0;
}

operators::site_corr::~site_corr()
{
    if (!conn_slave) {
        for (int i = 0; i < system.n * system.ncell; i++)
            delete[] connection[i];
        delete[] connection;
    }
    delete[] val_mat;
}

void operators::site_corr::refresh(int *ndiv) {
    int n_qpts = 1;
    int *nidxspc = new int[system.dim];
    int *qidx    = new int[system.dim];
    // Allocate memory.
    for (int i = system.dim - 1; i >= 0; i--) {
        n_qpts *= ndiv[i];
        nidxspc[i] = n_qpts;
    }
    delete[] values; values = new double[n_qpts];
    for (int i = 0; i < n_qpts; i++) values[i] = 0;

    // Execute
    mini_fourier(0, ndiv, nidxspc, qidx, system, val_mat, values);

    delete[] nidxspc;
    delete[] qidx;
}

void operators::site_corr::refresh(char mode, int rc_count)
{
    delete[] values; values = new double[rc_count];
    for (int i = 0; i < rc_count; i++) values[i] = 0;
    stat_rc(values, system, rc_count, mode, val_mat);
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

    val_mat[system.idx_rij(ri, rk)] += spin_part * x;
}

void operators::charge_struct::measure(int ri, int si, int rj, int sj,
                                       int rk, int sk, int rl, int sl, double x)
{
    if (ri != rj || rk != rl ||
        si != sj || sk != sl)
        return;

    // Regardless spin direction, add up.
    val_mat[system.idx_rij(ri, rk)] += x;
}

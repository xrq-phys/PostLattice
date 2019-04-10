#define _USE_MATH_DEFINES
#include "operators.hh"
#include <sstream>
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
                                   int i, int si, int j, int sj, std::complex<double> x)
{ abort(); }
// }

void operators::doublon::measure(int i, int si, int j, int sj,
                                 int k, int sk, int l, int sl, std::complex<double> x)
{ values[0] += (i == j && k == l && i == k && si == sj && sk == sl && si != sk) ? x / 2. : 0; }

operators::sc_corr::sc_corr(lattice::lattice &system_i, const int rc_count_i, const int use_p_i, std::string form_s)
: operators::operators(), system(system_i), use_p(use_p_i),
  rc_count(rc_count_i < system_i.rc_n ? rc_count_i : system_i.rc_n)
{
    double xre, xim;
    char form_name;
    // R-based correlation.
    values = new std::complex<double>[rc_count];
    for (int i = 0; i < rc_count; i++)
        values[i] = 0;

    // Real-space correlation.
    val_mat = new std::complex<double>[system_i.n * system_i.ncell];
    for (int i = 0; i < system_i.n * system_i.ncell; i++)
        val_mat[i] = 0;

    // Waveform.
    form = new std::complex<double>[system_i.nbond];
    for (int i = 0; i < system_i.nbond; i++)
        form[i] = 0;
    std::stringstream form_ss(form_s);
    form_ss >> form_name;
    if (form_name == 's') {
        for (int i = 0; i < system_i.nbond; i++)
            form[i] = 1;
    } else if (form_name == 'd') {
        form[0] = 1;
        if (system_i.nbond > 1)
            form[1] = 1;
        for (int i = 2; i < system_i.nbond; i++)
            form[i] = -1;
    } else if (form_name < '9' && form_name > '0' /* P in N directions */) {
        int bond_sel = form_name - '1';
        form[bond_sel * 2] = 1;
        form[bond_sel * 2 + 1] = -1;
    } else if (form_name == 'f') {
        for (int i = 0; i < 3; i++) {
            form[2 * i] = 1;
            form[2 * i + 1] = -1;
        }
    } else if (form_name == '+' /* Customized wave form, see document. */) {
        for (int i = 0; i < system_i.nbond; i++) {
            form_ss >> xre >> xim;
            form[i].real(xre);
            form[i].imag(xim);
        }
    } else { /* Not supported. */ }
}

void operators::sc_corr::measure(int a, int sa, int b, int sb, 
                                 int i, int si, int j, int sj, std::complex<double> x)
{
    std::complex<double> form_loc;
    int rmin;
    
    if (validate(a, sa, b, sb, i, si, j, sj)) {
        form_loc = 1;

        if (sb != si) 
            form_loc *= -1;

        // Waveform factor only affects this part.
        int alpha_c = (std::abs(system.nn[a][b]) - 1) * 2 + (system.nn[a][b] > 0 ? 0 : 1);
        int alpha_a = (std::abs(system.nn[i][j]) - 1) * 2 + (system.nn[i][j] > 0 ? 0 : 1);
        form_loc *= form[alpha_a] * form[alpha_c];

        val_mat[system.idx_rij(i, b)] += x * form_loc / (2. * system.n);
    }
}

void operators::sc_corr::refresh()
{ refresh('S'); }

void operators::sc_corr::refresh(char mode)
{
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
                    values[j] += val_mat[i];
        }
}

operators::spin_struct::spin_struct(lattice::lattice &system_i, const int *ndiv)
: operators::operators(), system(system_i)
{
    // Total number of q-points.
    n_points = 1;
    for (int i = 0; i < system.dim; i++) 
        n_points *= ndiv[i];

    values = new std::complex<double>[n_points];
    for (int i = 0; i < n_points; i++)
        values[i] = 0;
    points = new int*[n_points];
    int *qpt_buff = new int[system.dim];
    allocate_qpoints(points, ndiv, qpt_buff, 0, system.dim);
    delete[] qpt_buff;
}

void operators::spin_struct::measure(int ri, int si, int rj, int sj,
                                     int rk, int sk, int rl, int sl, std::complex<double> x)
{
    if (ri != rj || rk != rl)
        return;
    
    // Si \cdot Sk.
    double spin_part = (pauli_x[si][sj] * pauli_x[sk][sl] 
                      + pauli_z[si][sj] * pauli_z[sk][sl] 
                      - pauli_y[si][sj] * pauli_y[sk][sl]) / 4.;

    double *dr = new double[system.dim],
           *r1 = new double[system.dim],
           *r2 = new double[system.dim];
    system.r(r1, ri);
    system.r(r2, rk);
    for (int i = 0; i < system.dim; i++)
        dr[i] = r1[i] - r2[i];
    for (int i = 0; i < n_points; i++) {
        std::complex<double> space_part = 1 / (3. * system.n * system.n) *
            std::exp(std::complex<double>(0, 1) * inner_r_qidx(dr, points[i], system.dim));
        values[i] += spin_part * space_part * x;
    }
}

#include "measure.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <cmath>
#include <cassert>
#include <omp.h>

using namespace std;

// {
// TOOL ROUTINES
// TODO: Put it somewhere else.
long idx_4so(int i, int si, int j, int sj, int a, int sa, int b, int sb, int n)
{ return i * n * n * n * 16 + j * n * n * 16 + a * n * 16 + b * 16 + si * 8 + sj * 4 + sa * 2 + sb; }
long idx_1so(int i, int si)
{ return i * 2 + si; }
// }

void measure(lattice::lattice &physics, operator_options &options)
{
    fstream fid_g1e = fstream(options.g1e_fnm);
    fstream fid_g2e = fstream(options.g2e_fnm);
    assert(!(fid_g1e.fail() || fid_g2e.fail()));

    // Check & modify some options.
    if (options.af_n == nullptr) {
        options.af_n = new int[physics.dim];
        for (int i = 0; i < physics.dim; i++)
            options.af_n[i] = 1;
    }
    operators::doublon opr_db;
    operators::spin_struct opr_af(physics, options.af_n);
    operators::sc_corr opr_sc(physics, options.sc_n, options.sc);

    map <long, int> seen_ijab, seen_iiii;
    int *ri, *si, *ra, *sa, *rj, *sj, *rb, *sb;
    int nproc, iproc;
    int iline;
    bool *exec_this;
    double *x, **x1e, dummy;

    nproc = omp_get_num_threads();
    ri = new int[nproc]; si = new int[nproc];
    rj = new int[nproc]; sj = new int[nproc];
    ra = new int[nproc]; sa = new int[nproc];
    rb = new int[nproc]; sb = new int[nproc];
    x = new double[nproc];
    exec_this = new bool[nproc];
    x1e = new double *[physics.n * 2];
    for (int i = 0; i < physics.n * 2; i++) {
        x1e[i] = new double[physics.n * 2];
        for (int j = 0; j < physics.n * 2; j++)
            x1e[i][j] = 0.;
    }

    // Load in one-body part.
    while (true) {
        fid_g1e >> ra[0] >> sa[0] >> ri[0] >> si[0] >> x[0] >> dummy;
        if (fid_g1e.eof())
            break;
        x1e[ri[0] * 2 + si[0]][ra[0] * 2 + sa[0]] = x[0];
    }
    fid_g1e.close();

    iline = 0;
    while (!fid_g2e.eof()) {
        for (int i = 0; i < nproc; i++) {
            fid_g2e >> rb[i] >> sb[i] >> rj[i] >> sj[i]
                    >> ra[i] >> sa[i] >> ri[i] >> si[i] >> x[i] >> dummy;
            iline += 1;
            if (fid_g2e.eof()) {
                for (int j = i; j < nproc; j++)
                    exec_this[j] = false;
                break;
            }
            long fld4 = idx_4so(ri[i], si[i], rj[i], sj[i], ra[i], sa[i], rb[i], sb[i], physics.n);
            if (seen_ijab.count(fld4)) {
                if (options.verbose)
                    cout << "Duplicate found at " << iline << " with " << seen_ijab[fld4];
                exec_this[i] = false;
                break;
            }
            seen_ijab[fld4] = iline;
            exec_this[i] = true;
        }

#pragma omp parallel for default(shared) private(iproc, dummy)
        for (int i = 0; i < nproc; i++)
            if (exec_this[i] && sa[i] == si[i]) {
#pragma omp critical
                if (options.verbose)
                    cout << ' ' << ri[i] << ' ' << si[i] << ' ' << rj[i] << ' ' << sj[i]
                         << ' ' << ra[i] << ' ' << sa[i] << ' ' << ra[i] << ' ' << sb[i] << ' ' << endl;
                if (options.sc != '-') {
                    dummy = x[i] + int(ra[i] == ri[i] && rb[i] == rj[i]) / 2.
                            - x1e[rj[i] * 2 + sj[i]][rb[i] * 2 + sb[i]] * int(ra[i] == ri[i]) / 2.
                            - x1e[ri[i] * 2 + si[i]][ra[i] * 2 + sa[i]] * int(rb[i] == rj[i]) / 2.;
                    opr_sc.measure(rb[i], sb[i], ra[i], sa[i], rj[i], sj[i], ri[i], si[i], -dummy);
                    opr_sc.measure(rb[i], sb[i], ra[i], sa[i], ri[i], si[i], rj[i], sj[i],  dummy);
                }
                if (options.af) {
                    opr_af.measure(rb[i], sb[i], rj[i], sj[i], ra[i], sa[i], ri[i], si[i], x[i]);
                    if (si[i] != sj[i])
                        opr_af.measure(rb[i], sb[i], ri[i], si[i], ra[i], sa[i], rj[i], sj[i],
                                       -x[i] + (ri[i] == ra[i] ? x1e[rb[i] + sb[i]][rj[i] + sj[i]] : 0));
                }
                if (options.db)
                    opr_db.measure(rb[i], sb[i], rj[i], sj[i], ra[i], sa[i], ri[i], si[i], x[i]);
            }
    }
    fid_g2e.close();

    fstream fid_out_sc(options.sc_fnm, fstream::out);
    for (int i = 0; i < opr_sc.rc_count; i++)
        fid_out_sc << setw(8) << sqrt(double(physics.r_c[i])) << ' ' << opr_sc.values[i] << endl;
    fid_out_sc.close();

    fstream fid_out_af(options.af_fnm, fstream::out);
    for (int i = 0; i < opr_af.n_points; i++)
        fid_out_af << opr_af.points[i][0] << ' '
                   << opr_af.points[i][1] << ' ' << opr_af.values[i] << endl;
    fid_out_af.close();

    fstream fid_out_db(options.db_fnm, fstream::out);
    fid_out_db << opr_db.values[0] << endl;
    fid_out_db.close();

    delete[] ri; delete[] si; delete[] rj; delete[] sj;
    delete[] ra; delete[] sa; delete[] rb; delete[] sb;
    delete[] x;  delete[] exec_this;
    for (int i = 0; i < physics.n * 2; i++)
        delete[] x1e[i];
    delete[] x1e;
}

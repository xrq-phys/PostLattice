#include "measure.hh"
#include "plot.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <cmath>
#include <cassert>

using namespace std;

// {
// TOOL ROUTINES
// TODO: Put it somewhere else.
long idx_4so(int i, int si, int j, int sj, int a, int sa, int b, int sb, int n)
{ return i * n * n * n * 16 + j * n * n * 16 + a * n * 16 + b * 16 + si * 8 + sj * 4 + sa * 2 + sb; }
long idx_1so(int i, int si)
{ return i * 2 + si; }
// }

void convert(lattice::lattice &physics, operator_options &options)
{
    fstream fid_mat(options.af_fnm);
    operators::site_corr opr(physics);

    while (!fid_mat.eof()) {
        int i, j;
        double x;

        fid_mat >> i >> j >> x;
        if (fid_mat.eof())
            break;
        opr.val_mat[physics.idx_rij(i, j)] = x;
    }
    fid_mat.close();

    fstream fid_rc(options.af_rc_fnm);
    opr.refresh(options.rc_af_stat, physics.rc_n);
    for (int i = 0; i < physics.rc_n; i++)
        fid_rc << setw( 8) << fixed << sqrt(double(physics.r_c[i])) << ' '
               << setw(18) << scientific << opr.values[i] << ' '
               << setw( 8) << (options.rc_af_stat == 'M' ||
                               options.rc_af_stat == 'm' ? 1 : physics.r_n[i]) << endl;
    fid_rc.close();
}

void measure(lattice::lattice &physics, operator_options &options)
{
    fstream fid_g1e(options.g1e_fnm);
    fstream fid_g2e(options.g2e_fnm);
    fstream fid_g1e_idx(options.g1e_idx_fnm);
    fstream fid_g2e_idx(options.g2e_idx_fnm);
    assert(!(fid_g1e.fail() || fid_g2e.fail()));

    // Check & modify some options.
    if (options.af_n == nullptr) {
        options.af_n = new int[physics.dim];
        for (int i = 0; i < physics.dim; i++)
            options.af_n[i] = 1;
    }
    operators::doublon opr_db;
    operators::spin_struct opr_af(physics);
    operators::charge_struct opr_st(opr_af);
    operators::sc_corr opr_sc(physics, options.sc_n, options.sc, options.sc_use_p);

    map <long, int> seen_ijab;
    int ri, si, ra, sa, rj, sj, rb, sb;
    int iline;
    double x, **x1e, dummy;

    x1e = new double *[physics.n * 2];
    for (int i = 0; i < physics.n * 2; i++) {
        x1e[i] = new double[physics.n * 2];
        for (int j = 0; j < physics.n * 2; j++)
            x1e[i][j] = 0.;
    }

    // Skip GF definition file headers.
    if (!options.ca_verbose) {
        char cbuffer[200];
        for (int i = 0; i < 5; i++) {
            fid_g1e_idx.getline(cbuffer, 200, '\n');
            fid_g2e_idx.getline(cbuffer, 200, '\n');
        }
    } // Yet another non-empty check.
    assert(!(fid_g1e.eof() || fid_g2e.eof()));

    // Load in one-body part.
    while (true) {
        if (options.ca_verbose)
            fid_g1e >> ra >> sa >> ri >> si >> x >> dummy;
        else {
            if (options.ca_legacy) {
                fid_g1e_idx >> iline >> ra >> ri >> sa;
                si = sa;
            } else
                fid_g1e_idx >> ra >> sa >> ri >> si;
            fid_g1e >> x;
        }
        if (fid_g1e.eof())
            break;
        x1e[ri * 2 + si][ra * 2 + sa] = x;
    }
    fid_g1e.close();

    iline = 0;
    while (!fid_g2e.eof()) {
        if (options.ca_verbose)
            fid_g2e >> rb >> sb >> rj >> sj
                    >> ra >> sa >> ri >> si >> x >> dummy;
        else {
            if (options.ca_legacy) {
                fid_g2e_idx >> rb >> rj >> sb >> ra >> ri >> sa;
                sj = sb;
                si = sa;
            } else
                fid_g2e_idx >> rb >> sb >> rj >> sj
                             >> ra >> sa >> ri >> si;
            fid_g2e >> x;
        }
        iline += 1;
        if (fid_g2e.eof())
            break;
        if (options.chk_duplicate) {
            long fld4 = idx_4so(ri, si, rj, sj, ra, sa, rb, sb, physics.n);
            if (seen_ijab.count(fld4)) {
                if (options.verbose)
                    cout << "Duplicate found at " << iline << " with " << seen_ijab[fld4] << endl;
                continue;
            }
            seen_ijab[fld4] = iline;
        }

        // if (options.verbose)
        //     cout << ' ' << ri << ' ' << si << ' ' << rj << ' ' << sj
        //          << ' ' << ra << ' ' << sa << ' ' << ra << ' ' << sb << ' ' << endl;
        if (options.sc != '-' && si == sa &&
            (options.sc_simple ? ra == 0 : 1)) { // From HPhi, we don't use DUUD/UDDU terms.
            if (!opr_sc.use_p) {
                opr_sc.measure(rb, sb, ra, sa, rj, sj, ri, si, -x);
                opr_sc.measure(rb, sb, ra, sa, ri, si, rj, sj,  x);
            } else
                opr_sc.measure(rb, sb, ra, sa, rj, sj, ri, si,
                               -x + (ra == rj ? x1e[ri * 2 + si][rb * 2 + sb] : 0));
        }
        if (options.af) {
            opr_af.measure(rb, sb, rj, sj, ra, sa, ri, si, x);
            if (si != sj && !options.direct) // For HPhi, use flipping terms instead.
                opr_af.measure(rb, sb, ri, si, ra, sa, rj, sj,
                               -x + (ri == ra ? x1e[rj * 2 + sj][rb * 2 + sb] : 0));
        }
        if (options.st)
            opr_st.measure(rb, sb, rj, sj, ra, sa, ri, si, x);
        if (options.db)
            opr_db.measure(rb, sb, rj, sj, ra, sa, ri, si, x);
    }
    fid_g2e.close();

    if (options.sc != '-') {
        opr_sc.refresh(options.rc_sc_stat);
        fstream fid_out_sc(options.sc_fnm, fstream::out);
        fstream fid_out_sc_rc(options.sc_rc_fnm, fstream::out);

        for (int k = 0; k < physics.n * physics.ncell; k++) {
            int j = k % physics.n,
                i = k / physics.n;
            if (options.corr_r) {
                if (i == 0) {
                    double r[2];
                    physics.r(r, j);
                    fid_out_sc << scientific << r[0] << ' ' << r[1] << ' '
                               << opr_sc.val_mat[k] << endl;
                }
            } else
                fid_out_sc << scientific << j << ' ' << i << ' '
                           << opr_sc.val_mat[k] << endl;
        }
        fid_out_sc.close();

        for (int i = 0; i < opr_sc.rc_count; i++)
            fid_out_sc_rc << setw( 8) << fixed << sqrt(double(physics.r_c[i])) << ' '
                          << setw(18) << scientific << opr_sc.values[i] * (options.sc_simple ? physics.n : 1) << ' '
                          << setw( 8) << (options.rc_sc_stat == 'M' ||
                                          options.rc_sc_stat == 'm' ? 1 : physics.r_n[i]) << endl;
        fid_out_sc_rc.close();
    }

    if (options.af) {
        opr_af.refresh(options.rc_af_stat, physics.rc_n); ///< [TODO] Stop using params for SC calc.
        fstream fid_out_af(options.af_fnm, fstream::out);
        fstream fid_out_af_rc(options.af_rc_fnm, fstream::out);

        for (int i = 0; i < physics.n * physics.ncell; i++) {
            if (options.corr_r) {
                if (opr_af.connection[i][1] == 0) {
                    double r[2];
                    physics.r(r, opr_af.connection[i][0]);
                    fid_out_af << scientific << r[0] << ' ' << r[1] << ' '
                               << opr_af.val_mat[i] << endl;
                }
            } else
                fid_out_af << scientific << opr_af.connection[i][0] << ' ' << opr_af.connection[i][1] << ' '
                                         << opr_af.val_mat[i] << endl;
        }
        fid_out_af.close();

        for (int i = 0; i < physics.rc_n; i++)
            fid_out_af_rc << setw( 8) << fixed << sqrt(double(physics.r_c[i])) << ' '
                          << setw(18) << scientific << opr_af.values[i] << ' '
                          << setw( 8) << (options.rc_af_stat == 'M' ||
                                          options.rc_af_stat == 'm' ? 1 : physics.r_n[i]) << endl;
        fid_out_af_rc.close();
    }

    if (options.st) {
        opr_st.refresh();
        fstream fid_out_st(options.st_fnm, fstream::out);
        for (int i = 0; i < physics.n * physics.ncell; i++) {
            if (options.corr_r) {
                if (opr_st.connection[i][1] == 0) {
                    double r[2];
                    physics.r(r, opr_st.connection[i][0]);
                    fid_out_st << scientific << r[0] << ' ' << r[1] << ' '
                               << opr_st.val_mat[i] << endl;
                }
            } else
                fid_out_st << scientific << opr_st.connection[i][0] << ' ' << opr_st.connection[i][1] << ' '
                           << opr_st.val_mat[i] << endl;
        }
        fid_out_st.close();
    }

    if (options.db) {
        opr_af.refresh();
        fstream fid_out_db(options.db_fnm, fstream::out);
        fid_out_db << scientific << opr_db.values[0] << endl;
        fid_out_db.close();
    }

    if (options.af && options.st)
        plot_lattice(options.mp_fnm.c_str(), physics, false, true, opr_st.val_mat, true, opr_af.val_mat);

    for (int i = 0; i < physics.n * 2; i++)
        delete[] x1e[i];
    delete[] x1e;
}

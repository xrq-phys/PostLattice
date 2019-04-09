#include "measure.hh"
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
    operators::spin_struct opr_af(physics);
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

    // Load in one-body part.
    while (true) {
        fid_g1e >> ra >> sa >> ri >> si >> x >> dummy;
        if (fid_g1e.eof())
            break;
        x1e[ri * 2 + si][ra * 2 + sa] = x;
    }
    fid_g1e.close();

    iline = 0;
    while (!fid_g2e.eof()) {
        fid_g2e >> rb >> sb >> rj >> sj
                >> ra >> sa >> ri >> si >> x >> dummy;
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
        if (options.db)
            opr_db.measure(rb, sb, rj, sj, ra, sa, ri, si, x);
    }
    fid_g2e.close();

    if (options.sc != '-') {
        opr_sc.refresh(options.sc_stat);
        fstream fid_out_sc(options.sc_fnm, fstream::out);
        for (int i = 0; i < opr_sc.rc_count; i++)
            fid_out_sc << setw( 8) << fixed << sqrt(double(physics.r_c[i])) << ' '
                       << setw(18) << scientific << opr_sc.values[i] * (options.sc_simple ? physics.n : 1) << ' '
                       << setw( 8) << (options.sc_stat == 'M' || 
                                       options.sc_stat == 'm' ? 1 : physics.r_n[i]) << endl;
        fid_out_sc.close();
    }

    if (options.af) {
        int n_tmp = 0;
        double v_tmp = 0.0;
        double v_err = 0.0;
        opr_af.refresh();
        fstream fid_out_af(options.af_fnm, fstream::out);
        for (int i = 0; i < physics.n * physics.ncell; i++) {
            fid_out_af << scientific << opr_af.connection[i][0] << ' ' << opr_af.connection[i][1] << ' '
                                     << opr_af.val_mat[i] << endl;
            if (physics.nn[opr_af.connection[i][0]][opr_af.connection[i][1]]) {
                n_tmp++;
                v_tmp += opr_af.val_mat[i]*4;
                v_err += pow(opr_af.val_mat[i]*4, 2);
            }
        }
        fid_out_af.close();
        v_tmp /= n_tmp;
        v_err /= n_tmp;
        printf("NN Spin Correlation: Cos(Theta) = %lf (Std_Err = %lf)\n", v_tmp, sqrt(v_err - v_tmp * v_tmp));
    }

    if (options.db) {
        opr_af.refresh();
        fstream fid_out_db(options.db_fnm, fstream::out);
        fid_out_db << scientific << opr_db.values[0] << endl;
        fid_out_db.close();
    }

    for (int i = 0; i < physics.n * 2; i++)
        delete[] x1e[i];
    delete[] x1e;
}

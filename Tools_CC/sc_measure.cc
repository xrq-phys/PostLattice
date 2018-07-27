// Only 2D Square Lattices Now, Sorry.
// NB IMPORTANT: all r lenghs are stored as r^2 to avoid usage of floating points.

#include <iomanip>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <map>
#include <omp.h>
#include "sc_classes.hh"

using namespace std;

long idx_4so(int i, int si, int j, int sj, int a, int sa, int b, int sb, int n)
{ return i * n * n * n * 16 + j * n * n * 16 + a * n * 16 + b * 16 + si * 8 + sj * 4 + sa * 2 + sb; }

long idx_1so(int i, int si)
{ return i * 2 + si; }

int main(const int argc, const char *argv[])
{
    assert(argc > 4);
    fstream fid_g1e(argv[2]);
    fstream fid_g2e(argv[3]);
    lattice physics(atoi(argv[1]));
    sc_corr quantity(physics, atoi(argv[4]), argc == 6);
    assert(!(fid_g1e.fail() || fid_g2e.fail()));
    map <long, int> seen_ijab, seen_iiii;
    int *ri, *si, *ra, *sa, *rj, *sj, *rb, *sb;
    int nproc, iproc;
    int iline;
    bool *exec_this;
    double *x, dummy;

    nproc = omp_get_num_threads();
    ri = new int[nproc];
    si = new int[nproc];
    rj = new int[nproc];
    sj = new int[nproc];
    ra = new int[nproc];
    sa = new int[nproc];
    rb = new int[nproc];
    sb = new int[nproc];
    x = new double[nproc];
    exec_this = new bool[nproc];

    // Contribution from 2-body.
    iline = 0;
    while (!fid_g2e.eof()) {
        for (int i = 0; i < nproc; i++) {
            fid_g2e >> rb[i] >> sb[i] >> rj[i] >> sj[i] >> ra[i] >> sa[i] >> ri[i] >> si[i] >> x[i] >> dummy;
            iline += 1;
            if (fid_g2e.eof()) {
                for (int j = i; j < nproc; j++)
                    exec_this[j] = false;
                break;
            }
            int fld1 = 0;
            int fld4 = idx_4so(ri[i], si[i], rj[i], sj[i], ra[i], sa[i], rb[i], sb[i], physics.n);
            bool use_fld1 = ri[i] == rj[i] && ri[i] == ra[i] && ri[i] == rb[i] && si[i] == sa[i] && si[i] != sj[i];
            if (use_fld1) fld1 = idx_1so(ri[i], si[i]);
            if (seen_ijab.count(fld4)) { 
                cout << "Duplicate found at " << iline << " with " << seen_ijab[fld4];
                if (use_fld1)
                    cout << " and " << seen_iiii[fld1];
                cout << endl;
                exec_this[i] = false; 
                break; 
            }
            if (use_fld1)
                if (!seen_iiii.count(fld1))
                    seen_iiii[fld1] = iline;
                else
                    seen_ijab[fld4] = iline;
            else
                seen_ijab[fld4] = iline;
            exec_this[i] = true;
        }

#pragma omp parallel for default(shared) private(iproc)
        for (int i = 0; i < nproc; i++)
            if (exec_this[i] && sa[i] == si[i]) {
                quantity.measure(rb[i], sb[i], ra[i], sa[i], rj[i], sj[i], ri[i], si[i], -x[i]);
                quantity.measure(rb[i], sb[i], ra[i], sa[i], ri[i], si[i], rj[i], sj[i],  x[i]);
            }
    }
    fid_g2e.close();

    /* Contribution from 1-body.
    while (!fid_g1e.eof()) {
        for (int i = 0; i < nproc; i++) {
            fid_g1e >> ra[i] >> sa[i] >> ri[i] >> si[i] >> x[i] >> dummy;
            if (fid_g1e.eof())
                for (int j = i; j < nproc; j++)
                    exec_this[j] = false;
        }
        
#pragma omp parallel for default(shared) private(iproc)
        for (int i = 0; i < nproc; i++)
            if (exec_this[i])
                for (int rk = 0; rk < physics.n; rk++)
                    for (int sk = 0; sk < 2; sk++)
                        quantity.measure(ra[i], sa[i], rk, sk, rk, sk, ri[i], si[i], x[i]);
    }*/
    fid_g1e.close();

    fstream fid_out("sc.txt", fstream::out);
    for (int i = 0; i < atoi(argv[4]); i++)
        fid_out << setw(8) << sqrt(double(quantity.rc_lst[i])) << ' ' << quantity.values[i] << endl;
    fid_out.close();

    return 0;
}

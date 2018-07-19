// Only 2D Square Lattices Now, Sorry.
// NB IMPORTANT: all r lenghs are stored as r^2 to avoid usage of floating points.

#include <iomanip>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <omp.h>
#include "sc_classes.hh"

using namespace std;

int main(const int argc, const char *argv[])
{
    assert(argc == 5);
    fstream fid_g1e(argv[2]);
    fstream fid_g2e(argv[3]);
    lattice physics(atoi(argv[1]));
    sc_corr quantity(physics, atoi(argv[4]));
    assert(!(fid_g1e.fail() || fid_g2e.fail()));
    double *x, dummy;
    int *ri, *si, *ra, *sa, *rj, *sj, *rb, *sb;
    int nproc, iproc;
    bool *exec_this;

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
    for (auto i = 0; i < nproc; i++)
        exec_this[i] = true;

    // Contribution from 2-body.
    while (!fid_g2e.eof()) {
        for (auto i = 0; i < nproc; i++) {
            fid_g2e >> rb[i] >> sb[i] >> rj[i] >> sj[i] >> ra[i] >> sa[i] >> ri[i] >> si[i] >> x[i] >> dummy;
            if (fid_g2e.eof())
                for (auto j = i; j < nproc; j++)
                    exec_this[j] = false;
        }
        
#pragma omp parallel for default(shared) private(iproc)
        for (auto i = 0; i < nproc; i++)
            if (exec_this[i])
                quantity.measure(rb[i], sb[i], ra[i], sa[i], rj[i], sj[i], ri[i], si[i], -x[i]);
    }
    fid_g2e.close();

    for (auto i = 0; i < nproc; i++)
        exec_this[i] = true;

    // Contribution from 1-body.
    while (!fid_g1e.eof()) {
        for (auto i = 0; i < nproc; i++) {
            fid_g1e >> ra[i] >> sa[i] >> ri[i] >> si[i] >> x[i] >> dummy;
            if (fid_g1e.eof())
                for (auto j = i; j < nproc; j++)
                    exec_this[j] = false;
        }
        
#pragma omp parallel for default(shared) private(iproc)
        for (auto i = 0; i < nproc; i++)
            if (exec_this[i])
                for (auto rk = 0; rk < physics.n; rk++)
                    for (auto sk = 0; sk < 2; sk++)
                        quantity.measure(ra[i], sa[i], rk, sk, rk, sk, ri[i], si[i], x[i]);
    }
    fid_g1e.close();

    fstream fid_out("sc.txt", fstream::out);
    for (auto i = 0; i < atoi(argv[4]); i++)
        fid_out << setw(8) << sqrt(double(quantity.rc_lst[i])) << ' ' << quantity.values[i] << endl;
    fid_out.close();

    return 0;
}

#include "plot.hh"
#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>

// Predefined image offset / scaling
#define x_offset 3
#define y_offset 3
#define x_scal 10
#define y_scal 10

using namespace std;

void plot_lattice(const char *mp_fnm, lattice::lattice &system, const bool label_on,
                  double (*size_func)(double *, int), double *size_mat,
                  double (*spin_func)(double *, int), double *spin_mat)
{
    assert(system.dim == 2);
    double r [2],
           rj[2];
    double dsize;
    fstream mp_fid(mp_fnm, fstream::out);
    // Metapost file header.
    mp_fid << "beginfig(1);" << endl;
    mp_fid << "pickup pencircle scaled 1pt;" << endl;

    // PBC bonds
    for (int i = 0; i < system.n; i++) {
        system.r(r, i);

        for (int j = i; j < system.n; j++)
            if (system.nn[i][j]) {
                system.r(rj, j);

                if (pow(r[0] - rj[0], 2) + pow(r[1] - rj[1], 2) >= 1. / 4)
                    mp_fid << "draw ("
                        << r [0] * x_scal + x_offset << "cm," 
                        << r [1] * y_scal + y_offset << "cm)--(" 
                        << rj[0] * x_scal + x_offset << "cm," 
                        << rj[1] * y_scal + y_offset << "cm)"
                        << " withcolor (1.,.7,.0);" << endl;
            }
    }

    // In-supercell bonds
    for (int i = 0; i < system.n; i++) {
        system.r(r, i);

        for (int j = i; j < system.n; j++)
            if (system.nn[i][j]) {
                system.r(rj, j);

                if (pow(r[0] - rj[0], 2) + pow(r[1] - rj[1], 2) < 1. / 4)
                    mp_fid << "draw ("
                        << r [0] * x_scal + x_offset << "cm," 
                        << r [1] * y_scal + y_offset << "cm)--(" 
                        << rj[0] * x_scal + x_offset << "cm," 
                        << rj[1] * y_scal + y_offset << "cm);" << endl;
            }
    }

    // Site
    for (int i = 0; i < system.n; i++) {
        system.r(r, i);
        // Label
        if (label_on)
            mp_fid << "label.bot(btex" << i << "," << system.calc_rmin(0, i) << "etex,("
                   << r[0] * x_scal + x_offset << "cm,"
                   << r[1] * y_scal + y_offset << "cm));" << endl;
        // Lattice point or density correlation
        dsize = size_func(size_mat, i);
        mp_fid << "pickup pencircle scaled " << dsize << "pt;" << endl << "drawdot ("
               << r[0] * x_scal + x_offset << "cm,"
               << r[1] * y_scal + y_offset << "cm) withcolor (.0,.0,1.);" << endl;
        // Spin correlation
        dsize = spin_func(spin_mat, i);
        mp_fid << "pickup pencircle scaled 1pt;" << endl;
        if (abs(dsize) < 1) {
            dsize *= 20 * system.n;
            mp_fid << "drawarrow ("
                   << r[0] * x_scal + x_offset - dsize << "cm,"
                   << r[1] * y_scal + y_offset + dsize << "cm) withcolor ("
                   << int(dsize > 0) << ".,.0," << int(dsize < 0) << ".);" << endl;
        }
    }

    // Metapost file ending.
    mp_fid << "endfig;" << endl 
           << "end"     << endl;
}
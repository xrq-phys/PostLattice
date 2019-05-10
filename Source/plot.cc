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

double find_min(const double *arr, int n)
{ double m = 0; for (int i = 0; i < n; i++) if (arr[i] > m) m = arr[i]; return m; }

void plot_lattice(const char *mp_fnm, lattice::lattice &system, bool label_on,
                  bool size_on, const double *size_mat,
                  bool spin_on, const double *spin_mat)
{
    assert(system.dim == 2);
    double r [2],
           rj[2];
    double dsize, resize_s, resize_c, min_c;
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
    resize_c = 1. / 1;
    resize_s = 1. / 10;
    min_c = find_min(size_mat, system.n);
    if (size_on && system.n > 8)
        resize_c = 0.5 / (size_mat[8] - min_c);
    if (spin_on && system.n > 8)
        resize_s = 0.1 / abs(spin_mat[8]);
    for (int i = 1; i < system.n; i++) {
        system.r(r, i);
        // Label
        if (label_on)
            mp_fid << "label.bot(btex" << i << "," << system.calc_rmin(0, i) << "etex,("
                   << r[0] * x_scal + x_offset << "cm,"
                   << r[1] * y_scal + y_offset << "cm));" << endl;
        // Lattice point or density correlation
        dsize = size_on ? (size_mat[i] - min_c) * resize_c + 0.1 : 4.0;
        mp_fid << "pickup pencircle scaled " << dsize << "cm;" << endl << "drawdot ("
               << r[0] * x_scal + x_offset << "cm,"
               << r[1] * y_scal + y_offset << "cm) withcolor (.0,1.,.0);" << endl;
        // Spin correlation
        if (spin_on) {
            dsize *= spin_mat[i] * resize_s;
            mp_fid << "pickup pencircle scaled 2pt;" << endl;
            mp_fid << "draw ("
                   << r[0] * x_scal + x_offset << "cm,"
                   << r[1] * y_scal + y_offset - dsize << "cm)--("
                   << r[0] * x_scal + x_offset << "cm,"
                   << r[1] * y_scal + y_offset + dsize << "cm) withcolor ("
                   << int(dsize > 0) << ".,.0," << int(dsize < 0) << ".);" << endl;
        }
    }

    // Metapost file ending.
    mp_fid << "endfig;" << endl 
           << "end"     << endl;
}

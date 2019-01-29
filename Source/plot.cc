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

void plot_lattice(const char *mp_fnm, lattice::lattice &system)
{
    assert(system.dim == 2);
    double r[2], rj[2];
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
        mp_fid << "pickup pencircle scaled 4pt;" << endl << "drawdot (" 
               << r[0] * x_scal + x_offset << "cm," 
               << r[1] * y_scal + y_offset << "cm) withcolor (.0,.0,1.);" << endl;
    }

    // Metapost file ending.
    mp_fid << "endfig;" << endl 
           << "end"     << endl;
}
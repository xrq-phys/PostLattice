#include "green_gen.hh"
#include <cassert>
#include <fstream>
#include <iomanip>
#include <string>

// Buffer Filename
#ifndef tm_fnm
#define tm_fnm ("/tmp/greentwo.tmp")
#endif

using namespace std;

void write_ca(const char *ca_fnm, lattice::lattice &system)
{
    fstream ca_fid(ca_fnm, fstream::out);

    ca_fid << "=================================" << endl
           << " CiAj " << system.n * system.n * 2 << endl
           << "=================================" << endl
           << " CiAj Full ======================" << endl
           << "=================================" << endl;
    for (int i = 0; i < system.n; i++)
        for (int j = 0; j < system.n; j++) {
            ca_fid << i << ' ' << 0 << ' ' << j << ' ' << 0 << endl
                   << i << ' ' << 1 << ' ' << j << ' ' << 1 << endl;
        }
    ca_fid.close();
}

void write_caca(const char *caca_fnm, lattice::lattice &system, bool sc)
{
    fstream tmop_fid(tm_fnm, fstream::out);
    int nopr = 0;
    
    for (int i = 0; i < system.n; i++)
        for (int j = 0; j < system.n; j++)
            for (int k = 0; k < system.n; k++)
                for (int l = 0; l < system.n; l++) {
                    if (i == j && k == l) {
                        tmop_fid << i << ' ' << 0 << ' ' << j << ' ' << 0 << ' ' << k << ' ' << 0 << ' ' << l << ' ' << 0 << endl
                                 << i << ' ' << 1 << ' ' << j << ' ' << 1 << ' ' << k << ' ' << 1 << ' ' << l << ' ' << 1 << endl
                                 << i << ' ' << 1 << ' ' << j << ' ' << 1 << ' ' << k << ' ' << 0 << ' ' << l << ' ' << 0 << endl
                                 << i << ' ' << 0 << ' ' << j << ' ' << 0 << ' ' << k << ' ' << 1 << ' ' << l << ' ' << 1 << endl;
                        nopr += 4;
                        if (i != k) {
                            tmop_fid << i << ' ' << 1 << ' ' << l << ' ' << 1 << ' ' << k << ' ' << 0 << ' ' << j << ' ' << 0 << endl
                                     << i << ' ' << 0 << ' ' << l << ' ' << 0 << ' ' << k << ' ' << 1 << ' ' << j << ' ' << 1 << endl;
                            nopr += 2;
                        }
                    } else if (sc && system.nn[i][k] && system.nn[j][l]) {
                        if (i != l || j != k) {
                            tmop_fid << i << ' ' << 1 << ' ' << j << ' ' << 1 << ' ' << k << ' ' << 0 << ' ' << l << ' ' << 0 << endl
                                     << i << ' ' << 0 << ' ' << j << ' ' << 0 << ' ' << k << ' ' << 1 << ' ' << l << ' ' << 1 << endl;
                            nopr += 2;
                        }
                    }
                }
    tmop_fid.close();
    
    fstream tmin_fid(tm_fnm);
    fstream caca_fid(caca_fnm, fstream::out);
    caca_fid << "=============================" << endl
             << "   CiAjCkAlSCAF      " << nopr << endl
             << "=============================" << endl
             << "== SC/AF Correlation ========" << endl
             << "=============================" << endl;
    string io_buffer;
    while (true) {
        getline(tmin_fid, io_buffer);
        if (tmin_fid.eof()) 
            break;
        caca_fid << io_buffer << endl;
    }
    tmin_fid.close();
    caca_fid.close();
}

#include "green_gen.hh"
#include <cassert>
#include <fstream>
#include <iomanip>
#include <string>

// Buffer Filename
#ifndef tm_fnm
#define tm_fnm ("/tmp/greentwo.tmp")
#endif

void write_ca(const char *ca_fnm, lattice::lattice &system)
{
    std::fstream ca_fid(ca_fnm, std::fstream::out);

    ca_fid << "=================================" << std::endl
           << " CiAj " << system.n * system.n * 2 << std::endl
           << "=================================" << std::endl
           << " CiAj Full ======================" << std::endl
           << "=================================" << std::endl;
    for (int i = 0; i < system.n; i++)
        for (int j = 0; j < system.n; j++) {
            ca_fid << i << ' ' << 0 << ' ' << j << ' ' << 0 << std::endl
                   << i << ' ' << 1 << ' ' << j << ' ' << 1 << std::endl;
        }
    ca_fid.close();
}

void write_caca(const char *caca_fnm, lattice::lattice &system)
{
    std::fstream tmop_fid(tm_fnm, std::fstream::out);
    int nopr = 0;
    
    for (int i = 0; i < system.n; i++)
        for (int j = 0; j < system.n; j++)
            for (int k = 0; k < system.n; k++)
                for (int l = 0; l < system.n; l++) {
                    if (i == j && k == l) {
                        tmop_fid << i << ' ' << 0 << ' ' << j << ' ' << 0 << ' ' << k << ' ' << 0 << ' ' << l << ' ' << 0 << std::endl 
                                 << i << ' ' << 1 << ' ' << j << ' ' << 1 << ' ' << k << ' ' << 1 << ' ' << l << ' ' << 1 << std::endl
                                 << i << ' ' << 1 << ' ' << j << ' ' << 1 << ' ' << k << ' ' << 0 << ' ' << l << ' ' << 0 << std::endl
                                 << i << ' ' << 0 << ' ' << j << ' ' << 0 << ' ' << k << ' ' << 1 << ' ' << l << ' ' << 1 << std::endl;
                        nopr += 4;
                        if (i != k) {
                            tmop_fid << i << ' ' << 1 << ' ' << l << ' ' << 1 << ' ' << k << ' ' << 0 << ' ' << j << ' ' << 0 << std::endl
                                     << i << ' ' << 0 << ' ' << l << ' ' << 0 << ' ' << k << ' ' << 1 << ' ' << j << ' ' << 1 << std::endl;
                            nopr += 2;
                        }
                    } else if (system.nn[i][k] && system.nn[j][l]) {
                        tmop_fid << i << ' ' << 1 << ' ' << j << ' ' << 1 << ' ' << k << ' ' << 0 << ' ' << l << ' ' << 0 << std::endl
                                 << i << ' ' << 0 << ' ' << j << ' ' << 0 << ' ' << k << ' ' << 1 << ' ' << l << ' ' << 1 << std::endl;
                        nopr += 2;
                    }
                }
    tmop_fid.close();
    
    std::fstream tmin_fid(tm_fnm);
    std::fstream caca_fid(caca_fnm, std::fstream::out);
    caca_fid << "=============================" << std::endl
             << "   CiAjCkAlSCAF      " << nopr << std::endl
             << "=============================" << std::endl
             << "== SC/AF Correlation ========" << std::endl
             << "=============================" << std::endl;
    std::string io_buffer;
    while (true) {
        std::getline(tmin_fid, io_buffer);
        if (tmin_fid.eof()) 
            break;
        caca_fid << io_buffer << std::endl;
    }
    tmin_fid.close();
    caca_fid.close();
}

/**
 * @file core.cc
 * Main method. Parses input from ini and calls specific routines.
 */
#include "INIReader.h"
#include "plot.hh"
#include "green_gen.hh"
#include "measure.hh"
#include "lattice.hh"
#include "trig_lattice.hh"
#include "honeycomb.hh"
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 2)
    { cerr << "Usage: Post <Input>.ini" << endl; return 0; }
    INIReader inp_fid(argv[1]);
    if (inp_fid.ParseError() < 0)
    { cerr << "Bad or non-existing input file " << argv[1] << endl; abort(); }

    lattice::lattice *physics;
    string lattice_n = inp_fid.Get("Physics", "System", "UNDEFINED");
    long W = inp_fid.GetInteger("Physics", "W", 1);
    long L = inp_fid.GetInteger("Physics", "L", W);
    long a0W = inp_fid.GetInteger("Physics", "a0W", W);
    long a1L = inp_fid.GetInteger("Physics", "a1L", L);
    long a1W = inp_fid.GetInteger("Physics", "a1W", 0);
    if (lattice_n == "Square")
        physics = new lattice::square2d(a0W, a1L);
    else if (lattice_n == "Triangular")
        physics = new lattice::trig2d(a0W, a1L, a1W);
    else if (lattice_n == "Honeycomb")
        physics = new lattice::honeycomb(a0W, a1L, a1W);
    else
    { cerr << "System " << lattice_n << " is not supported." << endl; abort(); }

    operator_options options;
    options.direct = inp_fid.GetBoolean("Control", "Direct", false);
    options.verbose = inp_fid.GetBoolean("Control", "Verbose", false);
    options.g1e_fnm = inp_fid.Get("Control", "GreenOne", "Green1.txt");
    options.g2e_fnm = inp_fid.Get("Control", "GreenTwo", "Green2.txt");
    string mode_n = inp_fid.Get("Control", "Mode", "UNDEFINED");
    // Allow command-line overriding of control parameters.
    if (argc > 3) {
        options.g1e_fnm = argv[2];
        options.g2e_fnm = argv[3];
        mode_n = "Measure";
        if (argc == 5 && argv[4][0] == '-' && argv[4][1] == 'v')
            options.verbose = true;
    }

    // Execute
    if (mode_n == "Gen") {
        write_ca  (options.g1e_fnm.c_str(), *physics);
        write_caca(options.g2e_fnm.c_str(), *physics,
                   inp_fid.Get("Operator", "SC", "-").c_str()[0] != '-');
    } else if (mode_n == "Measure") {
        options.chk_duplicate = inp_fid.GetBoolean("Control", "Check_Duplicate", true);
        options.db = inp_fid.GetBoolean("Operator", "Doublon", true);
        options.af = inp_fid.GetBoolean("Operator", "Spin_Structure", true);
        options.sc = inp_fid.Get("Operator", "SC", "-").c_str()[0];
        if (options.sc != '-') {
            options.sc_n = inp_fid.GetInteger("Operator", "SC_NumR", 6);
            options.sc_stat = inp_fid.Get("Operator", "SC_Stat", "S").c_str()[0];
        }
        if (options.af)
            switch (physics->dim) {
            case 2:
                options.af_n = new int[2];
                options.af_n[0] = inp_fid.GetInteger("Operator", "AF_N_X", 4);
                options.af_n[1] = inp_fid.GetInteger("Operator", "AF_N_Y", 4);
                break;
            case 3:
                options.af_n = new int[3];
                options.af_n[0] = inp_fid.GetInteger("Operator", "AF_N_X", 4);
                options.af_n[1] = inp_fid.GetInteger("Operator", "AF_N_Y", 4);
                options.af_n[2] = inp_fid.GetInteger("Operator", "AF_N_Z", 4);
                break;
            default:
                break;
            }
        options.db_fnm = inp_fid.Get("Operator", "Doublon_Out", "doublon.txt").c_str();
        options.sc_fnm = inp_fid.Get("Operator", "SC_Out", "sc.txt").c_str();
        options.af_fnm = inp_fid.Get("Operator", "AF_Out", "af.txt").c_str();
        measure(*physics, options);
    } else if (mode_n == "Plot") {
        options.mp_fnm = inp_fid.Get("Control", "Image", "lattice.mp");
        plot_lattice(options.mp_fnm.c_str(), *physics);
    } else
    { cerr << "Unsupported job name " << mode_n << endl; abort(); }

    return 0;
}

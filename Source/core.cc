/**
 * @file core.cc
 * Main method. Parses input from ini and calls specific routines.
 */
#include "INIReader.h"
#include "green_gen.hh"
#include "measure.hh"
#include "lattice.hh"
#include "trig_lattice.hh"
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
    if (lattice_n == "Square") {
        long W = inp_fid.GetInteger("Physics", "W", 1);
        physics = new lattice::square2d(W);
    } else if (lattice_n == "Triangular") {
        long a0W = inp_fid.GetInteger("Physics", "a0W", 1),
             a1L = inp_fid.GetInteger("Physics", "a1L", 1),
             a1W = inp_fid.GetInteger("Physics", "a1W", 0);
        physics = new lattice::trig2d(a0W, a1L, a1W);
    } else
    { cerr << "System " << lattice_n << " is not supported." << endl; abort(); }

    operator_options options;
    options.direct = inp_fid.GetBoolean("Control", "Direct", false);
    options.verbose = inp_fid.GetBoolean("Control", "Verbose", false);
    options.g1e_fnm = inp_fid.Get("Control", "GreenOne", "Green1.txt").c_str();
    options.g2e_fnm = inp_fid.Get("Control", "GreenTwo", "Green2.txt").c_str();
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
        write_ca  (options.g1e_fnm, *physics);
        write_caca(options.g2e_fnm, *physics,
                   inp_fid.Get("Operator", "SC", "-").c_str()[0] != '-');
    } else if (mode_n == "Measure") {
        options.db = inp_fid.GetBoolean("Operator", "Doublon", true);
        options.af = inp_fid.GetBoolean("Operator", "Spin_Structure", true);
        options.sc = inp_fid.Get("Operator", "SC", "-").c_str()[0];
        if (options.sc != '-')
            options.sc_n = inp_fid.GetInteger("Operator", "SC_NumR", 6);
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
    } else
    { cerr << "Unsupported job name " << mode_n << endl; abort(); }

    return 0;
}

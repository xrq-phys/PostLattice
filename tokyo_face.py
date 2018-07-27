#!/usr/bin/env python

import re
import sys
import numpy as np
from os import path

import operators
import lattice

#=======================================================================
# Loading Phase

runmode = -1
flip_index = [ -1 ]
a_length = [ -1, -1 ]
sys_type = "unavailable"
prefix = "unavailable"
suffix = "unavailable"
system_fnm = "system.py"
if (len(sys.argv) > 1):
    system_fnm = sys.argv[1]
if (not path.isfile(system_fnm)):
    print("System file does not exist.")
    sys.exit()
exec(open(system_fnm).read())

if (runmode == 'wick'):
    runmode = 1
elif (runmode == 'direct'):
    runmode = 0
else:
    print("Bad Run Mode")
    sys.exit()

if (sys_type != "square2d"):
    print("This system is not supported.")
    sys.exit()
if (a_length[0] < 0 or a_length[1] < 0):
    print("Invalid A.")
    sys.exit()
phys = lattice.square_2d(a_length)

if (len(sys.argv) > 3):
    greenone_fnm = sys.argv[2]
    greentwo_fnm = sys.argv[3]
else:
    greenone_fnm = "output/" + prefix + "_cisajs" + suffix + ".dat"
    greentwo_fnm = "output/" + prefix + "_cisajscktalt" + suffix + ".dat"
if (not path.isfile(greenone_fnm) or not path.isfile(greentwo_fnm)):
    print("Green function file does not exist: %s and %s" % (greenone_fnm, greentwo_fnm))
    sys.exit()

# Load 1-body part
greenone_dat = np.zeros([ 2 * phys.n, 2 * phys.n ], float)
greenone_fid = open(greenone_fnm, 'r')
greenone_fln = greenone_fid.readline()
while(greenone_fln.strip() != ''):
    greenone_arr = re.sub(' +', ' ', greenone_fln.strip()).split(' ')
    greenone_dat[int(greenone_arr[0]) * 2 + int(greenone_arr[1])]\
                [int(greenone_arr[2]) * 2 + int(greenone_arr[3])] = float(greenone_arr[4])
    greenone_fln = greenone_fid.readline()
greenone_fid.close()

#=======================================================================
# Measurement Phase

doublon = operators.double_occ(phys)
sstruct = operators.spin_struct(phys, [[ 2. * np.pi * (i % phys.a[1]) / phys.a[0], \
                                         2. * np.pi * (i / phys.a[1]) / phys.a[1] ]\
                                         for i in range(0, phys.n)])
if (len(sys.argv) < 5 or sys.argv[4][0] != '-'):
    sc = operators.sc_corr(phys, [0,   # [ 0, 0 ]
                                  1,   # [ 1, 0 ]
                                  2,   # [ 1, 1 ]
                                  4,   # [ 2, 0 ]
                                  5,   # [ 2, 1 ]
                                  8,   # [ 2, 2 ]
                                  9,   # [ 3, 0 ]
                                  10]) # [ 3, 1 ] 

seen_table = {}
greentwo_lnm = 0
greentwo_fid = open(greentwo_fnm, 'r')
greentwo_fln = greentwo_fid.readline()
while(greentwo_fln.strip() != ''):
    greentwo_arr = re.sub(' +', ' ', greentwo_fln.strip()).split(' ')
    i  =   int(greentwo_arr[0])
    si =   int(greentwo_arr[1])
    j  =   int(greentwo_arr[2])
    sj =   int(greentwo_arr[3])
    k  =   int(greentwo_arr[4])
    sk =   int(greentwo_arr[5])
    l  =   int(greentwo_arr[6])
    sl =   int(greentwo_arr[7])
    x  = float(greentwo_arr[8])
    if (not (i, si, j, sj, k, sk, l, sl) in seen_table):
        seen_table[(i, si, j, sj, k, sk, l, sl)] = True
        doublon.measure(i, si, j, sj, k, sk, l, sl, x)
        if (not runmode or greentwo_lnm > len(flip_index) or not flip_index[greentwo_lnm]):
            sstruct.measure(i, si, j, sj, k, sk, l, sl, x)
        else:
            if (k == l):
                sstruct.measure(i, si, l, sl, k, sk, j, sj, greenone_dat[i + si][j + sj] - x)
            else:
                sstruct.measure(i, si, l, sl, k, sk, j, sj, -x)
        if (len(sys.argv) < 5 or sys.argv[4][0] != '-'):
            sc.measure(i, si, k, sk, j, sj, l, sl, -x)
    greentwo_fln = greentwo_fid.readline()
    greentwo_lnm += 1
greentwo_fid.close()

#=======================================================================
# Add contribution from one-body Green function to SC or other operators 
# that consists of 2-body operators like: c+ c+ c c.

if (len(sys.argv) > 4 and sys.argv[4][0] == '+'):

    print("Adding contributions from 1-body Green function...")
    print("WARNING: This procedure is extremely time-consuming and unparallelized.")
    print("WARNING: So it's strongly recommanded that you use an OMP-C++ implementation instead.")
    for i in range(2 * phys.n):
        for j in range(2 * phys.n): 
            ri = int(i / 2)
            rj = int(j / 2)
            si = i % 2
            sj = j % 2
            # To avoid measuring zero:
            if (np.abs(greenone_dat[i][j]) > 1e-4):
                # Add the contribution:
                # < c+_i,s c+_k,s' c_l,s' c_j,s >.
                for k in range(phys.n):
                    for s in range(2):
                        sc.measure(ri, si, k, s, k, s, rj, sj, greenone_dat[i][j])

#=======================================================================
# Output Phase

doublon_fid = open("doublon.txt", 'w')
doublon_fid.write("%.10e\n" % doublon.value)
doublon_fid.close()

if (len(sys.argv) < 5 or sys.argv[4][0] != '-'):
    sc_fid = open("sc.txt", 'w')
    for ir in range(np.size(sc.rc_list)):
        sc_fid.write("%f %.10e\n" % (np.sqrt(sc.rc_list[ir]), sc.values[ir]))
    sc_fid.close()

sstruct_fid = open("sstruct.txt", 'w')
for y in range(phys.a[0]):
    for x in range(phys.a[1]):
        sstruct_fid.write("%d %d %lf %lf\n" % (x, y, np.real(sstruct.values[x + y * phys.a[1]]), \
                                                     np.imag(sstruct.values[x + y * phys.a[1]])))
sstruct_fid.close()

spn_nn_fid = open("spn_nn.txt", 'w')
spn_nn_fid.write("%.10e\n" % (sstruct.val_nn / phys.n / 4.))
spn_nn_fid.close()

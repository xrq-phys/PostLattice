#/usr/bin/env python

import re
import sys
import numpy as np
from os import path

import operators
import lattice

system_fnm = "system.py"
if (len(sys.argv) > 1):
    system_fnm = sys.argv[1]
if (not path.isfile(system_fnm)):
    print("System file does not exist.")
    sys.exit()
exec(open(system_fnm).read())

if (sys_type != "square2d"):
    print("This system is not supported.")
    sys.exit()
phys = lattice.square_2d(a_length)

if (not path.isfile(rdm_fnm)):
    print("Green Function file does not exist.")
    sys.exit()

doublon = operators.double_occ(phys)
sstruct = operators.spin_struct(phys, [[ 2. * np.pi * (i % phys.a[1]) / phys.a[0], \
                                         2. * np.pi * (i / phys.a[1]) / phys.a[1] ]\
                                         for i in range(0, phys.n)])
if (min(phys.a) >= 6):
    sc = operators.sc_corr(phys, 3)

rdm_fid = open(rdm_fnm, 'r')
rdm_fln = rdm_fid.readline()
while(rdm_fln.strip() != ''):
    rdm_arr = re.sub(' +', ' ', rdm_fln.strip()).split(' ')
    p  = int((int(rdm_arr[0]) - 1) / 2)
    sp = int((int(rdm_arr[0]) - 1) % 2)
    r  = int((int(rdm_arr[1]) - 1) / 2)
    sr = int((int(rdm_arr[1]) - 1) % 2)
    s  = int((int(rdm_arr[2]) - 1) / 2)
    ss = int((int(rdm_arr[2]) - 1) % 2)
    q  = int((int(rdm_arr[3]) - 1) / 2)
    sq = int((int(rdm_arr[3]) - 1) % 2)
    x  = float(rdm_arr[4])
    doublon.measure(p, sp, s, ss, r, sr, q, sq, -x)
    sstruct.measure(p, sp, s, ss, r, sr, q, sq, -x)
    if (min(phys.a) >= 6):
        sc.measure(p, sp, r, sr, s, ss, q, sq, x)
    rdm_fln = rdm_fid.readline()
rdm_fid.close()

doublon_fid = open("doublon.txt", 'w')
doublon_fid.write("%lf" % doublon.value)
doublon_fid.close()
if (min(phys.a) >= 6):
    sc_fid = open("sc.txt", 'w')
    sc_fid.write("%lf" % sc.value)
    sc_fid.close()
sstruct_fid = open("sstruct.txt", 'w')
for y in range(phys.a[0]):
    for x in range(phys.a[1]):
        sstruct_fid.write("%d %d %lf %lf\n" % (x, y, np.real(sstruct.values[x + y * phys.a[1]]), \
                                                     np.imag(sstruct.values[x + y * phys.a[1]])))
sstruct_fid.close()


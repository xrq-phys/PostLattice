#!/usr/bin/env python

import re
import numpy as np
from glob import glob
from operator import itemgetter
from matplotlib import pyplot as plt

sc_rnum   = 7
total_bin = 8.
calc_list = glob('U*')

delta_list   = np.zeros([len(calc_list)], float)
energy_list  = np.zeros([len(calc_list)], float)
doublon_list = np.zeros([len(calc_list)], float)
sstruct_list = np.zeros([len(calc_list)], complex)
sc_list      = np.zeros([len(calc_list)], float)
energy_err   = np.zeros([len(calc_list)], float)
doublon_err  = np.zeros([len(calc_list)], float)
sstruct_err  = np.zeros([len(calc_list)], float)
sc_err       = np.zeros([len(calc_list)], float)

for i in range(len(calc_list)):

    #===================
    # System Information
    calc_dir = calc_list[i]
    exec(open(calc_dir + '/model.def'))
    delta_list[i] = float(W * L - nelec) / (W * L)

    #====================
    # Energy and Doublon 
    energy_dat  = np.zeros([2], float)
    doublon_dat = np.zeros([2], float)

    for energy_fnm in glob(calc_dir + '/output/' + 'zvo_out_*.dat'):
        energy_fid = open(energy_fnm)
        energy_fln = re.sub(' +', ' ', energy_fid.readline().strip()).split(' ')
        energy_dat[0] += float(energy_fln[0])
        energy_dat[1] += float(energy_fln[0]) ** 2
        energy_fid.close()
    energy_list[i] = energy_dat[0]/total_bin
    energy_err [i] = np.sqrt(energy_dat[1]/total_bin - energy_list[i] ** 2)

    for doublon_fnm in glob(calc_dir + '/doublon.txt.*.dat.regular'):
        doublon_fid = open(doublon_fnm)
        doublon_fln = doublon_fid.readline()
        doublon_dat[0] += float(doublon_fln)
        doublon_dat[1] += float(doublon_fln) ** 2
        doublon_fid.close()
    doublon_list[i] = doublon_dat[0]/total_bin
    doublon_err [i] = np.sqrt(doublon_dat[1]/total_bin - doublon_list[i] ** 2)

    #======================
    # Spin and SC Structure
    sstruct_psz = 0
    sstruct_tmp = open(glob(calc_dir + '/sstruct.txt.*.dat.regular')[0])
    while(sstruct_tmp.readline().strip() != ''):
        sstruct_psz += 1
    sstruct_lkp = np.zeros([sstruct_psz, 2], int)
    sstruct_dat = np.zeros([sstruct_psz, 2], complex)
    sc_rc       = np.zeros([7             ], float)
    sc_dat      = np.zeros([7,           2], float)

    for sstruct_fnm in glob(calc_dir + '/sstruct.txt.*.dat.regular'):
        sstruct_fid = open(sstruct_fnm)
        for j in range(sstruct_psz):
            sstruct_fln = re.sub(' +', ' ', sstruct_fid.readline().strip()).split(' ')
            sstruct_lkp[j][0] = int(sstruct_fln[0])
            sstruct_lkp[j][1] = int(sstruct_fln[1])
            sstruct_dat[j][0] += np.sqrt(float(sstruct_fln[2]) ** 2 + float(sstruct_fln[3]) ** 2)
            sstruct_dat[j][1] += float(sstruct_fln[2]) ** 2 + float(sstruct_fln[3]) ** 2
        sstruct_fid.close()
    sstruct_mid = max(enumerate(np.real(sstruct_dat[:, 0])), key=itemgetter(1))[0]
    sstruct_list[i] = sstruct_dat[sstruct_mid][0]/total_bin
    sstruct_err [i] = np.sqrt(sstruct_dat[sstruct_mid][1]/total_bin - np.abs(sstruct_list[i]) ** 2)
    print("i = %d, Peak of spin structure factor @ the %dth q-point." % (i, sstruct_mid))

    for sc_fnm in glob(calc_dir + '/sc.txt.*.dat.sc'):
        sc_fid = open(sc_fnm)
        for j in range(sc_rnum):
            sc_fln = re.sub(' +', ' ', sc_fid.readline().strip()).split(' ')
            if (sc_rc[j] < 1e-6):
                sc_rc[j] = float(sc_fln[0])
            sc_dat[j][0] += float(sc_fln[1])
            sc_dat[j][1] += float(sc_fln[1]) ** 2
        sc_fid.close()
    sc_list[i] = sc_dat[j][0]/total_bin
    sc_err [i] = np.sqrt(sc_dat[j][1]/total_bin - sc_list[i] ** 2)

#=====
# Sort
energy_list  = [ x for _, x in sorted(zip(delta_list, energy_list )) ]
energy_err   = [ x for _, x in sorted(zip(delta_list, energy_err  )) ]
doublon_list = [ x for _, x in sorted(zip(delta_list, doublon_list)) ]
doublon_err  = [ x for _, x in sorted(zip(delta_list, doublon_err )) ]
sstruct_list = [ x for _, x in sorted(zip(delta_list, sstruct_list)) ]
sstruct_err  = [ x for _, x in sorted(zip(delta_list, sstruct_err )) ]
sc_list      = [ x for _, x in sorted(zip(delta_list, sc_list     )) ]
sc_err       = [ x for _, x in sorted(zip(delta_list, sc_err      )) ]
delta_list   = sorted(delta_list)

#========================
# Plotting Energy/Doublon
fig_en, ax_e = plt.subplots()

color_e = 'tab:blue'
ax_e.set_xlabel('Doping Rate')
ax_e.set_ylabel('Energy', color=color_e)
ax_e.errorbar(delta_list, energy_list, energy_err, color=color_e)
ax_e.tick_params(axis='y')

ax_n = ax_e.twinx()

color_n = 'tab:green'
ax_n.set_ylabel('Doublon', color=color_n)
ax_n.errorbar(delta_list, doublon_list, doublon_err, color=color_n)
ax_n.tick_params(axis='y')

fig_en.tight_layout()

#===============
# Plotting AF/SC
fig_as, ax_a = plt.subplots()

color_a = 'tab:blue'
ax_a.set_xlabel('Doping Rate')
ax_a.set_ylabel('AF Parameter', color=color_a)
ax_a.errorbar(delta_list, sstruct_list, sstruct_err, color=color_a)
ax_a.tick_params(axis='y')

ax_s = ax_a.twinx()

color_s = 'tab:red'
ax_s.set_ylabel('SC Parameter', color=color_s)
ax_s.errorbar(delta_list, sc_list, sc_err, color=color_s)
ax_s.tick_params(axis='y')

plt.show()

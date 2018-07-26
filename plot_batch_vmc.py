#!/usr/bin/env python

import sys
import re
import numpy as np
from glob import glob
from operator import itemgetter
import matplotlib.pyplot as plt
import matplotlib

sc_rc     = 2.
total_bin = 0.
calc_list = glob('U*')
gf_suffix = '*'

if (len(sys.argv) > 1):
    calc_list = glob(sys.argv[1])
if (len(sys.argv) > 2):
    gf_suffix = sys.argv[2]

delta_list   = np.zeros([len(calc_list)], float)
energy_list  = np.zeros([len(calc_list)], float)
doublon_list = np.zeros([len(calc_list)], float)
sstruct_list = np.zeros([len(calc_list)], float)
sc_list      = np.zeros([len(calc_list)], float)
energy_err   = np.zeros([len(calc_list)], float)
doublon_err  = np.zeros([len(calc_list)], float)
sstruct_err  = np.zeros([len(calc_list)], float)
sc_err       = np.zeros([len(calc_list)], float)

for i in range(len(calc_list)):

    #===================
    # System Information
    calc_dir = calc_list[i]
    W = 0; L = 0; nelec = -1
    exec(open(calc_dir + '/model.def'))
    delta_list[i] = float(W * L - nelec) / (W * L)
    if (total_bin < 1E-1):
        total_bin = len(glob(calc_dir + '/output/' + 'zvo_out_' + gf_suffix + '.dat'))

    #====================
    # Energy and Doublon 
    energy_dat  = np.zeros([2], float)
    doublon_dat = np.zeros([2], float)

    for energy_fnm in glob(calc_dir + '/output/' + 'zvo_out_' + gf_suffix + '.dat'):
        energy_fid = open(energy_fnm)
        energy_fln = re.sub(' +', ' ', energy_fid.readline().strip()).split(' ')
        energy_dat[0] += float(energy_fln[0])
        energy_dat[1] += float(energy_fln[0]) ** 2
        energy_fid.close()
    energy_list[i] = energy_dat[0]/total_bin
    energy_err [i] = np.sqrt(energy_dat[1]/total_bin - energy_list[i] ** 2)/np.sqrt(total_bin - 1)

    for doublon_fnm in glob(calc_dir + '/doublon.txt.' + gf_suffix + '.dat'):
        doublon_fid = open(doublon_fnm)
        doublon_fln = doublon_fid.readline()
        doublon_dat[0] += float(doublon_fln)
        doublon_dat[1] += float(doublon_fln) ** 2
        doublon_fid.close()
    doublon_list[i] = doublon_dat[0]/total_bin
    doublon_err [i] = np.sqrt(doublon_dat[1]/total_bin - doublon_list[i] ** 2)/np.sqrt(total_bin - 1)

    #======================
    # Spin and SC Structure
    sstruct_psz = 0
    sstruct_tmp = open(glob(calc_dir + '/sstruct.txt.' + gf_suffix + '.dat')[0])
    while(sstruct_tmp.readline().strip() != ''):
        sstruct_psz += 1
    sstruct_tmp.close()
    sstruct_lkp = np.zeros([sstruct_psz, 2], int)
    sstruct_dat = np.zeros([sstruct_psz, 2], complex)
    sc_dat = np.zeros([2], float)

    for sstruct_fnm in glob(calc_dir + '/sstruct.txt.' + gf_suffix + '.dat'):
        sstruct_fid = open(sstruct_fnm)
        for j in range(sstruct_psz):
            sstruct_fln = re.sub(' +', ' ', sstruct_fid.readline().strip()).split(' ')
            sstruct_lkp[j][0] = int(sstruct_fln[0])
            sstruct_lkp[j][1] = int(sstruct_fln[1])
            sstruct_dat[j][0] += float(sstruct_fln[2]) * 1. + float(sstruct_fln[3]) * 1j
            sstruct_dat[j][1] += float(sstruct_fln[2]) ** 2 + float(sstruct_fln[3]) ** 2
        sstruct_fid.close()
    sstruct_mid = max(enumerate(np.abs(sstruct_dat[:, 0])), key=itemgetter(1))[0]
    sstruct_list[i] = np.abs(sstruct_dat[sstruct_mid][0])/total_bin
    sstruct_err [i] = np.sqrt(sstruct_dat[sstruct_mid][1]/total_bin - np.abs(sstruct_list[i]) ** 2)/\
                      np.sqrt(total_bin - 1)
    print("delta = %f, Peak of spin structure factor @ q-point (%fPi, %fPi)." % \
          (delta_list[i], sstruct_lkp[sstruct_mid][0] * 2./W, sstruct_lkp[sstruct_mid][1] * 2./W))

    for sc_fnm in glob(calc_dir + '/sc.txt.' + gf_suffix + '.dat'):
        sc_fid = open(sc_fnm)
        sc_cur = 0.
        sc_cnt = 0
        sc_fln = sc_fid.readline().strip()
        while (sc_fln != ''):
            sc_fln = re.sub(' +', ' ', sc_fln).split(' ')
            if (float(sc_fln[0]) - sc_rc > -1E-2):
                sc_cnt += 1
                sc_cur += float(sc_fln[1])
            sc_fln = sc_fid.readline().strip()
        sc_fid.close()
        sc_cur /= sc_cnt
        sc_dat[0] += sc_cur
        sc_dat[1] += sc_cur ** 2
    sc_list[i] = sc_dat[0]/total_bin
    sc_err [i] = np.sqrt(sc_dat[1]/total_bin - sc_list[i] ** 2)/np.sqrt(total_bin - 1)

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

#==========================
# Write Value/Error To File
energy_ofid  = open("Energy_vs_Delta.dat", "w")
doublon_ofid = open("Doublon_vs_Delta.dat", "w")
sstruct_ofid = open("SpinStructure_vs_Delta.dat", "w")
sc_ofid      = open("SC_Correlation_vs_Delta.dat", "w")
for i in range(len(calc_list)):
    energy_ofid.write("%lf %lf %lf\n" % (delta_list[i], energy_list[i], energy_err[i]))
    doublon_ofid.write("%lf %lf %lf\n" % (delta_list[i], doublon_list[i], doublon_err[i]))
    sstruct_ofid.write("%lf %lf %lf\n" % (delta_list[i], sstruct_list[i], sstruct_err[i]))
    sc_ofid.write("%lf %lf %lf\n" % (delta_list[i], sc_list[i], sc_err[i]))
energy_ofid.close()
doublon_ofid.close()
sstruct_ofid.close()
sc_ofid.close()

#========================
# Plotting Energy/Doublon
if (matplotlib.__version__ >= '2.0.0'):
    plt.style.use('classic')
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

#!/usr/bin/env python

import sys
import re
import numpy as np
from glob import glob
from operator import itemgetter
import matplotlib.pyplot as plt
import matplotlib

sc_rc     = 2.
calc_list = glob('U*')

if (len(sys.argv) > 1):
    calc_list = glob(sys.argv[1])

delta_list   = np.zeros([len(calc_list)], float)
energy_list  = np.zeros([len(calc_list)], float)
doublon_list = np.zeros([len(calc_list)], float)
sstruct_list = np.zeros([len(calc_list)], float)
sc_list      = np.zeros([len(calc_list)], float)

for i in range(len(calc_list)):

    #===================
    # System Information
    calc_dir = calc_list[i]
    exec(open(calc_dir + '/model.def'))
    delta_list[i] = float(W * L - nelec) / (W * L)

    #====================
    # Energy and Doublon 
    energy_fid = open(calc_dir + '/output/' + 'zvo_energy.dat')
    energy_fln = re.sub(' +', ' ', energy_fid.readline().strip()).split(' ')
    energy_fid.close()
    energy_list[i] = float(energy_fln[1])

    doublon_fid = open(calc_dir + '/doublon.txt')
    doublon_fln = doublon_fid.readline()
    doublon_fid.close()
    doublon_list[i] = float(doublon_fln)

    #======================
    # Spin and SC Structure
    sstruct_psz = 0
    sstruct_fid = open(calc_dir + '/sstruct.txt')
    while(sstruct_fid.readline().strip() != ''):
        sstruct_psz += 1
    sstruct_lkp = np.zeros([sstruct_psz, 2], int)
    sstruct_dat = np.zeros([sstruct_psz   ], float)

    sstruct_fid.seek(0)
    for j in range(sstruct_psz):
        sstruct_fln = re.sub(' +', ' ', sstruct_fid.readline().strip()).split(' ')
        sstruct_lkp[j][0] = int(sstruct_fln[0])
        sstruct_lkp[j][1] = int(sstruct_fln[1])
        sstruct_dat[j]   += np.sqrt(float(sstruct_fln[2]) ** 2 + float(sstruct_fln[3]) ** 2)
    sstruct_fid.close()
    sstruct_mid = max(enumerate(np.abs(sstruct_dat)), key=itemgetter(1))[0]
    sstruct_list[i] = sstruct_dat[sstruct_mid]
    print("delta = %f, Peak of spin structure factor @ q-point (%fPi, %fPi)." % \
          (delta_list[i], sstruct_lkp[sstruct_mid][0] * 2./W, sstruct_lkp[sstruct_mid][1] * 2./W))

    sc_fid = open(calc_dir + '/sc.txt')
    sc_cnt = 0
    sc_cur = 0.
    sc_fln = sc_fid.readline().strip()
    while (sc_fln != ''):
        sc_fln = re.sub(' +', ' ', sc_fln).split(' ')
        if (float(sc_fln[0]) - sc_rc > -1E-2):
            sc_cnt += 1
            sc_cur += float(sc_fln[1])
        sc_fln = sc_fid.readline().strip()
    sc_fid.close()
    sc_cur /= sc_cnt
    sc_list[i] = sc_cur

#=====
# Sort
energy_list  = [ x for _, x in sorted(zip(delta_list, energy_list )) ]
doublon_list = [ x for _, x in sorted(zip(delta_list, doublon_list)) ]
sstruct_list = [ x for _, x in sorted(zip(delta_list, sstruct_list)) ]
sc_list      = [ x for _, x in sorted(zip(delta_list, sc_list     )) ]
delta_list   = sorted(delta_list)

#==========================
# Write Value/Error To File
energy_ofid  = open("Energy_vs_Delta.dat", "w")
doublon_ofid = open("Doublon_vs_Delta.dat", "w")
sstruct_ofid = open("SpinStructure_vs_Delta.dat", "w")
sc_ofid      = open("SC_Correlation_vs_Delta.dat", "w")
for i in range(len(calc_list)):
    energy_ofid.write("%lf %lf\n" % (delta_list[i], energy_list[i]))
    doublon_ofid.write("%lf %lf\n" % (delta_list[i], doublon_list[i]))
    sstruct_ofid.write("%lf %lf\n" % (delta_list[i], sstruct_list[i]))
    sc_ofid.write("%lf %lf\n" % (delta_list[i], sc_list[i]))
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
ax_e.plot(delta_list, energy_list, color=color_e)
ax_e.tick_params(axis='y')

ax_n = ax_e.twinx()

color_n = 'tab:green'
ax_n.set_ylabel('Doublon', color=color_n)
ax_n.plot(delta_list, doublon_list, color=color_n)
ax_n.tick_params(axis='y')

fig_en.tight_layout()

#===============
# Plotting AF/SC
fig_as, ax_a = plt.subplots()

color_a = 'tab:blue'
ax_a.set_xlabel('Doping Rate')
ax_a.set_ylabel('AF Parameter', color=color_a)
ax_a.plot(delta_list, sstruct_list, color=color_a)
ax_a.tick_params(axis='y')

ax_s = ax_a.twinx()

color_s = 'tab:red'
ax_s.set_ylabel('SC Parameter', color=color_s)
ax_s.plot(delta_list, sc_list, color=color_s)
ax_s.tick_params(axis='y')

plt.show()

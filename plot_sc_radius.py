#!/usr/bin/env python

import re
import numpy as np
from glob import glob
from matplotlib import pyplot as plt

fnm_list = glob('sc.txt.*.dat.sc')
total_bin = len(fnm_list)
r_set = 0
r_list = []
sc_list = []

for fnm in fnm_list:
    i = 0
    fid = open(fnm)
    fln = fid.readline().strip()
    sc_tmp = np.zeros(r_set)
    while (fln != ''):
        fln = re.sub(' +', ' ', fln).split(' ')
        if (r_set == 0):
            r_list.append(float(fln[0]))
            sc_list.append([ float(fln[1]), 0 ])
        else:
            sc_tmp[i] = float(fln[1])
            i += 1
        fln = fid.readline().strip()
    if (r_set == 0):
        r_set = len(r_list)
        sc_list = np.array(sc_list)
        sc_tmp = np.array(sc_list[:, 0])
        sc_list[:, :] = 0
    for j in range(r_set - 1):
        sc_list[j, 0] +=  sc_tmp[j] - sc_tmp[j + 1]
        sc_list[j, 1] += (sc_tmp[j] - sc_tmp[j + 1]) ** 2
sc_list = sc_list/total_bin
sc_list[:, 1] = np.sqrt(sc_list[:, 1] - sc_list[:, 0] ** 2) #/np.sqrt(total_bin)

plt.errorbar(r_list, sc_list[:, 0], sc_list[:, 1])
plt.show()
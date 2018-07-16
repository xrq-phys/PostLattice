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
    while (fln != ''):
        fln = re.sub(' +', ' ', fln).split(' ')
        if (r_set == 0):
            r_list.append(float(fln[0]))
            sc_list.append([ float(fln[1]), \
                             float(fln[1]) ** 2 ])
        else:
            sc_list[i] += [ float(fln[1]), \
                            float(fln[1]) ** 2 ]
            i += 1
        fln = fid.readline().strip()
    if (r_set == 0):
        r_set = 1
        sc_list = np.array(sc_list)
sc_list = sc_list/total_bin
sc_list[:, 1] = np.sqrt(sc_list[:, 1] - sc_list[:, 0] ** 2) #/np.sqrt(total_bin)

plt.errorbar(r_list, sc_list[:, 0], sc_list[:, 1])
plt.show()

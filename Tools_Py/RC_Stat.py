#!/usr/bin/env python

import re
import sys
import numpy as np
from glob import glob

rc_lst    = []
rn_lst    = []
sc_lst    = []
sc_err    = []
calc_list = glob('*.dat')

if (len(sys.argv) > 1):
    calc_list = glob(sys.argv[1])
print(calc_list)

for calc_fnm in calc_list:
    calc_fid = open(calc_fnm)
    calc_fln = re.sub(' +', ' ', calc_fid.readline().strip()).split(' ')
    i = 0
    while (len(calc_fln) > 1):
        if (len(rc_lst) == i):
            rc_lst.append(float(calc_fln[0]))
            rn_lst.append(float(calc_fln[2]))
            sc_lst.append(float(calc_fln[1]))
            sc_err.append(float(calc_fln[1]) ** 2)
        else:
            sc_lst[i] += float(calc_fln[1])
            sc_err[i] += float(calc_fln[1]) ** 2
        calc_fln = re.sub(' +', ' ', calc_fid.readline().strip()).split(' ')
        i += 1

rc_lst = np.array(rc_lst)
sc_lst = np.array(sc_lst)
sc_err = np.array(sc_err)

sc_lst /= len(calc_list)
sc_err /= len(calc_list)
sc_err -= sc_lst ** 2
sc_err = np.sqrt(sc_err / (len(calc_list) - 1))

sc_out = open("rc_avg.dat", "w")
for i in xrange(len(rc_lst)):
    sc_out.write("%lf %lf %lf\n" % (rc_lst[i], sc_lst[i] / rn_lst[i],\
                                               sc_err[i] / rn_lst[i]))
sc_out.close()


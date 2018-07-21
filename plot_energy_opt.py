#!/usr/bin/env python

import re
import sys
from os import path
import matplotlib.pyplot as plt
import matplotlib

prefix = 'zvo'
suffix = '001'

if (len(sys.argv) > 2):
    prefix = sys.argv[1]
    prefix = sys.argv[2]

out_fnm = 'output/' + prefix + '_out_' + suffix + '.dat'
if (not path.isfile(out_fnm)):
    print("Out file not found: " + out_fnm)
    sys.exit()

eng_lst = []
out_fid = open(out_fnm)
out_fln = out_fid.readline()
while(out_fln.strip() != ''):
    out_arr = re.sub(' +', ' ', out_fln.strip()).split(' ')
    eng_lst.append(float(out_arr[0]))
    out_fln = out_fid.readline()

if (matplotlib.__version__ >= '2.0.0'):
    plt.style.use('classic')
plt.plot(range(len(eng_lst)), eng_lst)
plt.show()

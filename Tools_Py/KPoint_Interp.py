#!/usr/bin/env python

import sys
from os import path
from scipy import interpolate
import numpy as np

if (not path.isfile(sys.argv[1]) or not path.isfile(sys.argv[2])):
    print('Usage: KPoint_Interp.py SPIN_STRUCT_DATA_FILE QPOINTS_FILE.')
    print('Note: Only for 2D systems currently.')
    sys.exit()

dat_fid = open(sys.argv[1], 'r')
opt_fid = open(sys.argv[2], 'r')
out_fid = open('kline.dat', 'w')
nod_fid = open('knodes.gp', 'w')

x = []
y = []
z = []

# Read data file.
buff = dat_fid.readline().strip()
while (buff != ''):
    [ xi, yi, zi ] = np.array(buff.split(' '), float)
    x.append(xi)
    y.append(yi)
    z.append(zi)
    buff = dat_fid.readline().strip()

# Do interpolation.
f = interpolate.interp2d(x, y, z, kind='linear')

# Read options file and write output.
[ n_node, n_sample ] = np.array(opt_fid.readline().strip().split(' '), int)
buff = opt_fid.readline().strip().split()
[ xi, yi ] = np.array(buff[1:3], float)
nod_fid.write('set xtics (\'' + buff[0] + '\' 0.0')
stot = 0
for i in xrange(n_node - 1):
    xk = xi
    yk = yi
    buff = opt_fid.readline().strip().split()
    [ xi, yi ] = np.array(buff[1:3], float)
    si = np.sqrt((xk - xi) ** 2 + (yk - yi) ** 2) / n_sample
    for j in xrange(n_sample):
        xj = xk + (xi - xk) * j / n_sample
        yj = yk + (yi - yk) * j / n_sample
        out_fid.write("%f %lf\n" % (stot, f(xj, yj)))
        stot += si
    nod_fid.write(', \'' + buff[0] + '\' %lf' % stot)
nod_fid.write(')\n')
nod_fid.write('set grid xtics lt 1 lc 0')
out_fid.write("%f %lf\n" % (stot, f(xi, yi)))

dat_fid.close()
opt_fid.close()
nod_fid.close()
out_fid.close()


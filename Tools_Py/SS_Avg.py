#!/usr/bin/env python
import re
import sys
import numpy as np
from glob import glob

kX = 2
kY = 2
calc_list = glob('*.dat')
if (len(sys.argv) > 1):
    calc_list = glob(sys.argv[1])
if (len(sys.argv) == 4):
    kX = int(sys.argv[2])
    kY = int(sys.argv[3])

tot = 0.0
err = 0.0
val = 0.0
tmp = 0.0
num = len(calc_list)

for fnm in calc_list:
	fid = open(fnm)
	rdy = 0
	while(rdy < 1):
		fln = fid.readline().split()
		if(len(fln) < 3):
			continue
		if(int(fln[0]) == kX and int(fln[1]) == kY):
			tmp = float(fln[2])
			if (tmp > val):
				val = tmp
			rdy += 1
		#if(int(fln[0]) == W / 2 and int(fln[1]) == 0):
		#	tmp = float(fln[2])
		#	if (tmp > val):
		#		val = tmp
		#	rdy += 1
		#if(int(fln[0]) == 0 and int(fln[1]) == L / 2):
		#	tmp = float(fln[2])
		#	if (tmp > val):
		#		val = tmp
		#	rdy += 1
	tot += val
	err += val**2
	val = 0.0
	fid.close()
err -= tot**2 / num

print("%lf %lf" % (tot / num, np.sqrt(err) / num))

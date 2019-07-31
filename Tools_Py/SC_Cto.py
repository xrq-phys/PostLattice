#!/usr/bin/env python
import sys
import numpy as np
num = 0.0
tot = 0.0
err = 0.0
cto = float(sys.argv[2])
fid = open(sys.argv[1])
fln = [ '0', '0', '0' ]
while (len(fln) == 3):
	if (float(fln[0]) > cto):
		num += 1.0
		tot += float(fln[1])
		err += float(fln[2])**2
	fln = fid.readline().strip().split(' ')
fid.close()
print("%lf %lf" % (tot / num, np.sqrt(err) / num))


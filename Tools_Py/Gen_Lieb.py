#!/usr/bin/env python
import sys

class Lieb:
    W = 1
    L = 1
    Wsub = 1
    Lsub = 1
    N = 3

    def __init__(self, W, L, Wsub, Lsub):
        self.W = W
        self.L = L
        self.Wsub = Wsub
        self.Lsub = Lsub
        self.N = W * L * 3

    def trans_sites(self):
        pairs = []
        W = self.W
        L = self.L
        
        for xi in range(W):
            for yi in range(L):
                isite0 = (xi + yi * W) * 3
                pairs.append((isite0, isite0 + 1))
                pairs.append((isite0, isite0 + 2))

                yj = yi
                xj = (xi + 1) % W
                jsite0 = (xj + yj * W) * 3
                pairs.append((isite0 + 1, jsite0))

                xj = xi
                yj = (yi + 1) % L
                jsite0 = (xj + yj * W) * 3
                pairs.append((isite0 + 2, jsite0))
        return pairs

    def orbitalidx(self):
        idx = {}
        nidx = 0
        W = self.W
        L = self.L
        Wsub = self.Wsub
        Lsub = self.Lsub
        N = self.N

        for xi in range(W):
            for yi in range(L):
                for ki in range(3):
                    i = (xi + yi * W) * 3 + ki
                    for xj in range(W):
                        for yj in range(L):
                            for kj in range(3):
                                j = (xj + yj * W) * 3 + kj
                                if not ((i, j) in idx.keys()):
                                    for xp in range(W / Wsub):
                                        for yp in range(L / Lsub):
                                            ip = (((xi + Wsub * xp) % W) + \
                                                  ((yi + Lsub * yp) % W) * W) * 3 + ki
                                            jp = (((xj + Wsub * xp) % W) + \
                                                  ((yj + Lsub * yp) % W) * W) * 3 + kj
                                            idx[(ip, jp)] = nidx
                                    nidx += 1
        return (nidx, idx)

    def jastrowidx(self):
        idx = {}
        nidx = 0
        W = self.W
        L = self.L
        Wsub = self.Wsub
        Lsub = self.Lsub
        N = self.N

        for xi in range(W):
            for yi in range(L):
                for ki in range(3):
                    i = (xi + yi * W) * 3 + ki
                    for j in range(i + 1, N):
                        kj = j % 3
                        rj = j //3
                        xj =rj % W
                        yj =rj //W
                        if not ((i, j) in idx.keys()):
                            for xp in range(W / Wsub):
                                for yp in range(L / Lsub):
                                    ip = (((xi + Wsub * xp) % W) + \
                                          ((yi + Lsub * yp) % W) * W) * 3 + ki
                                    jp = (((xj + Wsub * xp) % W) + \
                                          ((yj + Lsub * yp) % W) * W) * 3 + kj
                                    idx[(ip, jp)] = nidx
                                    idx[(jp, ip)] = nidx
                            nidx += 1
        return (nidx, idx)

    def momprojidx(self):
        idx = []
        W = self.W
        L = self.L
        Wsub = self.Wsub
        Lsub = self.Lsub

        for xr in range(Wsub):
            for yr in range(Lsub):
                proj = list()
                for xi in range(W):
                    for yi in range(L):
                        for k in range(3):
                            i = (xi + yi * W) * 3 + k
                            xj = (xi + xr) % W
                            yj = (yi + yr) % L
                            j = (xj + yj * W) * 3 + k
                            proj.append((i, j))
                idx.append(proj)
        return idx

if len(sys.argv) < 7:
    print("Usage: Gen [options] {W} {L} {Wsub} {Lsub}, {t}, {J}")
    print(" Options")
    print("  -a: apply periodic-antiperiodic boundary condition (not ready)")
    print("  -t: write legacy input for tensor network version of mVMC instead")
    sys.exit()

flag_ap = 0
flag_tn = 0
for i in range(1, len(sys.argv)):
    if sys.argv[i][0] == '-':
        switch = sys.argv.pop(i)
        for c in switch:
            if c == 'a':
                flag_ap = 1
            if c == 't':
                flag_tn = 1
        break

latt = Lieb(int(sys.argv[1]), \
            int(sys.argv[2]), \
            int(sys.argv[3]), \
            int(sys.argv[4]))
vart = float(sys.argv[5])
varJ = float(sys.argv[6])

trans = latt.trans_sites()
(norb, orbidx) = latt.orbitalidx()
(ncor, coridx) = latt.jastrowidx()
qptidx = latt.momprojidx()

fout = open('orbitalidx.def', 'w')
fout.write("====================\n")
fout.write(" NOrb    %d         \n" % norb)
fout.write(" Complex 0          \n")
fout.write("====================\n")
fout.write("====================\n")
keys = list(orbidx.keys())
keys.sort()
for orb in keys:
    fout.write(("%8d %8d        " % orb) + str(orbidx[orb]) + '\n')
for i in range(norb):
    fout.write(str(i) + " 1\n")
fout.close()

fout = open('jastrowidx.def', 'w')
fout.write("====================\n")
fout.write(" Jastrow %d         \n" % ncor)
fout.write(" Complex 0          \n")
fout.write("====================\n")
fout.write("====================\n")
keys = list(coridx.keys())
keys.sort()
for orb in keys:
    fout.write(("%8d %8d        " % orb) + str(coridx[orb]) + '\n')
for i in range(ncor):
    fout.write(str(i) + " 1\n")
fout.close()

fout = open('qptransidx.def', 'w')
fout.write("====================\n")
fout.write(" NQpt    %d         \n" % len(qptidx))
fout.write("====================\n")
fout.write("====================\n")
fout.write("====================\n")
for i in range(len(qptidx)):
    fout.write(str(i) + " 1.0000\n")
for i in range(len(qptidx)):
    for rq in qptidx[i]:
        fout.write(str(i) + " %d %d\n" % rq)
fout.close()

fout = open('transfer.def', 'w')
fout.write("====================\n")
fout.write(" NTrans   %d        \n" % (len(trans) * 4))
fout.write("====================\n")
fout.write("====================\n")
fout.write("====================\n")
for nn in trans:
    fout.write(("%d 0 %d 0 " % nn) + str(vart) + '\n')
    fout.write(("%d 1 %d 1 " % nn) + str(vart) + '\n')
    fout.write(("%d 0 %d 0 " % (nn[1], nn[0])) + str(vart) + '\n')
    fout.write(("%d 1 %d 1 " % (nn[1], nn[0])) + str(vart) + '\n')
fout.close()

fout = open('integralJ.def', 'w')
fout.write("====================\n")
fout.write(" NInt     %d        \n" % len(trans))
fout.write("====================\n")
fout.write("====================\n")
fout.write("====================\n")
for nn in trans:
    fout.write(("%d %d " % nn) + str(varJ) + '\n')
fout.close()


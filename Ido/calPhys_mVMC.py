#!/usr/bin/env python
import sys, glob
from math import *
import numpy as np

T1={(1,0): 1.0, (-1,0): 1.0, (0,1): -1.0, (0,-1): -1.0
}
Tx={(1,0): 1.0, (-1,0): 1.0
}
Ty={(0,1): 1.0, (0,-1): 1.0
}

def position(r):
    x=(r[0]+Lx)%Lx
    y=(r[1]+Ly)%Ly
    return (x,y)

def InvVec(r):
    return (-r[0],-r[1])

def indexToPosition(i):
    x=i%Lx
    y=i/Lx
    return (x,y)

def positionToIndex(r):
    x=(r[0]+Lx)%Lx
    y=(r[1]+Ly)%Ly
    return x+Lx*y

def neighborIndex(i,dr):
    r=indexToPosition(i)
    x=r[0]+dr[0]
    y=r[1]+dr[1]
    return positionToIndex([x,y])

def direction(i,j):
    ri=indexToPosition(i)
    rj=indexToPosition(j)
    dx=(rj[0]-ri[0]+Lx)%Lx
    dy=(rj[1]-ri[1]+Ly)%Ly
    return (dx,dy)

def locgrnIdx(i,j,s):
    return (Nsite*i+j)+s*Nsite*Nsite

def VecSum(ri,rj):
    return (ri[0]+rj[0], ri[1]+rj[1]) 

def VecInv(ri):
    return (-ri[0], -ri[1]) 

def VecDiff(ri,rj):
    return (ri[0]-rj[0], ri[1]-rj[1]) 

def subIndex(i):
    r=indexToPosition(i)
    sx=r[0]%Sx
    sy=r[1]%Sy
    return sx+Sx*sy

def subIndexToIndex(s):
    sx = s%Sx 
    sy = s/Sx 
    i = positionToIndex([sx,sy])
    return i 
  
def distance(r):
    # the square of the distance between site i and j
    rx=r[0]
    ry=r[1]
    if rx>Lx/2:
        while rx>Lx/2:
            rx -= Lx
    if ry>Ly/2:
        while ry>Ly/2:
            ry -= Ly
    return sqrt(rx*rx+ry*ry)

def sgnAP(i,dr):
    if APFlag == 0:
      return 1.0
    else:
      r=indexToPosition(i)
      x=r[0]+dr[0]
      if (x < Lx-0.5 and x > -0.5):
          return 1.0
      else:
          return -1.0
      #y=r[1]+dr[1]
      #if (y < Ly-0.5 and y > -0.5):
      #  return 1.0
      #else:
      #  return -1.0
      
def WaveNumberList(Lx,Ly,APFlag):
  wave = [] 
  if(APFlag == 1):  
    for mx in xrange(-Lx/2,Lx/2):
      qx = 2.0*pi*float(0.5+mx)/float(Lx)
      for my in xrange(-Ly/2,Ly/2+1):
        qy = 2.0*pi*float(my)/float(Ly)
        wave.append((qx,qy))
  else:
    for mx in range(-Lx/2,Lx/2+1):
      qx = 2.0*pi*float(mx)/float(Lx)
      for my in xrange(-Ly/2,Ly/2+1):
        qy = 2.0*pi*float(my)/float(Ly)
        wave.append((qx,qy))
  return wave

def CosSinList(Lx,Ly,APFlag):
  cosList = {}
  sinList = {}
  Nsite = Lx*Ly
  for q in WaveNumberList(Lx,Ly,APFlag):
    (qx,qy) = q
    for ix in range(-Lx,Lx+1):
      for iy in range(-Ly,Ly+1):
        ri = (ix,iy)
        cosList[(q,ri)] = cos(qx*float(ix) + qy*float(iy))
        sinList[(q,ri)] = sin(qx*float(ix) + qy*float(iy))
      
      #for j in xrange(Nsite):
      #  jx = j%Lx
      #  jy = j/Lx
      #  cosList[(q,i,j)] = cos(qx*float(ix-jx) + qy*float(iy-jy))
      #  sinList[(q,i,j)] = sin(qx*float(ix-jx) + qy*float(iy-jy))
  return [cosList,sinList]

def ReadPara(para):
  ifile = open(para, 'r')
  for itr in xrange(5):
    ifile.readline().split()
  [name,DataFileHead] = ifile.readline().split()
  for itr in xrange(5):
    ifile.readline().split()
  [name,NDataIdxStart] = ifile.readline().split()
  [name,NDataQtySmp]   = ifile.readline().split()
  ifile.readline().split()
  [name,Nsite] = ifile.readline().split()
  [name,Ne] = ifile.readline().split()
  for itr in xrange(3):
    ifile.readline().split()
  [name,SRItrStep] = ifile.readline().split()
  [name,SRItrSmp]  = ifile.readline().split()
  for itr in xrange(3):
    ifile.readline().split()
  [name,SRStepDt] = ifile.readline().split()
  ifile.close()

  Nsite = int(Nsite)
  Lx = int(sqrt(Nsite))
  Ly = int(sqrt(Nsite))
  Ne = int(Ne)
  NDataIdxStart = int(NDataIdxStart)
  NDataQtySmp   = int(NDataQtySmp)
  SRItrStep = int(SRItrStep)
  SRItrSmp  = int(SRItrSmp)
  SRStepDt  = float(SRStepDt)

  return [DataFileHead,NDataIdxStart,NDataQtySmp,Nsite,Ne,Lx,Ly,SRItrStep,SRItrSmp,SRStepDt]

def GetFileList(flist):
  ifile = open(flist, 'r')
  lines = []
  for line in ifile:
    line = line.rstrip("\n")
    lines.append(line)
  ifile.close()
  return lines

def ReadCisAjs(flist):
  ifile = open(flist, 'r')
  ifile.readline().split()
  [name,NCisAjs] = ifile.readline().split()
  NCisAjs = int(NCisAjs)
  ciajList = []
  ciajListAppend = ciajList.append
  for itr in xrange(3):
    ifile.readline().split()
  for itr in xrange(NCisAjs):
    tmp = ifile.readline().split()
    tmp = map(int, tmp)
    [i,s,j,t] = tmp
    #[idx,i,s,j,t] = tmp
    ciajListAppend((i,s,j,t))
  ifile.close()
  return [NCisAjs, list(ciajList)]

def ReadCisAjsCktAlt(flist):
  ifile = open(flist, 'r')
  ifile.readline().split()
  [name,NCisAjsCktAlt] = ifile.readline().split()
  NCisAjsCktAlt = int(NCisAjsCktAlt)
  ciajckalList = []
  ciajckalListAppend = ciajckalList.append
  for itr in xrange(3):
    ifile.readline().split()
  for itr in xrange(NCisAjsCktAlt):
    tmp = ifile.readline().split()
    tmp = map(int, tmp)
    [idx,idx2,i,j,s,k,l,t] = tmp
    ciajckalListAppend((i,j,s,k,l,t))
  ifile.close()
  return [NCisAjsCktAlt, list(ciajckalList)]

def ReadCisAjsCktAltDC(flist):
  ifile = open(flist, 'r')
  ifile.readline().split()
  [name,NCisAjsCktAlt] = ifile.readline().split()
  NCisAjsCktAlt = int(NCisAjsCktAlt)
  ciajckalList = []
  ciajckalListAppend = ciajckalList.append
  for itr in xrange(3):
    ifile.readline().split()
  for itr in xrange(NCisAjsCktAlt):
    tmp = ifile.readline().split()
    tmp = map(int, tmp)
    [i,si,j,sj,k,sk,l,sl] = tmp
    ciajckalListAppend((i,si,j,sj,k,sk,l,sl))
  ifile.close()
  return [NCisAjsCktAlt, list(ciajckalList)]

################################################################
#
#   Main:Calculate Physical Quantities
#
################################################################

if len(sys.argv)<9:
    print "./calPhys_mVMC.py modepara.def greenOne.def greenTwo.def Mode DCFlag APFlag NkFlag SCFlag[num..]"
    #print "./average.py xnamelist.def Mode DCFlag APFlag NkFlag [num..]"
    print "usage: Mode == 0 -> For Ground states"
    print "       Mode == 1 -> For Nonequilibrium states\n"
    print "       DCFlag==0 -> Calculate two-body green functions by using one-body green functions"
    print "       DCFlag==1 -> Calculate two-body green functions directly\n"
    print "       APFlag==0 -> Periodic-Periodic boundary condition"
    print "       APFlag==1 -> AntiPeriodic-Periodic boundary condition\n"
    print "       NkFlag==0 -> Not Calculate momentum distribution n(k)"
    print "       NkFlag==1 -> Calculate momentum distribution n(k)\n"
    print "       SCFlag==0 -> Not Calculate Supercoducting correlations Pd(r)"
    print "       SCFlag==1 -> Calculate Superconducting correlations Pd(r)\n"
    sys.exit()

Mode   = int(sys.argv[4])
DCFlag = int(sys.argv[5])
APFlag = int(sys.argv[6])
NkFlag = int(sys.argv[7])
SCFlag = int(sys.argv[8])
    
Sx = 2
Sy = Sx
Nsub = Sx*Sy 

#fileList = GetFileList(sys.argv[1])
paraList = ReadPara(sys.argv[1])
[NCisAjs, caList] = ReadCisAjs(sys.argv[2])
#if DCFlag == 0:
#  [NCisAjsCktAlt, cacaList] = ReadCisAjsCktAlt(sys.argv[3])
#else:
#  [NCisAjsCktAlt, cacaList] = ReadCisAjsCktAltDC(sys.argv[3])
[NCisAjsCktAlt, cacaList] = ReadCisAjsCktAltDC(sys.argv[3])

pre = paraList[0]
Nsite = paraList[3]
Ne = paraList[4]
Lx = paraList[5]
Ly = paraList[6]
if Mode == 0:
  Ntime = 1
else:
  Ntime = paraList[7]
dt = paraList[9]
U = 1.0
pre2 = 'phy'
pre3 = 'ave'

WaveVector = WaveNumberList(Lx, Ly, APFlag)
[CosList, SinList] = CosSinList(Lx,Ly,APFlag)

if len(sys.argv)>9:
    numList = sys.argv[9:]
else:
    numList = [ n[len(pre)+5:len(pre)+8] for n in glob.glob(pre+"_out_*.dat")]
numList.sort()
nNum = len(numList)

for num in numList:
# Read Green Functions
    ca    = open(pre+'_cisajs_'        +num+'.dat','r')
    #if DCFlag == 0:
    #  caca  = open(pre+'_cisajscktalt_'  +num+'.dat','r')
    #else:
    #  caca  = open(pre+'_cisajscktaltdc_'+num+'.dat','r')
    caca  = open(pre+'_cisajscktalt_'+num+'.dat','r')

    time = [[] for i in xrange(Ntime)]
    CisAjs = {}
    CisAjsCktAlt = {}
 
    for k in xrange(Ntime):
        time[k] = (k+1)*dt 
        data = ca.readline().split()
        # read greenFunc
        for idx in xrange(NCisAjs):
          CisAjs[caList[idx]] = float(data[idx])
        data = caca.readline().split()
        # read greenFunc2
        for idx in xrange(NCisAjsCktAlt):
          CisAjsCktAlt[cacaList[idx]] = float(data[idx])

    ca.close()
    caca.close()

# Cal Physical Quantities
    # Local Density
    ofile = open(pre2+'_local_'+num+'.dat','w')
    ofile.write("# i  j  S_iS_j\n")
    for rx in xrange(Lx):
      for ry in xrange(Ly):
        i = rx + ry*Lx
        ofile.write("{0}\t{1}\t{2: .18f}\t{3: .18f}\n".format(rx, ry, CisAjs[(i,0,i,0)].real, CisAjs[(i,1,i,1)].real))
    ofile.close()

    # Spin Correlation
    print num, "ss"
    ssdc = [[[[] for j in range(Nsite)] for i in xrange(Nsite)] for k in xrange(Ntime)]
    ofile = open(pre2+'_ss_'+num+'.dat','w')
    ofile.write("# i  j  S_iS_j\n")
    for k in xrange(Ntime):
      for i in xrange(Nsite):
        for j in xrange(Nsite):
        #for jsub in xrange(Nsub):
        #  j = subIndexToIndex(jsub)
          s = 0.0
          if(i==j):
              s += 0.5*(CisAjs[(i,0,i,0)] + CisAjs[(i,1,i,1)])
          s += (-0.5)*(CisAjsCktAlt[(i,0,j,0,j,1,i,1)]+CisAjsCktAlt[(j,0,i,0,i,1,j,1)])
          s += 0.25*(CisAjsCktAlt[(i,0,i,0,j,0,j,0)] + (CisAjsCktAlt[(i,1,i,1,j,1,j,1)]))
          s -= 0.25*(CisAjsCktAlt[(i,0,i,0,j,1,j,1)] + (CisAjsCktAlt[(i,1,i,1,j,0,j,0)]))
          ofile.write("{0}\t{1}\t{2: .18f}\t{3: .18f}\n".format(i,j,s.real,s.imag))
          ssdc[k][i][j] = s
    ofile.close()

    # Spin Structure Function
    print num, "sq"
    ofile = open(pre2+'_sq_'+num+'.dat','w')
    ofile.write("# qx/pi qy/pi S(q).real  S(q).imag\n")
    for k in xrange(Ntime):
      mat = [[0.0 for rx in xrange(Lx)] for ry in xrange(Ly)]
      for rx in xrange(Lx):
        for ry in xrange(Ly):
            for j in xrange(Nsite):
            #for jsub in xrange(Nsub):
            #    j = subIndexToIndex(jsub)
                jx = j%Lx
                jy = j/Lx
                ix = (jx + rx)%Lx
                iy = (jy + ry)%Ly
                i = ix + iy*Lx
                mat[rx][ry] += ssdc[k][i][j]
      sq = np.fft.fft2(mat)
      for mx in xrange(Lx+1):
        qx = 2.0*pi*float(mx)/float(Lx)
        for my in xrange(Ly+1):
            qy = 2.0*pi*float(my)/float(Ly)

            sqReal = sq[mx%Lx][my%Ly].real/ (3.0*float(Nsite))
            sqImag = sq[mx%Lx][my%Ly].imag/ (3.0*float(Nsite))
            #sqReal = sq[mx%Lx][my%Ly].real/ (3.0*float(Nsub))
            #sqImag = sq[mx%Lx][my%Ly].imag/ (3.0*float(Nsub))
            ofile.write("{0}\t{1}\t{2: .18e}\t{3: .18e}\n".format(qx/pi,qy/pi,sqReal,sqImag))
        ofile.write("\n")
    ofile.close()

    # Momentum Distribution
    print num, "nk"
    if(NkFlag == 1):
      ofile = open(pre2+'_nk_'+num+'.dat','w')
      ofile.write("# kx/pi  ky/pi  n(k).real  n(k).imag\n")
      for k in xrange(Ntime):
        for q in WaveVector:
          (qx,qy) = q
          nkReal = 0.0
          nkImag = 0.0
          for i in xrange(Nsite):
            ri = indexToPosition(i)
            for j in xrange(Nsite):
              rj = indexToPosition(j)
              rij = VecDiff(ri,rj)

              nkReal += CosList[(q,rij)] * (CisAjs[(i,0,j,0)] + CisAjs[(i,1,j,1)])
              nkImag += SinList[(q,rij)] * (CisAjs[(i,0,j,0)] + CisAjs[(i,1,j,1)])
          nkReal /= 2.0*float(Nsite)
          nkImag /= 2.0*float(Nsite)
          ofile.write("{0: .10f}\t{1: .10f}\t{2: .18e}\t{3: .18e}\n".\
                        format(qx/pi,qy/pi,nkReal.real,nkImag.imag))
        ofile.write("\n")
      ofile.close()
 
    # Charge Structure Factor 
    print num, "nq"
    ofile = open(pre2+'_nq_'+num+'.dat','w')
    ofile.write("# qx qy  N(q).real  N(q).imag\n")
    for k in range(Ntime):
      mat = [[0.0 for rx in xrange(Lx)] for ry in xrange(Ly)]
      for rx in xrange(Lx):
        for ry in xrange(Ly):
            for j in xrange(Nsite):
            #for jsub in xrange(Nsub):
            #    j = subIndexToIndex(jsub)
                jx = j%Lx
                jy = j/Lx
                ix = (jx + rx)%Lx
                iy = (jy + ry)%Ly
                i = ix + iy*Lx
                mat[rx][ry] += CisAjsCktAlt[(i,0,i,0,j,0,j,0)] + CisAjsCktAlt[(i,1,i,1,j,1,j,1)]
                mat[rx][ry] += CisAjsCktAlt[(i,0,i,0,j,1,j,1)] + CisAjsCktAlt[(i,1,i,1,j,0,j,0)]
                mat[rx][ry] -= (2.0*Ne/Nsite)*(2.0*Ne/Nsite)  
      nq = np.fft.fft2(mat)

      for mx in xrange(Lx+1):
        qx = 2.0*pi*float(mx)/float(Lx)
        for my in xrange(Ly+1):
            qy = 2.0*pi*float(my)/float(Ly)

            nqReal = nq[mx%Lx][my%Ly].real / float(Nsite)
            nqImag = nq[mx%Lx][my%Ly].imag / float(Nsite)
            ofile.write("{0: .10f}\t{1: .10f}\t{2: .18e}\t{3: .18e}\n".\
                            format(qx/pi,qy/pi,nqReal,nqImag))
        ofile.write("\n")
    ofile.close()

    # Double Occupancy
    print num, "do"
    ofile = open(pre2+'_do_'+num+'.dat','w')
    ofile.write("# nn\n")
    for k in xrange(Ntime):
      nn = 0.0
      for isub in xrange(Nsub):
        i = subIndexToIndex(isub)
        nn += CisAjsCktAlt[(i,0,i,0,i,1,i,1)]
      nn /= float(Nsub)
      ofile.write("{0: .18e}\t{1: .18e}\n".format(nn.real,nn.imag))
    ofile.close()


    # Super Conducting Correlation Pd(ro)
    if SCFlag == 1:
     print num, "sc"
     ofile = open(pre2+'_sc_'+num+'.dat','w')
     ofile.write("# i  j  Delta+_i_Delta_j\n")
     ofile2 = open(pre2+'_sc2_'+num+'.dat','w')
     ofile2.write("# i  j  Delta+_i_Delta_j\n")
     for t in xrange(Ntime):
      mat = [[0.0 for rx in xrange(Lx)] for ry in xrange(Ly)]
      SCmat = {}
      for rx in xrange(Lx):
        for ry in xrange(Ly):
         ro = [rx, ry]
         o = positionToIndex(ro)
          #i=0
         scAll = 0.0
         #for isub in xrange(Nsub):#summation ri
         # i = subIndexToIndex(isub)
         for i in xrange(Nsite):#ri summation
         #for i in xrange(50,51):#ri summation
          sc = 0.0
          io = neighborIndex(i,ro)
          sgn_io=sgnAP(i,ro)
          for dr_i,sgn_form_i in T1.iteritems():
            i_dr = neighborIndex(i,dr_i)
            sgn_i_dr=sgnAP(i,dr_i)
            for dr_io,sgn_form_io in T1.iteritems():
              io_dr = neighborIndex(io,dr_io)
              ro_dr = [rx+dr_io[0], ry+dr_io[1]]
              sgn_io_dr = sgnAP(i,ro_dr)
              sgn_dsc   = sgn_form_i*sgn_form_io
              sgn_ap    = sgn_io*sgn_i_dr*sgn_io_dr
              sgn       = sgn_dsc*sgn_ap
              #sgn = sgn_dsc

              sc += CisAjsCktAlt[(i,    0, io,    0, i_dr, 1, io_dr, 1)]*sgn
              sc += CisAjsCktAlt[(i,    0, io_dr, 0, i_dr, 1, io,    1)]*sgn
              sc += CisAjsCktAlt[(i_dr, 0, io,    0, i,    1, io_dr, 1)]*sgn
              sc += CisAjsCktAlt[(i_dr, 0, io_dr, 0, i,    1, io,    1)]*sgn
              
              sc += CisAjsCktAlt[(io,    0, i,    0, io_dr, 1, i_dr, 1)]*sgn
              sc += CisAjsCktAlt[(io,    0, i_dr, 0, io_dr, 1, i,    1)]*sgn
              sc += CisAjsCktAlt[(io_dr, 0, i,    0, io,    1, i_dr, 1)]*sgn
              sc += CisAjsCktAlt[(io_dr, 0, i_dr, 0, io,    1, i,    1)]*sgn

              if(i == io):
                sc -= (CisAjs[(io_dr, 0, i_dr,0)] + CisAjs[(io_dr, 1, i_dr, 1)])*sgn
                if(i_dr == io_dr):
                  sc += 2.0*sgn
              if(i_dr == io_dr):
                sc -= (CisAjs[(io, 0, i, 0)] + CisAjs[(io, 1, i, 1)])*sgn
              if(i == io_dr):
                sc -= (CisAjs[(io, 0, i_dr,0)] + CisAjs[(io, 1, i_dr, 1)])*sgn
                if(io == i_dr):
                  sc += 2.0*sgn
              if(io == i_dr):
                sc -= (CisAjs[(io_dr, 0, i,0)] + CisAjs[(io_dr, 1, i,1)])*sgn
          SCmat[(rx,ry,i)] = sc/4.0
          scAll = scAll + sc
         #mat[rx][ry] = scAll/(4.0)
         #mat[rx][ry] += scAll/(4.0*Nsub)
         mat[rx][ry] += scAll/(4.0*Nsite)
      maxP = {}
      for rx in xrange(Lx):
        for ry in xrange(Ly):
          r = distance((rx,ry))
          if not r in maxP:
            maxP[r] = 0.0
          if fabs(mat[ry][rx]) > fabs(maxP[r]):
            #maxP[r] = fabs(mat[ry][rx])
            maxP[r] = (mat[ry][rx])
      numLine = 0
      for rx in xrange(Lx):
        for ry in xrange(Ly):
          r = distance((rx,ry))
          val = maxP[r]
          if(r <= sqrt(2.0)*float(Lx)*0.5):
            ofile.write("{0: .10e}\t{1: .10e}\t{2: .10e}\n".format(r,val.real,val.imag))
            numLine = numLine + 1
      maxP = {}
      for ix in xrange(Lx):
        for iy in xrange(Ly):
          i = positionToIndex((ix,iy))
          for rx in xrange(Lx):
            for ry in xrange(Ly):
              r = distance((rx,ry))
              if not (i,rx,ry) in maxP:
                maxP[(i,rx,ry)] = 0.0
              if fabs(SCmat[(rx,ry,i)]) > fabs(maxP[(i,rx,ry)]):
                #maxP[(i,r)] = fabs(SCmat[(rx,ry,i)])
                maxP[(i,rx,ry)] = (SCmat[(rx,ry,i)])
          for rx in xrange(Lx):
            for ry in xrange(Ly):
              val = maxP[(i,rx,ry)]
              #if(r <= sqrt(2.0)*float(Lx)*0.5):
              #if(sqrt(2.0)*float(Lx)*0.5-0.1 < r <= sqrt(2.0)*float(Lx)*0.5):
              ofile2.write("{0}\t{1}\t{2}\t{3}\t{4: .10e}\t{5: .10e}\n".format(ix,iy,rx,ry,val.real,val.imag))
     ofile.close()
     ofile2.close()

##############################
# Make Average Files
##############################
ifile = [[] for i in xrange(nNum)]
ifile2 = [[] for i in xrange(nNum)]
# energy
print "ave", "ene"
for i in range(nNum):
    ifile[i] = open(pre+"_out_"+numList[i]+".dat",'r')
ofile = open(pre3+"_ene.dat",'w')
ofile.write("# E  sigma\n")
for k in xrange(Ntime):
    sum = [0.0, 0.0]
    for i in xrange(nNum):
      data = ifile[i].readline().split()
      time[k] = float(data[0])
      #dtmp = float(data[1])/float(Nsite)
      dtmp = float(data[0])/float(Nsite)
      sum[0] += dtmp
      sum[1] += dtmp*dtmp
    ave = sum[0]/float(nNum)
    var = sum[1]/float(nNum) - ave*ave
    sigma = sqrt(abs(var/float(nNum)))
    
    ofile.write("{0: .10f}\t{1: .10f}\t{2: .10f}\n".format(time[k],ave,sigma))
for i in xrange(nNum):
    ifile[i].close()
ofile.close()

# local charge (spin) density
print "ave", "local"
for i in xrange(nNum):
    ifile[i] = open(pre2+"_local_"+numList[i]+".dat",'r')
    data = ifile[i].readline()
ofile = open(pre3+"_local.dat",'w')
ofile2 = open(pre3+"_local_x.dat",'w')
ofile.write("# i j S_iS_j sigma\n")
ave_up = 0.0
ave_dw = 0.0
spin = 0.0
sigma_up = 0.0
sigma_dw = 0.0
for k in xrange(Ntime):
  val = [[[] for i in xrange(Ly+1)]\
          for j in xrange(4)]
  for j in xrange(Nsite):
    sum = [0.0, 0.0]
    sum2 = [0.0, 0.0]
    for i in xrange(nNum):
        data = ifile[i].readline().split()
        itmp = int(data[0])
        jtmp = int(data[1])
        dtmp = float(data[2])
        dtmp2 = float(data[3])
        sum[0] += dtmp
        sum[1] += dtmp*dtmp
        sum2[0] += dtmp2
        sum2[1] += dtmp2*dtmp2
    ave = sum[0]/float(nNum)
    var = sum[1]/float(nNum) - ave*ave
    ave2 = sum2[0]/float(nNum)
    var2 = sum2[1]/float(nNum) - ave2*ave2
    sigma = sqrt(abs(var/float(nNum)))
    sigma2 = sqrt(abs(var2/float(nNum)))
    ofile.write("{0}\t{1}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\n".format(itmp,jtmp,ave,sigma,ave2,sigma2))
    if(jtmp != Ly-1):
      ave_up += ave
      ave_dw += ave2
      spin += abs((ave-ave2))
      #spin += (ave-ave2)*(-1)**(jtmp)
      sigma_up += sigma*sigma
      sigma_dw += sigma2*sigma2
    else:
      #ofile2.write("{0}\t{1: .10e}\t{2: .10e}\t{3: .10e}\n".format(itmp,ave_up/Ly+ave_dw/Ly,ave_up/Ly-ave_dw/Ly,sqrt(sigma_up+sigma_dw)/Ly))
      ofile2.write("{0}\t{1: .10e}\t{2: .10e}\t{3: .10e}\n".format(itmp,ave_up/Ly+ave_dw/Ly,spin/Ly,sqrt(sigma_up+sigma_dw)/Ly))
      ave_up = 0.0
      ave_dw = 0.0
      spin = 0.0
      sigma_up = 0.0
      sigma_dw = 0.0
    if (jtmp == Ly-1):
      ofile.write("{0}\t{1}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\n".format(itmp,jtmp+1,ave,sigma,ave2,sigma2))
      #if (itmp == Lx-1):
      #  ofile.write("{0}\t{1}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\n".format(itmp+1,jtmp+1,ave,sigma,ave2,sigma2))
      ofile.write("\n")
    if (itmp == Lx-1):
      val[0][jtmp] = ave
      val[1][jtmp] = sigma
      val[2][jtmp] = ave2
      val[3][jtmp] = sigma2
      if(jtmp == Ly-1):
        val[0][jtmp+1] = ave
        val[1][jtmp+1] = sigma
        val[2][jtmp+1] = ave2
        val[3][jtmp+1] = sigma2
  if (itmp == Lx-1):
    for j in xrange(Ly+1):
      ofile.write("{0}\t{1}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\n".format(itmp+1,j,val[0][j],val[1][j],val[2][j],val[3][j]))
  ofile.write("\n")
for i in xrange(nNum):
    ifile[i].close()
ofile.close()
ofile2.close()


# spin correlation ver.dc
print "ave", "ss"
for i in xrange(nNum):
    ifile[i] = open(pre2+"_ss_"+numList[i]+".dat",'r')
    data = ifile[i].readline()
ofile = open(pre3+"_ss.dat",'w')
ofile.write("# i j S_iS_j sigma\n")
for k in xrange(Ntime):
  for j in xrange(Nsite*Nsub):
    sum = [0.0, 0.0]
    for i in xrange(nNum):
        data = ifile[i].readline().split()
        itmp = int(data[0])
        jtmp = int(data[1])
        dtmp = float(data[2])
        sum[0] += dtmp
        sum[1] += dtmp*dtmp
    ave = sum[0]/float(nNum)
    var = sum[1]/float(nNum) - ave*ave
    sigma = sqrt(abs(var/float(nNum)))
    ofile.write("{0}\t{1}\t{2: .10e}\t{3: .10e}\n".format(itmp,jtmp,ave,sigma))
    ofile.write("\n\n")
for i in xrange(nNum):
    ifile[i].close()
ofile.close()

# spin structure factor ver. dc
print "ave", "sq"
for i in xrange(nNum):
    ifile[i] = open(pre2+"_sq_"+numList[i]+".dat",'r')
    data = ifile[i].readline()
ofile = open(pre3+"_sq.dat",'w')
ofile.write("# qx/pi qy/pi  S(q).real sigma  S(q).imag sigma\n")
ofile2 = open(pre3+"_sq_t.dat",'w')
ofile2.write("# time  S(pi,0).real sigma  S(pi,0).imag sigma S(pi,pi).real sigma S(pi,pi).imag sigma\n")
for k in xrange(Ntime):
  for mx in xrange(Lx+1):
    for my in xrange(Ly+1):
        sum = [0.0, 0.0, 0.0, 0.0]
        for i in xrange(nNum):
            data = ifile[i].readline().split()
            dtmpx = float(data[0])
            dtmpy = float(data[1])
            dtmp1 = float(data[2])
            dtmp2 = float(data[3])
            sum[0] += dtmp1
            sum[1] += dtmp1*dtmp1
            sum[2] += dtmp2
            sum[3] += dtmp2*dtmp2
        ave1 = sum[0]/float(nNum)
        var1 = sum[1]/float(nNum) - ave1*ave1
        sigma1 = sqrt(abs(var1/float(nNum)))
        ave2 = sum[2]/float(nNum)
        var2 = sum[3]/float(nNum) - ave2*ave2
        sigma2 = sqrt(abs(var2/float(nNum)))
        ofile.write("{0:.10f}\t{1:.10f}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\n"\
                        .format(dtmpx,dtmpy,ave1,sigma1,ave2,sigma2))
        if((dtmpx==1.0) and (dtmpy==0.0)):
            ofile2.write("{0:.10e}\t{1: .10e}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t"\
                        .format(time[k],ave1,sigma1,ave2,sigma2))
        if((dtmpx==1.0) and (dtmpy==1.0)):
            ofile2.write("{0:.10e}\t{1: .10e}\t{2: .10e}\t{3: .10e}\n"\
                        .format(ave1,sigma1,ave2,sigma2))
    for i in xrange(nNum):
        ifile[i].readline()
    ofile.write("\n")
  ofile.write("\n\n")
for i in xrange(nNum):
    ifile[i].close()
ofile.close()
ofile2.close()

# momentum distribution
if(NkFlag == 1):
  print "ave", "nk"
  for i in range(nNum):
    ifile[i] = open(pre2+"_nk_"+numList[i]+".dat",'r')
    data = ifile[i].readline()
  ofile = open(pre3+"_nk.dat",'w')
  ofile.write("# kx/pi ky/pi  n(k).real sigma  n(k).imag sigma\n")
  ofile2 = open(pre3+"_dnk.dat",'w')
  ofile2.write("# time  dnk(pi,0).real sigma\n")
  for k in xrange(Ntime):
    if (APFlag == 0):
       xList = xrange(-Lx/2,Lx/2+1)
    else:
       xList = xrange(-Lx/2,Lx/2)
    for mx in xList:
      for my in xrange(-Ly/2,Ly/2+1):
        sum = [0.0, 0.0, 0.0, 0.0]
        for i in xrange(nNum):
            data = ifile[i].readline().split()
            dtmpx = float(data[0])
            dtmpy = float(data[1])
            dtmp1 = float(data[2])
            dtmp2 = float(data[3])
            sum[0] += dtmp1
            sum[1] += dtmp1*dtmp1
            sum[2] += dtmp2
            sum[3] += dtmp2*dtmp2
        ave1 = sum[0]/float(nNum)
        var1 = sum[1]/float(nNum) - ave1*ave1
        sigma1 = sqrt(abs(var1/float(nNum)))
        ave2 = sum[2]/float(nNum)
        var2 = sum[3]/float(nNum) - ave2*ave2
        sigma2 = sqrt(abs(var2/float(nNum)))
        ofile.write("{0:.10f}\t{1:.10f}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\n"\
                        .format(dtmpx,dtmpy,ave1,sigma1,ave2,sigma2))
        if mx==int(Lx/2-1) and my==int(0):
          ave3 = ave1
          sigma3 = sigma1
        if mx==int(Lx/2-1) and my==int(1):
          ave4 = ave1
          sigma4 = sigma1
      ofile.write("\n")
    for i in xrange(nNum):
      ifile[i].readline()
    ofile.write("\n\n")
    ofile2.write("{0:.10f}\t{1:.10f}\t{2: .10e}\n"\
                        .format(time[k],ave3-ave4,sqrt(sigma3*sigma3+sigma4*sigma4)))
  for i in xrange(nNum):
    ifile[i].close()
  ofile.close()
  ofile2.close()

# charge structure factor ver. dc
print "ave", "nq"
for i in xrange(nNum):
    ifile[i] = open(pre2+"_nq_"+numList[i]+".dat",'r')
    data = ifile[i].readline()
ofile = open(pre3+"_nq.dat",'w')
ofile.write("# qx/pi qy/pi  N(q).real sigma  N(q).imag sigma\n")
ofile2 = open(pre3+"_nq_t.dat",'w')
ofile2.write("# time  N(pi,0).real sigma  N(pi,0).imag sigma N(pi,pi).real sigma N(pi,pi).imag sigma\n")
for k in xrange(Ntime):
  for mx in xrange(Lx+1):
    for my in xrange(Ly+1):
        sum = [0.0, 0.0, 0.0, 0.0]
        for i in xrange(nNum):
            data = ifile[i].readline().split()
            dtmpx = float(data[0])
            dtmpy = float(data[1])
            dtmp1 = float(data[2])
            dtmp2 = float(data[3])
            sum[0] += dtmp1
            sum[1] += dtmp1*dtmp1
            sum[2] += dtmp2
            sum[3] += dtmp2*dtmp2
        ave1 = sum[0]/float(nNum)
        var1 = sum[1]/float(nNum) - ave1*ave1
        sigma1 = sqrt(abs(var1/float(nNum)))
        ave2 = sum[2]/float(nNum)
        var2 = sum[3]/float(nNum) - ave2*ave2
        sigma2 = sqrt(abs(var2/float(nNum)))
        ofile.write("{0:.10f}\t{1:.10f}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t{5: .10e}\n"\
                        .format(dtmpx,dtmpy,ave1,sigma1,ave2,sigma2))
        if((dtmpx==1.0) and (dtmpy==0.0)):
            ofile2.write("{0:.10e}\t{1: .10e}\t{2: .10e}\t{3: .10e}\t{4: .10e}\t"\
                        .format(time[k],ave1,sigma1,ave2,sigma2))
        if((dtmpx==1.0) and (dtmpy==1.0)):
            ofile2.write("{0:.10e}\t{1: .10e}\t{2: .10e}\t{3: .10e}\n"\
                        .format(ave1,sigma1,ave2,sigma2))
    for i in xrange(nNum):
        ifile[i].readline()
    ofile.write("\n")
  ofile.write("\n\n")
for i in xrange(nNum):
    ifile[i].close()
ofile.close()
ofile2.close()


# double occupancy
print "ave", "do"
for i in xrange(nNum):
    ifile[i] = open(pre2+"_do_"+numList[i]+".dat",'r')
    data = ifile[i].readline()
ofile = open(pre3+"_do.dat",'w')
ofile.write("# nn sigma\n")
for k in xrange(Ntime):
  sum = [0.0, 0.0]
  for i in xrange(nNum):
    data = ifile[i].readline().split()
    dtmp = float(data[0])
    sum[0] += dtmp
    sum[1] += dtmp*dtmp
  ave = sum[0]/float(nNum)
  var = sum[1]/float(nNum) - ave*ave
  sigma = sqrt(abs(var/float(nNum)))
  ofile.write("{0: .10f}\t{1: .10f}\t{2: .10f}\n".format(time[k],ave,sigma))
for i in xrange(nNum):
    ifile[i].close()
ofile.close()


if SCFlag==1:
# super conducting correlation ver.dc
 print "ave", "sc"
 for i in xrange(nNum):
    ifile[i] = open(pre2+"_sc_"+numList[i]+".dat",'r')
    data = ifile[i].readline()
 ofile = open(pre3+"_sc.dat",'w')
 ofile.write("# r P_dsc(r) sigma\n")
 for k in xrange(Ntime):
    phy =[]
    print numLine
    #for j in range((Nsite-Lx)/2):
    for j in xrange(numLine):
      sum = [0.0, 0.0]
      for n in xrange(nNum):
        data  = ifile[n].readline().split()
        r   = float(data[0])
        tmp = float(data[1])+1J*float(data[2])
        sum[0] += tmp
        sum[1] += abs(tmp)*abs(tmp)
      ave = sum[0]/float(nNum)
      var = sum[1]/float(nNum) - abs(ave)*abs(ave)
      #var = sum[1]/float(nNum-1.0) - abs(ave)*abs(ave)*float(nNum)/float(nNum-1.0)
      sigma = sqrt(abs(var/float(nNum)))
      list = [r, abs(ave), sigma]
      phy.append(list)
      #ofile.write("{0: .10e}\t{1: .10e}\t{2: .10e}\n".format(r,abs(ave),sigma))
      #print sorted(phy)
    phy.sort(reverse=True)
    rtmp = -1.0
    sc_ave = 0.0
    sc_sigma = 0.0
    icount=0.0
    for r,ave,sigma in phy:
      if (r == rtmp):
        continue
      #if( (r > 3.0) and (r < Lx-3)):
      #if( (r > Lx/2-2) and (r < Lx/2+2)):
      if( (r > sqrt(2.0)*float(Lx)*0.5-2) and (r < sqrt(2.0)*float(Lx)*0.5)):
        icount = icount+1.0
        sc_ave = sc_ave + ave
        sc_sigma = sc_sigma + sigma*sigma
      rtmp = r
    ofile.write("{0: .10e}\t{1: .10e}\t{2: .10e}\n".format(0.0,sc_ave/icount,sqrt(sc_sigma)/icount))
    ofile.write("{0: .10e}\t{1: .10e}\t{2: .10e}\n".format(float(Lx)-1.0,sc_ave/icount,sqrt(sc_sigma)/icount))
    ofile.write("\n\n")
    rtmp = -1.0
    for r,ave,sigma in phy:
      if (r == rtmp):
        continue
      ofile.write("{0: .10e}\t{1: .10e}\t{2: .10e}\n".format(r,ave,sigma))
      rtmp = r
    #for list in phy:
    #  ofile.write("{0: .10e}\t{1: .10e}\t{2: .10e}\n".format(list[0],list[1],list[2]))
      #if((r == sqrt(2.0)*float(Lx)*0.5)):
      #  ofile2.write("{0: .10e}\t{1: .10e}\t{2: .10e}\n".format(time[k],abs(ave),sigma))
 for i in xrange(nNum):
    ifile[i].close()
 ofile.close()


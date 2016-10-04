#!/usr/bin/python

import gzip, numpy, sys

infile=sys.argv[1]
oufile=sys.argv[2]
nBlock=int(sys.argv[3])
if infile[-3:]=='.gz':
    inf=gzip.open(infile,'r')
else:
    inf=open(infile,'r')

inf.readline()
inf.readline()
inf.readline()


pressureData={}
iBlock=0
block=0
while 1:
    line=inf.readline()
    if line=='':
        break

    [timestep,ntmp]=[int(d) for d in line.split()]
    iBlock+=1
    for i in range(ntmp):
        data=inf.readline().split()
        z=float(data[1])
        density=float(data[3])
        pxx=-density*float(data[4])
        pyy=-density*float(data[5])
        pzz=-density*float(data[6])
        pxy=-density*float(data[7])
        pxz=-density*float(data[8])
        pyz=-density*float(data[9])
        if pressureData.has_key(z):
            pressureData[z][0]+=pxx
            pressureData[z][1]+=pyy
            pressureData[z][2]+=pzz
            pressureData[z][3]+=pxy
            pressureData[z][4]+=pxz
            pressureData[z][5]+=pyz
            pressureData[z][6]+=1.0
        else:
            pressureData[z]=[pxx,pyy,pzz,pxy,pxz,pyz,1.0]
            
    if iBlock==nBlock:
        block+=1

        for z in pressureData.keys():
            for i in range(6):
                pressureData[z][i]/=pressureData[z][6]


        zz=pressureData.keys()
        zz.sort()


        ouf=open(oufile+'.'+str(block),'w')
        for z in zz:
            gamma=pressureData[z][2]-0.5*(pressureData[z][0]+pressureData[z][1])
            print >> ouf, "%12.6f "*2 % (z,gamma),
            for i in range(6):
                print >> ouf, "%12.6f " % (pressureData[z][i]),
            print >> ouf

        pressureData={}
        iBlock=0

#!/usr/bin/python

import gzip, numpy, sys

## get command line arguments
inFile=sys.argv[1]
outFile=sys.argv[2]
area=float(sys.argv[3])
dz=float(sys.argv[4])

vslice=area*dz

## open input file
if inFile[-3:]=='.gz':
    inf=gzip.open(inFile,'r')
else:
    inf=open(inFile,'r')

## skip header info
inf.readline()
inf.readline()
inf.readline()
while 1:
    line=inf.readline()
    if line=='':
        break
    data=line.split()
    iStep=int(data[0])
    nLayer=int(data[1])
    print "%6d step %6d data points " % (iStep,nLayer)
    z=numpy.empty((nLayer))
    pxx=numpy.empty((nLayer))
    pyy=numpy.empty((nLayer))
    pzz=numpy.empty((nLayer))
    pxy=numpy.empty((nLayer))
    pxz=numpy.empty((nLayer))
    pyz=numpy.empty((nLayer))

    ## read in data for this timestep
    for i in range(nLayer):
        data=inf.readline().split()
        z[i]=float(data[1])
        ncount=float(data[2])
        pxx[i]=-(ncount/vslice)*float(data[4])
        pyy[i]=-(ncount/vslice)*float(data[5])
        pzz[i]=-(ncount/vslice)*float(data[6])
        pxy[i]=-(ncount/vslice)*float(data[7])
        pxz[i]=-(ncount/vslice)*float(data[8])
        pyz[i]=-(ncount/vslice)*float(data[9])
        
    gamma=pzz-0.5*(pxx+pyy)

    ouf=open(outFile+'.'+str(iStep),'w')
    for i in range(nLayer):
        print >> ouf, "%12.6f "*8 % (z[i],pxx[i],pyy[i],pzz[i],pxy[i],pxz[i],pyz[i],gamma[i])
    ouf.close()

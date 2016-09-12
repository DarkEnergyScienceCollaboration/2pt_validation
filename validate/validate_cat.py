#!/usr/bin/env python
#
# This routine validates full sky fastcat catalogs.
#
import os
import numpy as np
from optparse import OptionParser
from vc_sub import *
import sys
sys.path=["."]+sys.path
import fastcat as fc
import matplotlib.pyplot as plt
import cPickle


if os.environ.has_key('DESC_LSS_ROOT'):
    root=os.environ['DESC_LSS_ROOT']
else:
    root="/project/projectdirs/lsst/LSSWG"
#default out path
opath=root+"/LNMocks"

ljackexe="../LimberJack/LimberJack"

parser = OptionParser(usage="usage: %prog [options] mock_directory ")
parser.add_option("--opath", dest="opath", default=opath,
                  help="Path to colore output", type="string")
parser.add_option("--Nr", dest="Nr", default=10,
                  help="Number of realizations", type="int")

parser.add_option("--Nstart", dest="Nstart", default=0,
                  help="starting realization", type="int")
parser.add_option("--zstart", dest="zstart", default=0.4,
                  help="mean of first z bin", type="float")
parser.add_option("--deltaz", dest="deltaz", default=0.1,
                  help="z bin width (bins overlap 50%)", type="float")
parser.add_option("--minprobz", dest="iprobz", default=0.95,
                  help="int. probability to call galaxy in a bin", type="float")
parser.add_option("--Nz", dest="Nz", default=2,
                  help="Number of redshift bins", type="int")
parser.add_option("--zgrid_max", dest="zgrid_max", default=1.5,
                  help="max z in N(z) grid", type="float")
parser.add_option("--zgrid_dz", dest="zgrid_dz", default=0.01,
                  help="dz in N(z) grid", type="float")
parser.add_option("--Nside", dest="Nside", default=512,
                  help="Nside for healpix", type="int")
parser.add_option("--plotNz", dest="plotNz", default=False,
                  help="Plot N(z)s", action="store_true")
parser.add_option("--plotCl", dest="plotCl", default=False,
                  help="Plot C(ell)s", action="store_true")
parser.add_option("--lavg", dest="lavg", default=1.1,
                  help="ell averaging in factor", type="float")
parser.add_option("--lpow", dest="lpow", default=1.,
                  help="power in ell for plotting", type="float")
parser.add_option("--pickle", dest="picklefn", default="",
                  help="If this file exist, use it to save/load from pickle", type="string")
parser.add_option("--ljackexe", dest="ljackexe", default=ljackexe,
                  help="Limerjack executable", type="string")
parser.add_option("--getTheory", dest="getTheory", default=True,
                  help="calculate (and plot) theory")


(o, args) = parser.parse_args()

if (len(args)!=1):
    parser.error("Need to specify mock dir on commannd line")



doAna=True
if (len(o.picklefn)):
    doAna=False
    try:
        Nza,Cls,Clm,Clv,larr,zarr,prs=cPickle.load(open(o.picklefn))
        print "Loaded from pickle:",o.picklefn
        print "Reopening catalog for theory"
        file_name=o.opath+"/"+args[0]+"/catalog0.h5"
        print "Reading:",file_name
        cat=fc.Catalog()
        cat.readH5(file_name)

    except IOError:
        doAna=True
    

if doAna:
    ## DON'T DO THIS!!
    #Nza=[[]]*o.Nz
    Nza=[[] for i in range(o.Nz)]
    Cls=[[] for i in range(o.Nz*(o.Nz+1)/2)]
    prs=[[] for i in range(o.Nz*(o.Nz+1)/2)]

    for ri in range(o.Nstart,o.Nstart+o.Nr):
        file_name=o.opath+"/"+args[0]+"/catalog"+str(ri)+".h5"
        print "Reading:",file_name
        cat=fc.Catalog()
        cat.readH5(file_name)
        print "Read", len(cat.data), "galaxies."
        samps,zarr,Nzs,shotnoise=getBinnedSample(cat,o)
        for i in range(o.Nz):
            Nza[i].append(Nzs[i])
        grids=[]
        Ng=[]
        print "Making grids..."
        for s in samps:
            g=gridGals (s,o)
            N=len(s)
            grids.append(g)
            Ng.append(N)
        print "Correlating"
        k=0
        for i in range(o.Nz):
            for j in range(i,o.Nz):
                print i,j
                if (i==j):
                    Cl=healpy.anafast(grids[i])
                else:
                    Cl=healpy.anafast(grids[i],grids[j])
                Cl-=shotnoise[i,j]
                # skip monopole/quadrupole
                Cl=Cl[2:]
                Cls[k].append(Cl)
                prs[k]=(i,j)
                k+=1
                
    Nza=map(np.array,Nza)
    Cls=map(np.array,Cls)
    Clm=[Cl.mean(axis=0) for Cl in Cls]
    Clv=[Cl.var(axis=0)/o.Nr for Cl in Cls]

    print Cls[0][:,3],'X',Clm[0][3],np.sqrt(Clv[0][3])
    ## note we skipped two modes above
    larr=np.array(range(len(Clm[0])))+2

    if (len(o.picklefn)):
        print "Saving to pickle"
        cPickle.dump((Nza,Cls,Clm,Clv,larr,zarr,prs),open(o.picklefn,'w'))

def avgarr(arr,fac,divbydelta=False):
    n=[]
    llow=2.
    lmax=len(arr)+2
    while llow<lmax:
        lhigh=int(llow*fac)+1
        if divbydelta:
            n.append(arr[llow-2:lhigh-2].mean()/(lhigh-llow))           
        else:
            n.append(arr[llow-2:lhigh-2].mean())           
        llow=lhigh
    return np.array(n)

if (o.lavg>1):
    print larr[:10]
    larr=avgarr(larr,o.lavg)
    print larr[:10]
    Clm=map(lambda x:avgarr(x,o.lavg),Clm)
    Clv=map(lambda x:avgarr(x,o.lavg,divbydelta=True),Clv)


## get theory
if (o.getTheory):
    the=[]
    print prs
    for i,j in prs:
        lar,th=getLJtheory(o,cat,zarr,Nza[i].mean(axis=0), Nza[j].mean(axis=0))
        the.append(th)
    thlarr=lar
    
if (o.plotNz):
    plt.figure()
    for i in range(o.Nz):
        mz=Nza[i].mean(axis=0)
        mn=Nza[i].min(axis=0)
        mx=Nza[i].max(axis=0)
        plt.plot (zarr,mz,label="bin %i"%i)
        #plt.plot (zarr,mn,':')        
        #plt.plot (zarr,mx,':')        
    plt.xlabel("z")
    plt.ylabel("p(z)")
    plt.legend()
    plt.show()
                  

    

if (o.plotCl):
    plt.figure()
    #for Clx in Cls[0]:
    #    plt.plot(larr,Clx,'-',label="(%i,%i)"%prs[0])
    #    plt.errorbar(larr,Clm[0],yerr=np.sqrt(Clv[0]),lw=3)
        

    for k,Cl in enumerate(Clm):
        base_line=plt.errorbar(larr,Cl*larr**o.lpow,fmt='.',yerr=larr**o.lpow*np.sqrt(Clv[k]),label="(%i,%i)"%prs[k])    #    #
        if (o.getTheory):
            plt.plot(thlarr,the[k]*thlarr**o.lpow,color=base_line.lines[0].get_color())

    plt.semilogx()
    plt.legend()
    plt.show()

    
    

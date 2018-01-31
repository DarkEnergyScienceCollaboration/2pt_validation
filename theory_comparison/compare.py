#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import pyccl as ccl
import matplotlib.pyplot as plt
import sacc

nsim=100
ds=[]
ts=[]
tgs=[]
for i in np.arange(nsim) :
    print(i)
    fname_d="/global/cscratch1/sd/damonge/sims_LSST/sims_red_noshear_eh/fastcats/GaussPZ_0.02+FullSky+run%d+ztrue/twopoints%d_ns2048.sacc"%(i+1,i+1)
    fname_t="/global/cscratch1/sd/damonge/sims_LSST/sims_red_noshear_eh/fastcats/GaussPZ_0.02+FullSky+run%d+ztrue/theory_ln%d_ns2048.sacc"%(i+1,i+1)
    fname_g="/global/cscratch1/sd/damonge/sims_LSST/sims_red_noshear_eh/fastcats/GaussPZ_0.02+FullSky+run%d+ztrue/theory_ga%d_ns2048.sacc"%(i+1,i+1)
    sd=sacc.SACC.loadFromHDF(fname_d)
    st=sacc.SACC.loadFromHDF(fname_t)
    sg=sacc.SACC.loadFromHDF(fname_g)
    ds.append(sd.mean.vector)
    ts.append(st.mean.vector)
    tgs.append(sg.mean.vector)
ds=np.array(ds)
ts=np.array(ts)
tgs=np.array(tgs)
dmean=np.mean(ds,axis=0)
tmean=np.mean(ts,axis=0)
tgmean=np.mean(tgs,axis=0)
dstd=np.std(ds,axis=0)
tstd=np.std(ts,axis=0)
tgstd=np.std(tgs,axis=0)

nside=2048.
r_theta=1/(np.sqrt(3.)*nside)
def window_pixel(l) :
    x=l*r_theta
    f=0.532+0.006*(x-0.5)**2
    y=f*x*1.0
    return np.exp(-y**2/2)


cols=['r','g','b','y','c','m','k']

plt.figure()
s=sacc.SACC.loadFromHDF(fname_d)
for i,t in enumerate(s.tracers) :
    dz=t.z[1]-t.z[0]
    plt.plot(t.z,t.Nz/(dz*4*np.pi*(180*60/np.pi)**2),cols[i%7]+'-',lw=2)
plt.xlim([0,1.6])
plt.xlabel('$z$',fontsize=18)
plt.ylabel('$dN/dzd\\Omega\\,\\,[{\\rm arcmin}^{-2}]$',fontsize=18)
plt.savefig('nz.png',box_inches='tight')

plt.figure()
for t1i,t2i,ells,ndx in st.sortTracers() :
    if t1i==t2i :
        bm=window_pixel(ells)
        theo=tmean[ndx]*bm*bm
        data=dmean[ndx]
        unce=dstd[ndx]
        plt.errorbar(ells,(data-theo)/unce,
                     yerr=unce/unce,
                     fmt=cols[t1i%7]+'-')
plt.xlabel('$\\ell$',fontsize=18)
plt.ylabel('$(C^{\\rm data}_\\ell-C^{\\rm theory}_\\ell)/\\sigma(C_\\ell)$',fontsize=18)
plt.savefig('residuals.png',box_inches='tight')

plt.figure()
for t1i,t2i,ells,ndx in st.sortTracers() :
    if t1i==t2i :
        bm=window_pixel(ells)
        theo=tmean[ndx]*bm*bm
        data=dmean[ndx]
        unce=dstd[ndx]
        plt.errorbar(ells,data,yerr=unce,fmt=cols[t1i%7]+'-')
        plt.plot(ells,data,cols[t1i%7]+'-',lw=1)
        plt.plot(ells,theo,cols[t1i%7]+'--',lw=2)
plt.errorbar([-1,-1],[-1,-1],yerr=[1,1],fmt='k-',label='Data mean +- rms')
plt.plot([-1,-1],[-1,-1],'k-',lw=2,label='Log-normal theory')
plt.loglog()
plt.xlim([10,2000])
plt.ylim([1E-9,1E-3])
plt.xlabel('$\\ell$',fontsize=18)
plt.ylabel('$C^{ii}_\\ell$',fontsize=18)
plt.legend(loc='lower left',frameon=False)
plt.savefig('residuals_abs.png',box_inches='tight')


plt.figure()
for t1i,t2i,ells,ndx in st.sortTracers() :
    if (t1i==7) and (t1i==t2i) :
        bm=window_pixel(ells)
        theo=tmean[ndx]*bm*bm
        theog=tgmean[ndx]*bm*bm
        data=dmean[ndx]
        unce=dstd[ndx]
        plt.plot(ells,data ,'r-' ,lw=1,label='Data mean')
        plt.plot(ells,theo ,'r--',lw=2,label='Log-normal theory')
        plt.plot(ells,theog,'r-.',lw=2,label='Gaussian theory')
plt.legend(loc='lower left',frameon=False)
plt.loglog()
plt.xlim([10,2000])
plt.ylim([1E-9,1E-4])
plt.xlabel('$\\ell$',fontsize=18)
plt.ylabel('$C^{88}_\\ell$',fontsize=18)
plt.savefig('compare_ln.png',box_inches='tight')
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import pyccl as ccl
from scipy.integrate import quad
from scipy.interpolate import interp1d

nz=1024
zmax=2.5
nk=1024
lkmin=-5.
lkmax=2.

#Parameters for N(z) \propto z^\alpha \exp[-(z/z0)^\beta]
alpha_nz=2.
beta_nz=1.0
z0_nz=0.3
#Number density in arcmin^-2
ndens_amin2=40.

cosmo=ccl.Cosmology(Omega_c=0.266,Omega_b=0.049,h=0.69,sigma8=0.8,n_s=0.96)
karr=10.**(lkmin+(lkmax-lkmin)*np.arange(nk)/(nk-1.))
pklinarr=ccl.linear_matter_power(cosmo,1.,karr)
pknlinarr=ccl.nonlin_matter_power(cosmo,1.,karr)

zarr=zmax*np.arange(nz)/(nz-1.)
gzarr=ccl.growth_factor(cosmo,1./(1+zarr))
bzarr=0.95/gzarr
nzarr=zarr**alpha_nz*np.exp(-(zarr/z0_nz)**beta_nz)
nzf=interp1d(zarr,nzarr)
#Normalize nz
ntot=quad(nzf,0,zmax)[0]
nzarr*=ndens_amin2*60.**2/ntot
nzf=interp1d(zarr,nzarr)

print "#Gals : %lE"%(0.4*4*np.pi*(180/np.pi)**2*quad(nzf,0,zmax)[0])

np.savetxt("pk_planck.txt",np.transpose([karr,pklinarr]))
np.savetxt("pk_nlin_planck.txt",np.transpose([karr,pknlinarr]))
np.savetxt("bz_lsst.txt",np.transpose([zarr,bzarr]))
np.savetxt("nz_lsst.txt",np.transpose([zarr,nzarr]))

plt.figure();
plt.plot(karr,pklinarr,'r-',label='Linear P(k)')
plt.plot(karr,pknlinarr,'b-',label='Non-linear P(k)')
plt.legend(loc='upper right',frameon=False);
plt.loglog()
plt.figure(); plt.plot(zarr,nzarr); plt.xlabel('redshift'); plt.ylabel('#objects per sq-deg')
plt.figure(); plt.plot(zarr,bzarr); plt.xlabel('redshift'); plt.ylabel('clustering bias')
plt.show()

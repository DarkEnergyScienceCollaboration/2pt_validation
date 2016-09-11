#
# Helper routines for vc
#
import os, sys
import numpy as np
import h5py
import healpy
import pylab
import numpy.random as npr

def execCoLoRe(seed,o):
    if not o.noexec:
        try:
            os.mkdir("colore_tmp")
        except:
            pass
        writeCInis ("colore_tmp/",seed,o)
        os.system (o.cpath+"/CoLoRe colore_tmp/params.ini")
    data=np.array(h5py.File('colore_tmp/out%i_0.h5'%(seed),'r')['sources'])
    return data

def writeCInis(direct,seed,o):
    open (direct+"/params.ini",'w').write("""
prefix_out= {direct}/out{seed}
output_format= HDF5
pk_filename= {cpath}/test_files/Pk_CAMB_test.dat
nz_filename= {direct}/Nz.txt
bias_filename= {direct}/bz.txt
omega_M= 0.3
omega_L= 0.7
omega_B= 0.049
h= 0.67
w= -1.0
ns= 0.96
sigma_8= 0.8
z_min= {zmin}
z_max= {zmax}
r_smooth= 1.
n_grid= {ngrid}
    seed= {seed}""".format(direct=direct,seed=seed,zmin=o.zmin, zmax=o.zmax,
                           ngrid=o.Ngrid,cpath=o.cpath))
    with open (direct+"/Nz.txt",'w') as f:
        for z in np.linspace(o.zmin, o.zmax,100):
            f.write("%g %g\n"%(z,o.dndz))
    with open (direct+"/bz.txt",'w') as f:
        midz=(o.zmin+o.zmax)/2
        for z in np.linspace(o.zmin, o.zmax,100):
            b=o.bias+(z-midz)*o.dbiasdz
            f.write("%g %g\n"%(z,b))
            

def gridGals(gals,o):
    Ns=o.Nside
    mp=np.zeros(12*Ns*Ns)
    ra=gals['RA']
    dec=gals['DEC']
    ra*=np.pi/180.
    theta=np.pi/2-np.pi/180*dec
    iss=healpy.ang2pix(Ns,theta,ra)
    for i in iss:
        mp[i]+=1.
    return mp

def Npix(o):
    return 12*o.Nside*o.Nside
def Ngals(o):
    return (o.zmax-o.zmin)*o.dndz*(4*np.pi*(180/np.pi)**2)

def gridPoisson(o):
    Ns=o.Nside
    Np=Npix(o)
    ## total number of galaxies is
    Ng=Ngals(o)
    ng=Ng/Np
    vals=np.array([npr.poisson(ng) for i in range(Np)])
    return vals

def plotCls(Cls,err,o):
    if o.figure=='None':
        return
    els=np.array(range(len(Cls)))
    els=els+1
    pylab.plot(els,Cls,'bo')
    pylab.errorbar(els,Cls,yerr=err)
    da=np.loadtxt('/home/anze/atestcl.dat')
    tl=da[:,0]
    tt=da[:,1]
    norm=tl*(tl+1)/3./np.pi
    tt/=norm
    pylab.plot (tl,tt,'r-')
    pylab.loglog()
    pylab.xlabel('$\ell$')
    pylab.ylabel('$C_\ell \ell (\ell+1)/4\pi$')
    pylab.tight_layout()
    pylab.savefig('fig1.pdf')
    if (o.figure=='show'):
        pylab.show()
    else:
        pylab.savefig(o.figure)


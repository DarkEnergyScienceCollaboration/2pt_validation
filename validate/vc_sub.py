#
# Helper routines for vc
#
import os, sys
import os.path as path
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
    """
    Grids galaxies from RA/DEC onto healpix grid
    """
    Ns=o.Nside
    mp=np.zeros(12*Ns*Ns)
    if 'RA' in gals.dtype.names:
        ra=gals['RA']
        dec=gals['DEC']
    else:
        ra=gals['ra']
        dec=gals['dec']
        
    ra*=np.pi/180.
    theta=np.pi/2-np.pi/180*dec
    iss=healpy.ang2pix(Ns,theta,ra)
    for i in iss:
        mp[i]+=1.
    ## normalise.
    pixmean=len(gals)*1.0/Npix(o)
    mp-=pixmean
    mp/=pixmean
    return mp

def Npix(o):
    """
     Returns number of healpix pixels given Nside in options o
    """
    return 12*o.Nside*o.Nside

def Ngals(o):
    """ 
    Returns #gals in the zbin
    """
    return (o.zmax-o.zmin)*o.dndz*(4*np.pi*(180/np.pi)**2)

def getBinnedSample(cat,o):
    """ 
    Bins galaxies in a catalog given bin definitions
    in options.
    Returns tuple: 
      (list of samples, zgrid, list of Nzs)
    """
    samps=[]
    Nzs=[]
    chosl=[]
    for bi in range(o.Nz):
       zmin=o.zstart+(bi-1)*o.deltaz
       zmax=o.zstart+(bi+1)*o.deltaz
       probs=cat.photoz.iPofZ(cat.data,zmin,zmax)
       chos=np.where(probs>o.iprobz)
       chosl.append(chos)
       carr=cat.data[chos]
       print "Sample %i, z=%f..%f has %i gals."%(bi,zmin,zmax,len(carr))
       samps.append(carr)
       zgrid,Nz=cat.photoz.NofZ(carr,0.0,o.zgrid_max,o.zgrid_dz)
       Nzs.append(Nz)
    ## now make the shot-noise sample
    shotnoise=np.zeros((o.Nz,o.Nz))
    for bi in range(o.Nz):
        nbari=len(chosl[bi][0])
        for bj in range(bi,o.Nz):
            nbarj=len(chosl[bj][0])
            ## number of commons
            if (bi==bj):
                shotnoise[bi,bj]=1./nbari
            else:
                shotnoise[bi,bj]=1.*len(np.intersect1d(chosl[bi],chosl[bj]))/(nbari*nbarj)
            print bi,bj, nbari,nbarj,len(np.intersect1d(chosl[bi],chosl[bj]))
    # assuming full sphere, noise per steradian (1/nbar = 4pi/Ntot)
    shotnoise*=4*np.pi
    return samps,zgrid,Nzs,shotnoise
       
    


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
    ## this is wrong
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


def getLJtheory(o,cat,zarr,W1,W2):
    """
    Run LimberJack and derive theory given 
    options o
    catalog car
    Windows functions W1 and W2 defined over array of redshifts zarr.

    """
    a=cat.meta
    open ("ljack/params.ini",'w').write("""
#Cosmological parameters
omega_m= {Om}
omega_l= {Ol}
omega_b= {Ob}
w0= {w}
wa= 0.0
h= {h}
ns= {ns}
s8= {s8}

#Radial resolution
d_chi= 5.

#Maximum multipole
l_max= 1000

#Behavioural flags (include number counts? include lensing shear? include CMB lensing?)
do_nc= 1
has_nc_dens= 1
has_nc_rsd= 0
has_nc_lensing= 0 
do_shear= 0
do_cmblens= 0
do_isw= 0

#Angular correlation function
do_w_theta= 0
use_logbin= 0
theta_min= 0
theta_max= 10.
n_bins_theta= 15
n_bins_decade= 5

#File names (window function, bias, magnification bias, power spectrum)
window_1_fname= ljack/W1
window_2_fname= ljack/W2
bias_fname= ljack/bias.txt
sbias_fname= ljack/sz_blue.txt
pk_fname= {pkfn}

#Output prefix
prefix_out= ljack/output
""".format(Om=a["colore_omega_M"], Ol=a["colore_omega_L"], Ob=a["colore_omega_B"],
           w=a["colore_w"],h=a["colore_h"], ns=a["colore_ns"], s8=a["colore_sigma_8"],
           pkfn=a["colore_pk_filename"]))
    # write bias
    with open("ljack/bias.txt","w") as f:
        for z,bz in zip(cat.bz["z"],cat.bz["bz"]):
            f.write("%g %g \n"%(z,bz))
    for wf,fn in [(W1,"W1"),(W2,"W2")]:
        with open("ljack/"+fn,"w") as f:
            for z,w in zip(zarr,wf):
                f.write("%g %g \n"%(z,w))
    #if file exist, delete
    ddfile="ljack/output_cl_dd.txt"
    if path.isfile(ddfile):
        os.remove(ddfile)
    os.system(o.ljackexe+" ljack/params.ini")
    res=np.loadtxt(ddfile)
    zar=res[:,0]
    th=res[:,1]
    return zar,th
        


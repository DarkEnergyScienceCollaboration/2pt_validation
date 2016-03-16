#
# Helper routines for vc
#
import os, sys
import numpy as np
import h5py
import healpy
import pylab
import numpy.random as npr
from scipy.integrate import quad

def writeDists(o):
    """ Writes bias and dn/dz files for
        redmagic like sample.
    Eli Rykoff says
    
    Density is ~constant comoving density 1.5e-3 h^3 Mpc^-3, sigma_z/(1+z) 
    ~0.01, and bias is 2-2.5ish.
    """
    d=o.outpath
    fn=open (d+"/Nz.txt",'w')
    fb=open (d+"/bz.txt",'w')
    pi=np.pi
    for z in np.linspace(0.01, 1.2,1000):
            fb.write("%g %g\n"%(z,2.2))
            ## for density, need to convert Mpc/h into n/sqdeg/dz
            ## c over H(z) for Mpc/h units
            coHz=lambda z:3000./np.sqrt(0.3*(1+z)**3+0.7)
            r=quad(coHz,0.,z)[0]
            #volume of 1sq * dz
            # 4pir^2 * (1deg/rad)**2 * dr/dz
            # hMpc^3 per dz per 1sqd
            vrat=r**2 * (pi/180.)**2  * coHz(z) 
            dens=1.5e-3*vrat
            fn.write("%g %g\n"%(z,dens))

def execCoLoRe(i,o):
    dr=o.outpath+"/Set"+str(i)
    try:
        os.mkdir(dr)
    except:
        pass
    writeCInis (dr,i,o)
    if o.stype=="exec":
        exe=o.cpath+"/CoLoRe "+dr+"/params.ini"
        print exe
        os.system (exe)
    elif o.stype=="bnl":
        exe='wq sub -r "N:{cores}; threads:12; hostfile:auto; group:[new,new2]; job_name:CoLoRe" -c "OMP_NUM_THREADS=%threads% mpirun -hostfile %hostfile% {cpath}./CoLoRe {dr}/params.ini" '.format(
            cores=o.nodes*12, cpath=o.cpath, dr=dr)
        print exe
        os.system(exe)
    else:
        print "Unknown exe"

def writeCInis(direct,i,o):
    open (direct+"/params.ini",'w').write("""
prefix_out= {direct}/out
output_format= HDF5
pk_filename= {cpath}/test_files/Pk_CAMB_test.dat
nz_filename= {opath}/Nz.txt
bias_filename= {opath}/bz.txt
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
    seed= {seed}""".format(direct=direct,seed=o.seed+i,zmin=o.zmin, zmax=o.zmax,
                           ngrid=o.Ngrid,cpath=o.cpath,opath=o.outpath))

            

#
# Helper routines for vc
#
import os, sys
import numpy as np
try:
    from astropy.cosmology import Planck15 as co
except:
    from astropy.cosmology import Planck13 as co
    print "Warning: using Planck13 cosmology rather than Planck15!"
import astropy.constants as const
import astropy.units as u
from scipy.integrate import quad

def writeDists(o):
    """ Writes bias and dn/dz files for
        redmagic like sample.
    Eli Rykoff says

    Density is ~constant comoving density 1.5e-3 h^3 Mpc^-3, sigma_z/(1+z)
    ~0.01, and bias is 2-2.5ish.
    """

    zmin=0.01
    zmax=1.2
    Nz=1000
    rho_comoving=1.5e-3
    #Nz shaping
    zshape=1.0

    d=o.outpath
    fn=open (d+"/Nz.txt",'w')
    fb=open (d+"/bz.txt",'w')
    pi=np.pi
    for z in np.linspace(zmin, zmax,Nz):
            fb.write("%g %g\n"%(z,2.2))
            ## for density, need to convert Mpc/h into n/sqdeg/dz
            ## c over H(z) for Mpc/h units
            coHz=(const.c/co.H(z)).to(u.Mpc).value*co.h
            # radius in Mpc/h
            r=co.comoving_distance(z).to("Mpc").value*co.h
            #volume of 1sq * dz
            # 4pir^2 * (1deg/rad)**2 * dr/dz
            # hMpc^3 per dz per 1sqd
            vrat=r**2 * (pi/180.)**2  * coHz
            dens=rho_comoving*vrat
            ## shape distribution to avoid sharep cut
            if (z>zshape):
                sup=(z-zshape)/(zmax-zshape)
                dens*=np.exp(-10*sup**2)
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
    elif o.stype=="cori":
        open(dr+"/script.sm","w").write("""#!/bin/bash -l
#SBATCH --partition regular
#SBATCH --nodes {nodes}
#SBATCH --time={time}
#SBATCH --job-name=CoLoRe_{i}
#SBATCH -C haswell
#SBATCH --account=m1727
cd {dr}
srun -n {cores} {cpath}/CoLoRe ./params.cfg >slurm.log 2>slurm.err
""".format(nodes=o.nodes, cores=o.nodes*64, cpath=o.cpath, dr=dr,i=i,time=o.time))
        os.system("sbatch "+dr+"/script.sm")
    else:
        print "Unknown exe"

def writeCInis(direct,i,o):
    open (direct+"/params.cfg",'w').write("""
global:{
  prefix_out= {direct}/out
  output_format= HDF5
  output_density= false
  pk_filename= {cpath}/test_files/Pk_CAMB_test.dat
  z_min= {zmin}
  z_max= {zmax}
  seed = {seed}
  write_pred = true
  pred_dz = 0.1
}
field_par:{
  r_smooth= 1.
  n_grid= {ngrid}
  smooth_potential= true
  dens_type = {dens_type}
  lpt_buffer_fraction = 0.5
  lpt_interp_type = 1
  output_lpt = 0
}
cosmo_par:{
  omega_M= 0.3
  omega_L= 0.7
  omega_B= 0.049
  h= 0.67
  w= -1.0
  ns= 0.96
  sigma_8= 0.8
  }
srcs1:{
  nz_filename= {opath}/Nz.txt
  bias_filename= {opath}/bz.txt
  include_shear = true
}
kappa:{
  z_out = {z_out}
  nside = {nside}
}

""".format(direct=direct,seed=o.seed+i,zmin=o.zmin, zmax=o.zmax,
                           ngrid=o.Ngrid,cpath=o.cpath,opath=o.outpath,
                           dens_type=o.dens_type,z_out=o.z_out,nside=o.nside))

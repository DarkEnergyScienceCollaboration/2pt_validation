for i in {17..88}
do
cat >> param_files/param_colore_${i}.cfg << EOF
global:
{
  #Output prefix. Output will be in prefix_<node ID>.<fits/txt>
  prefix_out= "/PATH/TO/OUTPUT/Mock_${i}_$((1000+i))_3072";
  #Output format. Select HDF5, FITS or ASCII
  output_format= "HDF5";
  #Output Gaussian overdensity field at z=0?
  output_density= false
  #Path to power spectrum at z=0. Power spectrum file must
  #be in CAMB format: k (h/Mpc), P(k) (Mpc/h)^3.
  pk_filename= "/PATH/TO/PK/pk_planck.txt"
  #This redshift range also defines the size of the box
  z_min= 0.001
  z_max= 2.5
  #RNG seed note that output will depend on number of nodes, etc not only
  #on the RNG seed
  seed=$((1000+i))
  write_pred= true
  pred_dz=0.1
}
field_par:
{
  #Extra Gaussian smoothing scale [Mpc/h] (set to a
  #negative value if you don't want any smoothing)
  r_smooth= 0.
  #Do you want to smooth the Newtonian potential as well?
  smooth_potential= true
  #Will use a Cartesian grid with n_grid^3 cells
  n_grid= 3072
  #Density field type
  # 0-lognormal
  # 1-1LPT
  # 2-1LPT
  dens_type= 2
  #If dens_type==1 or 2, buffer size (fraction per particle)
  lpt_buffer_fraction= 0.5
  #If dens_type==1 or 2, scheme to interpolate particle
  #positions into a grid
  # 0-NGP
  # 1-CIC
  # 2-TSC
  lpt_interp_type= 1
  #Set to 1 if you want to output the LPT particle positions
  output_lpt= 0
}
cosmo_par:
{
  #Non-relativistic matter
  omega_M= 0.315
  #Dark energy
  omega_L= 0.685
  #Baryons
  omega_B= 0.049
  #Hubble parameter (in units of 100 km/s/Mpc)
  h= 0.69
  #Dark energy equation of state
  w= -1.0
  #Primordial scalar spectral index, used only to extrapolate
  #P(k) at low k end (-3 used at high k end)
  ns= 0.96
  #Power spectrum normalization. The input power spectrum will be
  #renormalized to this sigma8
  sigma_8= 0.8
}
#For each galaxy population, create a section called srcsX, starting with X=1
srcs1:
{
  #Path to N(z) file. Should contain two columns
  # 1-> z, 2-> dN(z)/dz*dOmega
  # with dN/dzdOmega in units of deg^-2
  # Include one name per population, separated by spaces
  nz_filename= "/PATH/TO/NZ/nz_lsst.txt"
  #Path to bias file. Should contain two columns
  # 1-> z, 2-> b(z)
  # Include one name per population, separated by spaces
  bias_filename= "/PATH/TO/BZ/bz_lsst.txt"
  #Do you want to include shear ellipticities?
  include_shear= true
}

#Kappa fields
kappa:
{
  z_out= [0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4]
  nside= 1024
}
EOF


cat >> run_colore_files/run_colore_${i}.sh << EOF
#!/bin/bash -l
#SBATCH --partition regular
####SBATCH --qos premium
#SBATCH --nodes 32
#SBATCH --time=00:15:00
#SBATCH --job-name=CoLoRe_WL_${i}
#SBATCH -C haswell
#module load gsl
#module load fftw
export OMP_NUM_THREADS=64
export LD_LIBRARY_PATH=/opt/cray/hdf5/1.8.16/INTEL/15.0/lib:${LD_LIBRARY_PATH}          # for bash
export LD_LIBRARY_PATH=/opt/cray/hdf5-parallel/1.8.16/INTEL/15.0/lib:${LD_LIBRARY_PATH}
srun -n 32 -c 64 /PATH_TO_COLORE/CoLoRe param_files/param_colore_${i}.cfg 
EOF

sbatch run_colore_files/run_colore_${i}.sh
done


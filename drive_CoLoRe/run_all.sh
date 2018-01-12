#!/bin/bash

ng=4096
rsm=2.0
nnod=8
timelim=10
predir="/global/cscratch1/sd/damonge/sims_LSST"
rundir=${predir}"/sims_red_noshear/"
exec_colore="/global/homes/d/damonge/Codes/CoLoRe/CoLoRe"
which_partition="regular"
mkdir -p ${rundir}
mkdir -p ${rundir}/param_files
mkdir -p ${rundir}/run_colore_files

for i in {1..100}
do
    parfile=${rundir}/param_files/param_colore_${i}.cfg
    cat > ${parfile} << EOF
global:
{
  #Output prefix. Output will be in prefix_<node ID>.<fits/txt>
  prefix_out= "${rundir}/Mock_${i}_$((1000+i))_${ng}";
  #Output format. Select HDF5, FITS or ASCII
  output_format= "HDF5";
  #Output Gaussian overdensity field at z=0?
  output_density= false
  #Path to power spectrum at z=0. Power spectrum file must
  #be in CAMB format: k (h/Mpc), P(k) (Mpc/h)^3.
  pk_filename= "${predir}/sample_run/Pk_CAMB_test.dat"
  #This redshift range also defines the size of the box
  z_min= 0.001
  z_max= 1.6
  #RNG seed note that output will depend on number of nodes, etc not only
  #on the RNG seed
  seed=$((1000+i))
  write_pred= false
  just_write_pred= false
  pred_dz=0.1
}

field_par:
{
  #Extra Gaussian smoothing scale [Mpc/h] (set to a
  #negative value if you don't want any smoothing)
  r_smooth= ${rsm}
  #Do you want to smooth the Newtonian potential as well?
  smooth_potential= true
  #Will use a Cartesian grid with n_grid^3 cells
  n_grid= ${ng}
  #Density field type
  # 0-lognormal
  # 1-1LPT
  # 2-1LPT
  dens_type= 0
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
  omega_M= 0.3
  #Dark energy
  omega_L= 0.7
  #Baryons
  omega_B= 0.05
  #Hubble parameter (in units of 100 km/s/Mpc)
  h= 0.7
  #Dark energy equation of state
  w= -1.0
  #Primordial scalar spectral index, used only to extrapolate
  #P(k) at low k end (-3 used at high k end)
  ns= 0.96
  #Power spectrum normalization. The input power spectrum will be
  #renormalized to this sigma8
  sigma_8= 0.803869
}

#For each galaxy population, create a section called srcsX, starting with X=1
srcs1:
{
  #Path to N(z) file. Should contain two columns
  # 1-> z, 2-> dN(z)/dz*dOmega
  # with dN/dzdOmega in units of deg^-2
  # Include one name per population, separated by spaces
  nz_filename= "${predir}/sample_LSST/NzRed.txt"
  #Path to bias file. Should contain two columns
  # 1-> z, 2-> b(z)
  # Include one name per population, separated by spaces
  bias_filename= "${predir}/sample_LSST/BzBlue.txt"
  #Do you want to include shear ellipticities?
  include_shear= false
  #Do you want to store line-of-sight skewers for each object?
  store_skewers= false
}

#Kappa fields
kappa0:
{
  z_out= [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6]
  nside= 1024
}
EOF

    runfile=${rundir}/run_colore_files/run_colore_${i}.sh
    cat > ${runfile} << EOF
#!/bin/bash -l
#SBATCH --partition ${which_partition}
####SBATCH --qos premium
#SBATCH --nodes ${nnod}
#SBATCH --time=00:${timelim}:00
#SBATCH --job-name=CoLoRe_test_${i}
#SBATCH -C haswell
#module load gsl
#module load fftw
export OMP_NUM_THREADS=64
#export LD_LIBRARY_PATH=/opt/cray/hdf5/1.8.16/INTEL/15.0/lib:${LD_LIBRARY_PATH}          # for bash
#export LD_LIBRARY_PATH=/opt/cray/hdf5-parallel/1.8.16/INTEL/15.0/lib:${LD_LIBRARY_PATH}
srun -n ${nnod} -c 64 ${exec_colore} ${parfile}
EOF

#${exec_colore}_nompi --test-memory ${parfile}
sbatch ${runfile}
cat ${runfile}
done

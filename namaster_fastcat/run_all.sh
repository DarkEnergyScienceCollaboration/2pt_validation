#!/bin/bash -l
#SBATCH --partition debug
##SBATCH --qos premium
#SBATCH --nodes 8
#SBATCH --time=00:30:00
#SBATCH --job-name=CoLoRe_WL
#SBATCH --account=m1727
#SBATCH -C haswell

srun -n 32 python namaster_interface.py --input-file /global/cscratch1/sd/jsanch87/CoLoRe_LN100/fastcat_outputs/170522+GaussPZ_0.05+hpix_nodither_0.1+catalog27/fastcat_catalog0.h5 --output-file catalog27.sacc --nz-bins-file bins_z.txt --templates none --delta-ell 25 --nside 2048 --mpi 



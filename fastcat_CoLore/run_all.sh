#!/bin/bash -l
#SBATCH --partition regular
##SBATCH --qos premium
#SBATCH --nodes 8
#SBATCH --time=00:08:00
#SBATCH --job-name=CoLoRe_WL
#SBATCH --account=m1727
#SBATCH -C haswell
module load python/2.7-anaconda
module load h5py-parallel
srun -n 32 python mkcat.py --params_file=../mass_mapping/param_test.cfg --opath=/global/cscratch1/sd/jsanch87/CoLoRe_LN100/fastcat_outputs --pz_type=gauss --pz_sigma=0.05 --mpi --oextra='20percent' 



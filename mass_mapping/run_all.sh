#!/bin/bash -l
#SBATCH --partition debug
####SBATCH --qos premium
#SBATCH --nodes 64
#SBATCH --time=00:30:00
#SBATCH --job-name=CoLoRe_WL
#SBATCH --account=m1727
#SBATCH -C haswell
#module load gsl
#module load fftw
export OMP_NUM_THREADS=64
srun -n 64 -c 64 /global/homes/j/jsanch87/CoLoRe/CoLoRe param_hsc.cfg 


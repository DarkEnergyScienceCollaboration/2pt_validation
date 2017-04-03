#!/bin/bash -l
#SBATCH --partition debug
#SBATCH --nodes 32
#SBATCH --time=00:29:59
#SBATCH --job-name=CoLoRe_WL
#SBATCH --account=m1727
#SBATCH -C haswell
#module load gsl
#module load fftw
srun -n 256 ./CoLoRe param_colore.cfg 


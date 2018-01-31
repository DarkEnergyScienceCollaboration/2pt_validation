#!/bin/bash -l

#Time limit in minutes
timelim=15
#Number of nodes
nnod=1
#Queue to use
which_partition="regular"

#Output paths
predir="/global/cscratch1/sd/damonge/sims_LSST"
#rundir=${predir}"/sims_red_noshear/"
rundir=${predir}"/sims_red_noshear_eh/"

#Gaussian photo_z widht /(1+z)
pz=0.02

#Launch all jobs
for i in {1..100}
do
    parfile=${rundir}/param_files/param_colore_${i}.cfg
    opath=${rundir}/fastcats
    oextra=run${i}

    runfile=${rundir}/run_colore_files/run_fastcat_${i}.sh
    runcmd="python mkcat.py --params_file=${parfile} --opath=${opath} --pz_type=gauss --pz_sigma=${pz} --oextra=${oextra} --wf_type=none --ztrue"

    cat > ${runfile} <<EOF
#!/bin/bash -l
#SBATCH --partition ${which_partition}
##SBATCH --qos premium
#SBATCH --nodes ${nnod}
#SBATCH --time=00:${timelim}:00
#SBATCH --job-name=fastcat_${i}
#SBATCH --account=m1727
#SBATCH -C haswell
module load python/2.7-anaconda
#module load h5py-parallel
srun -n ${nnod} ${runcmd}

EOF

    cat ${runfile}
    sbatch ${runfile}
#    echo ${runcmd}
done

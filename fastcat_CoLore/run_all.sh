#!/bin/bash -l

timelim=15
predir="/global/cscratch1/sd/damonge/sims_LSST"
rundir=${predir}"/sims_red_noshear/"
pz=0.02
nnod=1
which_partition="regular"

for i in {51..100}
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

#python mkcat.py --params_file=/global/cscratch1/sd/damonge/sims_LSST/sims_red_noshear/param_files/param_colore_1.cfg --opath=/global/cscratch1/sd/damonge/sims_LSST/sims_red_noshear/fastcats --pz_type=gauss --pz_sigma=0.02 --wf_type=none --oextra='testing'
#srun -n ${nnod} python mkcat.py --params_file=${parfile} --opath=${opath} --pz_type=gauss --pz_sigma=${pz} --mpi --oextra='20percent' 

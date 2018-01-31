#!/bin/bash -l

#Number of nodes
nnod=1
#Time limit in minutes
timelim=5
#Queue to use
which_partition="regular"

#Output paths
predir="/global/cscratch1/sd/damonge/sims_LSST"
#rundir=${predir}"/sims_red_noshear/"
rundir=${predir}"/sims_red_noshear_eh/"

#Map resolutions
nside=2048
#Photo-z width
pz=0.02
#Window used to generate the catalogs
wintype=FullSky

#Type of transfer function to use when computing the theory
#transfer_function="boltzmann"
transfer_function="eisenstein_hu"

#Launch all jobs
for i in {1..100}
do
    infile=${rundir}fastcats/GaussPZ_${pz}+${wintype}+run${i}+ztrue/twopoints${i}_ns${nside}.sacc
    #outfile=${rundir}fastcats/GaussPZ_${pz}+${wintype}+run${i}+ztrue/theory_ln${i}_ns${nside}.sacc
    outfile=${rundir}fastcats/GaussPZ_${pz}+${wintype}+run${i}+ztrue/theory_ga${i}_ns${nside}.sacc
    parfile=${rundir}param_files/param_colore_${i}.cfg

    runfile=${rundir}/run_colore_files/run_theory_${i}.sh
    #runcmd="python mk_theory.py -i ${infile} -o ${outfile} --param-file ${parfile} --include-rsd false --do-lognormal true --transfer-function ${transfer_function}"
    runcmd="python mk_theory.py -i ${infile} -o ${outfile} --param-file ${parfile} --include-rsd false --transfer-function ${transfer_function}"
    cat > ${runfile} <<EOF
#!/bin/bash -l
#SBATCH --partition ${which_partition}
##SBATCH --qos premium
#SBATCH --nodes ${nnod}
#SBATCH --time=00:${timelim}:00
#SBATCH --job-name=mktheo_${i}
#SBATCH --account=m1727
#SBATCH -C haswell

srun -n ${nnod} ${runcmd}

EOF
    cat ${runfile}
    sbatch ${runfile}
    #echo $runcmd
    #${runcmd}
done

#!/bin/bash -l

timelim=5
predir="/global/cscratch1/sd/damonge/sims_LSST"
rundir=${predir}"/sims_red_noshear/"
nnod=1
nside=2048
pz=0.02
wintype=FullSky
which_partition="regular"

for i in {1..50}
do
    infile=${rundir}fastcats/GaussPZ_${pz}+${wintype}+run${i}+ztrue/twopoints${i}_ns${nside}.sacc
    outfile=${rundir}fastcats/GaussPZ_${pz}+${wintype}+run${i}+ztrue/theory${i}_ns${nside}.sacc
    parfile=${rundir}param_files/param_colore_${i}.cfg

    runfile=${rundir}/run_colore_files/run_theory_${i}.sh
    runcmd="python mk_theory.py -i ${infile} -o ${outfile} --param-file ${parfile} --include-rsd false"
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
    #${runcmd}
done

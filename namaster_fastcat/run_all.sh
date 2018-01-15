#!/bin/bash -l

timelim=20
predir="/global/cscratch1/sd/damonge/sims_LSST"
rundir=${predir}"/sims_red_noshear/"
nnod=1
nside=2048
pz=0.02
delta_ell=25
wintype=FullSky
which_partition="regular"
binsfile=${predir}/sample_LSST/bins_z_red.txt

for i in {51..100}
do
    infile=${rundir}fastcats/GaussPZ_${pz}+${wintype}+run${i}+ztrue/fastcat_catalog0.h5
    outfile=${rundir}fastcats/GaussPZ_${pz}+${wintype}+run${i}+ztrue/twopoints${i}_ns${nside}.sacc

    runfile=${rundir}/run_colore_files/run_namaster_${i}.sh
    runcmd="python namaster_interface.py --input-file ${infile} --output-file ${outfile} --nz-bins-file ${binsfile} --templates none --delta-ell ${delta_ell} --nside ${nside}"
    cat > ${runfile} <<EOF
#!/bin/bash -l
#SBATCH --partition ${which_partition}
##SBATCH --qos premium
#SBATCH --nodes ${nnod}
#SBATCH --time=00:${timelim}:00
#SBATCH --job-name=namaster_${i}
#SBATCH --account=m1727
#SBATCH -C haswell

srun -n ${nnod} ${runcmd}

EOF
    cat ${runfile}
    sbatch ${runfile}
#    echo ${runcmd}
#    ${runcmd}
done

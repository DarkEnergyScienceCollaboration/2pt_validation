#!/bin/bash -l

timelim=20
predir="/global/cscratch1/sd/damonge/sims_LSST"
rundir=${predir}"/sims_red_noshear/"
nnod=1
nside=2048
which_partition="regular"

for i in {46..85}
do
    infile=${rundir}fastcats/180106+GaussPZ_0.02+FullSky+run${i}+ztrue/fastcat_catalog0.h5
    outfile=${rundir}fastcats/180106+GaussPZ_0.02+FullSky+run${i}+ztrue/twopoints${i}_ns${nside}.sacc
    binsfile=${predir}/sample_LSST/bins_z_red.txt

    runfile=${rundir}/run_colore_files/run_namaster_${i}.sh
    cmnd="python namaster_interface.py --input-file ${infile} --output-file ${outfile} --nz-bins-file ${binsfile} --templates none --delta-ell 25 --nside ${nside}"
    cat > ${runfile} <<EOF
#!/bin/bash -l
#SBATCH --partition ${which_partition}
##SBATCH --qos premium
#SBATCH --nodes ${nnod}
#SBATCH --time=00:${timelim}:00
#SBATCH --job-name=namaster_${i}
#SBATCH --account=m1727
#SBATCH -C haswell

srun -n ${nnod} ${cmnd}

EOF
    cat ${runfile}
    sbatch ${runfile}
#    ${cmnd}
done

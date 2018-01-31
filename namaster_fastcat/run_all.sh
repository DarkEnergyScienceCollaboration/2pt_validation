#!/bin/bash -l

#Number of nodes
nnod=1
#Time limit in minutes
timelim=20
#Queue to use
which_partition="regular"

#Output paths
predir="/global/cscratch1/sd/damonge/sims_LSST"
#rundir=${predir}"/sims_red_noshear/"
rundir=${predir}"/sims_red_noshear_eh/"
#Gaussian photo-z used when generating the catalog
pz=0.02
#Window function used when generating the catalog
wintype=FullSky

#Angular resolution of the maps used to compute power spectra
nside=2048
#Bandpower width for the output power spectra
delta_ell=25
#File describing the redshift bins
binsfile=${predir}/sample_LSST/bins_z_red.txt

#Launch all jobs
for i in {1..100}
do
    infile=${rundir}fastcats/GaussPZ_${pz}+${wintype}+run${i}+ztrue/fastcat_catalog0.h5
    outfile=${rundir}fastcats/GaussPZ_${pz}+${wintype}+run${i}+ztrue/twopoints${i}_ns${nside}.sacc
    runfile=${rundir}/run_colore_files/run_namaster_${i}.sh
    runcmd="python namaster_interface.py --input-file ${infile} --output-file ${outfile} --nz-bins-file ${binsfile} --templates none --delta-ell ${delta_ell} --nside ${nside} --subtract-shot-noise"
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

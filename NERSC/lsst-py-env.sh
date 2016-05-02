# you run this by saying source NERSC/lsst-py-env.sh
PS_SAVE=$PS1
source /project/projectdirs/lsst/lsstDM/Cori/conda-env/2016-04-15/setupVanillaStack.sh
export PS1=$PS_SAVE
source eups-setups.sh
setup sims_maf     

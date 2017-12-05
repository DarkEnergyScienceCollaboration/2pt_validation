## Shifter @ NERSC

Shifter at NERSC allows you to run docker images instead of
executables, including MPI process.

You can inquire, which images are available like this:

```
anzes@cori06:~/shifter > shifterimg images | grep desc_lss
cori       docker     READY    5f34ad81e1   2017-12-05T05:41:11 slosar/desc_lss:v0            
```

To run, you need to formulate the slurm script as follows:

```
#!/bin/bash
#SBATCH --image=docker:slosar/desc_lss:v0
#SBATCH --nodes=1
#SBATCH --partition=debug
#SBATCH -C haswell
#SBATCH --volume="/global/homes/a/anzes/shifter/output:/output;/global/homes/a/anzes/shifter/input:/input"
srun -n 4 shifter /home/lss/CoLoRe/runCoLoRe /input/param.cfg
```

If you stare into above script for sufficiently long, everything
should become obvious. The -n 4 runs 4 MPI nodes, so there are 64
threads available for each MPI node (all one a single physical node).



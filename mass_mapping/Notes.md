# Notes


## Memory consumption for CoLoRe
For an Nside=2048 simulation, in GB:
* Gaussian fields: 64
* 1LPT: 144
* 2LPT: 385
* Spherical shells: 186

The instantaneous memory that must be available is Mem_total = Gaussian + MAX(1LPT,2LPT,Spherical shells)

The total memory scales with Nside^3


## Run time tests
Javi has run a 2LPT Nside=4096 sim in ~30min using 32 nodes and 8 cores per node on Cori.
Also tried to run 2LPT Nside=6144 on 256 KNL nodes in Cori (now only available in the debug queue).
Each node has 272 threads but it ran out of time (30 minutes).

## Resolution
The following spatial resolutions (in Mpc/h) are achievable with different combinations of Nside and zmax:

| Nside \ zmax |  1.0  |  1.5  |  2.0  |  2.5  |
| ------------ |:-----:|:-----:|:-----:|:-----:|
|  2048        | 2.2   |  2.9  |  3.4  |  4.0  |
|  4096        | 1.1   |  1.5  |  1.7  |  2.0  |
|  6144        | 0.8   |  1.0  |  1.2  |  1.3  |
|  8192        | 0.6   |  0.7  |  0.9  |  1.0  |

# Notes

## Memory consumption for CoLoRe
For an Nside=2048 simulation, in GB:
* Gaussian fields: 64
* 1LPT: 144
* 2LPT: 385
* Spherical shells: 186

The instantaneous memory that must be available is Mem_total = Gaussian + MAX(1LPT,2LPT,Spherical shells)
The total memory scales with Nside^3

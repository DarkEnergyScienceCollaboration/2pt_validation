# Mass map simulation notes

## Introduction

In order to test the different mass mapping algorithms for the LSST-DESC effort
to create mass maps for HSC data, we generated one full-sky simulation which
tries to resemble the HSC-Wide survey characteristics.

This simulation has been generated using `CoLoRe`: http://github.com/damonge/CoLoRe.
The `CoLoRe` parameter file used to generate it is in this directory `param_hsc.cfg`.

The simulation is a 2LPT dark matter simulation with a linear bias poisson sampling
to generate the source catalog. It is based on a Planck LCDM cosmology. The power 
spectrum used to generated this simulation can be found at `pk_planck.txt`.
The sources are generated following the N(z) contained in `nz_hsc.txt` and with the
bias at `bz_lsst.txt`.

We used a box with 4096^3 cells from redshift $z_{min}=0.05$ to $z_{max}=2.5$. This
leads to a resolution of 1.96 Mpc/h. We generated 12 lens planes at $z=0.2, 0.4, 0.6,
0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4$.  

## N(z)



## Catalogs

The catalogs are `FITS` tables containing the following columns:

* `RA`: Right ascension in degrees of the source.
* `DEC`: Declination in degrees of the source.
* `Z_COSMO`: Cosmological redshift of the source.
* `DZ_RSD`: Redshift change due to redshift space distortions.
* `E1`: Noiseless $e_{1}$ ellipticity of the source.
* `E2`: Noiseless $e_{2}$ ellipticity of the source.

There are 64 files, each one corresponds roughly to 2 to 3 healpix pixels of Nside=4.
The correspondence between each file and the pixel number is:
pixels = i, i+64, i+2*64, ... while pixels<192

Where i is the index of the file (from 0 to 64).

## CPU requirements

The simulation was run using 64 Cori Haswell nodes, with 64 threads each (32 cores/node 2
threads). The total running time was 810.3 seconds.

## License

*TBD*      

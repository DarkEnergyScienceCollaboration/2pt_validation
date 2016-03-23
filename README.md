# 2pt_validation

This repo was created to track progress on 2PT validation (LSS 1.1.X tasks).

There are three pieces that make things up:

https://github.com/damonge/CoLoRe -- highly optimized C code that makes
lognormal mocks.

https://github.com/slosar/fastcat -- public repo that contains a Catalog
object for abstract storage of mock catalogs (used to be able to be use
randomfield as backend, but this has (temporarily?) been pushed into
attick)

https://github.com/DarkEnergyScienceCollaboration/2pt_validation -- this repo that contains various glues between CoLoRe and fastcat, etc.

To see how this works in action, look at the following ipython
notebook,

https://github.com/DarkEnergyScienceCollaboration/2pt_validation/blob/master/ipynb/treecorr_example.ipynb

inspired by Javier's hack last week. It loads a fastcat catalog and then
uses treecorr to calculate correlation function. As you know, I'm not a
fan of correlation functions, but it works great as an example. Note that
fastcat objects are not braindead, for example, when loading the HDF5
file it will work out that window function is a simple declination cut
and will set up and appropriate window object -- this can then be used
to generate random catalogs. However, at the same time HDF5 format was
written in a way that doesn't require you to go through the python route
so you can use any language that you want. The format is described in
this dir:

https://github.com/slosar/fastcat/tree/master/hdf5_format_doc

The 0.1 format is what was used for the 10 mocks last week and 0.2 is a
minor update since people expressed wish that the mocks should really be
unique, all the data mucking is done on the mock-creation level (i.e. I
already add gaussian PZ errors for you etc. 

There are three basic directories:
* ipyn - ipython notebooks
* drive_CoLoRe : helps submit CoLoRe jobs
* validate_CoLoRe : a very simply script for CoLoRe validation, a 0th order test so that we don't issue people junk
* fastcat_CoLoRe : produce fastcat catalogs from CoLoRe outputs. Cuts the data, adds PZ errors, etc.

# NERSC resources

Files relevant to this project are stored under 
```
/project/projectdirs/lsst/LSSWG
```
on NERSC systems (e.g. cori)

The following subdirectories are available:
* LNMocks -- actual mocks, created in subdirectories tagged by date (and maybe something else)
* HumnaWindowFunctions -- two maps provided by Humna, to be used with issue #1




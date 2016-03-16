# 2pt_validation

This repo was created to track progress on 2PT validation (LSS 1.1.X tasks).

It relies heavily on fastcat, which is on public github. The HDF5 format used in the mock catalog
is described there.

The best way to see what is going on is to inspect the notebook that
calculates correlation function using treecorr (inspired by
Francisco's Javier)

There are three basic directory:
* ipyn - ipython notebooks
* drive_CoLoRe : helps submit CoLoRe jobs
* validate_CoLoRe : a very simply script for CoLoRe validation, a 0th order test so that we don't issue people junk
* fastcat_CoLoRe : produce fastcat catalogs from CoLoRe outputs. Cuts the data, adds PZ errors, etc.


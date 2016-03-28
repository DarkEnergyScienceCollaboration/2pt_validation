#!/usr/bin/env python
import os
from optparse import OptionParser
from vc_sub import *
import numpy as np 

parser = OptionParser()
if (os.environ.has_key("COLORE_EXEC")):
    cpath=os.environ["COLORE_EXEC"]
else:
    cpath="../CoLoRe/"

opath="/project/projectdirs/lsst/LSSWG/colore_raw"

parser.add_option("--cpath", dest="cpath", default=cpath,
                  help="Path to CoLoRe (will add /CoLoRe for executable)", type="string")
parser.add_option("--outpath", dest="outpath", default=opath,
                  help="Path to output path", type="string")
parser.add_option("--stype", dest="stype", default="nersc",
                  help="Submission type (exec,wq,nersc)", type="string")
parser.add_option("--N", dest="Nr", default=10,
                  help="Number of realizations", type="int")
parser.add_option("--Nstart", dest="Nstart", default=0,
                  help="starting realization", type="int")
parser.add_option("--Ngrid", dest="Ngrid", default=128,
                  help="FFT size", type="int")
parser.add_option("--nodes", dest="nodes", default=12,
                  help="Number of nodes to be used", type="int")
parser.add_option("--time", dest="time", default="00:15:00",
                  help="time to put into slurm field", type="string")
parser.add_option("--seed", dest="seed", default=1000,
                  help="Random seed", type="int")
parser.add_option("--zmin", dest="zmin", default=0.05,
                  help="zmin", type="float")
parser.add_option("--zmax", dest="zmax", default=1.1,
                  help="zmax", type="float")

(o, args) = parser.parse_args()

#write distributions 
writeDists(o) 
# now loop over realizations
for i in range(o.Nstart,o.Nr):
    execCoLoRe(i,o)

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

parser.add_option("--cpath", dest="cpath", default=cpath,
                  help="Path to CoLoRe (will add /CoLoRe for executable)", type="string")
parser.add_option("--N", dest="Nr", default=10,
                  help="Number of realizations", type="int")
parser.add_option("--zmin", dest="zmin", default=0.5,
                  help="minz", type="float")
parser.add_option("--zmax", dest="zmax", default=0.6,
                  help="maxz", type="float")
parser.add_option("--bias", dest="bias", default=2.0,
                  help="bias at (minz+maxz)/2", type="float")
parser.add_option("--dbiasdz", dest="dbiasdz", default=0,
                  help="dbias/dz", type="float")
parser.add_option("--dndz", dest="dndz", default=200.0,
                  help="dn/dz in #/sqdeg/dz", type="float")
parser.add_option("--Ngrid", dest="Ngrid", default=128,
                  help="FFT grid size", type="int")
parser.add_option("--Nside", dest="Nside", default=256,
                  help="Healpix Nside", type="int")
parser.add_option("--poisson", dest="poisson", default=False,
                  help="Sample Poisson, rather than galaxies", action="store_true")
parser.add_option("--noexec", dest="noexec", default=False,
                  help="Don't execute color, use previous output", action="store_true")
parser.add_option("--figure", dest="figure", default='show',
                  help="Save figure to, None for nothing, show to show", type="string")

(o, args) = parser.parse_args()

Ngt=0
sCls,ssCls=None,None
Ng=Ngals(o)
Np=Npix(o)
ng=Ng/Np

for i in range(o.Nr):
    # execute CoLoRe
    if o.poisson:
        hpgrid=gridPoisson(o)
        Ngt+=hpgrid.sum()
    else:
        gals=execCoLoRe(i,o)
        hpgrid=gridGals (gals,o)
        # total number of galaxies
        Ngt+=len(gals)
    # grid onto healpix grid
    hpgrid=hpgrid/ng
    Cls=healpy.anafast(hpgrid)
    if (sCls==None):
        sCls=Cls
        ssCls=Cls*Cls
    else:
        sCls+=Cls
        ssCls+=Cls*Cls

sCls/=o.Nr
ssCls/=o.Nr
Ngt/=o.Nr
ssCls-=sCls*sCls
err=np.sqrt(ssCls)
sCls-=1/(Ng/(4*np.pi))
print "Avg Ngals=",Ngt, Ng, ng
plotCls(sCls,err,o)

#print sCls.sum(), sCls.mean()

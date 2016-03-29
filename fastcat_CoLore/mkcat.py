#!/usr/bin/env python
import sys
from subs import *
sys.path=["."]+sys.path
import fastcat as fc
import numpy as np 
from optparse import OptionParser
import datetime

ipath="/project/projectdirs/lsst/LSSWG/colore_raw"
opath="/project/projectdirs/lsst/LSSWG/LNMocks"
parser = OptionParser()
parser.add_option("--ipath", dest="ipath", default=ipath,
                  help="Path to colore output", type="string")
parser.add_option("--opath", dest="opath", default=opath,
                  help="Path to colore output", type="string")
parser.add_option("--N", dest="Nr", default=10,
                  help="Number of realizations", type="int")
parser.add_option("--Nstart", dest="Nstart", default=0,
                  help="starting realization", type="int")
parser.add_option("--decmin", dest="decmin", default=-70,
                  help="minimum declination", type="float")
parser.add_option("--decmax", dest="decmax", default=10,
                  help="maximum declination", type="float")
parser.add_option("--bcut", dest="bcut", default=5,
                  help="Cut in galactic b", type="float")
parser.add_option("--realspace", dest="realspace", default=False,
                  help="Realspace instead of redshift space",action="store_true")

(o, args) = parser.parse_args()

bz=np.genfromtxt(o.ipath+'/bz.txt', dtype=None, names=["z","bz"])
dNdz=np.genfromtxt(o.ipath+'/Nz.txt', dtype=None, names=["z","dNdz"])

for i in range(o.Nstart,o.Nr):
    print "Reading set ",i
    gals,inif=readColore(o.ipath+"/Set%i"%i)
    print len(gals)," galaxies read."
    if (len(gals)==0):
        print "No galaxies!"
        stop()

    # start catalog with known dNdz and bz
    N=len(gals)
    meta={}
    for k,v in inif.items():
        meta['colore_'+k]=v
    meta['mkcat_git_version']=get_git_revision_short_hash()
    meta['timestamp']=str(datetime.datetime.now())
    meta['realspace']=o.realspace
    cat=fc.Catalog(N, fields=['ra','dec','z'],dNdz=dNdz, bz=bz,
                    meta=meta)
    cat['ra']=gals['RA']
    cat['dec']=gals['DEC']
    if (o.realspace):
        cat['z']=gals['Z_COSMO']
    else:
        cat['z']=gals['Z_COSMO']+gals['DZ_RSD']

    ## first apply window function
    print "Applying window..."
    cat.setWindow(fc.window.WindowDecBcut(o.decmin, o.decmax, o.bcut),
                  apply_to_data=True)
    ## next apply photoz
    print "Applying photoz..."
    cat.setPhotoZ(fc.photoz.PhotoZGauss(0.01),apply_to_data=True)
    ## write the actual catalog
    print "Writing..."
    cat.writeH5(o.opath+'/catalog%i.h5'%(i))



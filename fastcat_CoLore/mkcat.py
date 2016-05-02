#!/usr/bin/env python
import sys
from subs import *
sys.path=["."]+sys.path
import fastcat as fc
import numpy as np 
from optparse import OptionParser
import datetime
import os
from mpi4py import MPI


#default in path
ipath="/project/projectdirs/lsst/LSSWG/colore_raw"
#default out path
opath="/project/projectdirs/lsst/LSSWG/LNMocks"
#path to humna depth maps
hpath="/project/projectdirs/lsst/LSSWG/HumnaDepthVariations"

parser = OptionParser()
parser.add_option("--ipath", dest="ipath", default=ipath,
                  help="Path to colore output", type="string")
parser.add_option("--opath", dest="opath", default=opath,
                  help="Path to colore output", type="string")
parser.add_option("--oextra", dest="oextra", default="",
                  help="Extra string to be put in output path", type="string")
parser.add_option("--hpath", dest="hpath", default=hpath,
                  help="Path to humna depth maps", type="string")
parser.add_option("--N", dest="Nr", default=10,
                  help="Number of realizations", type="int")
parser.add_option("--Nstart", dest="Nstart", default=0,
                  help="starting realization", type="int")
parser.add_option("--Ngals", dest="Ngals", default=0,
                  help="If non-zero, subsample to this number of gals", type="int")
## WF options

parser.add_option("--wftype",dest="wftype",type="string",
                  help="window func type [radecbcut,healpix]",
                  default="healpix")
parser.add_option("--wf_humnamap", dest="humnamap", type="string",
                  help="humna map type [nodither, reprandom]", default="nodither")
parser.add_option("--wf_dlogndmlim", dest="dlogndmlim", type="float",
                  help="change in number of sources as a function of depth",
                  default=0.1)
parser.add_option("--wf_decmin", dest="decmin", default=-70,
                  help="minimum declination", type="float")
parser.add_option("--wf_decmax", dest="decmax", default=10,
                  help="maximum declination", type="float")
parser.add_option("--wf_bcut", dest="bcut", default=5,
                  help="Cut in galactic b", type="float")

## PZ options

parser.add_option("--pztype",dest="pztype",type="string",
                  help="photo z type [none,gauss, doublegauss]",
                  default="gauss")
parser.add_option("--pz_sigma", dest="pz_sigma", default=0.01,
                  help="PZ: Guass sigma for (1+z)", type="float")
parser.add_option("--pz_Acat", dest="pz_Acat", default=0.2,
                  help="PZ: A catastrophic", type="float")
parser.add_option("--pz_zcat", dest="pz_zcat", default=0.4,
                  help="PZ: z catastrophic", type="float")
parser.add_option("--pz_sigmacat", dest="pz_sigmacat", default=0.05,
                  help="PZ: sigma catastrophic", type="float")


parser.add_option("--realspace", dest="realspace", default=True,
                  help="Realspace instead of redshift space",action="store_true")
parser.add_option("--ztrue", dest="ztrue", default=False,
                  help="Add column with true redshift ",action="store_true")

(o, args) = parser.parse_args()

bz=np.genfromtxt(o.ipath+'/bz.txt', dtype=None, names=["z","bz"])
dNdz=np.genfromtxt(o.ipath+'/Nz.txt', dtype=None, names=["z","dNdz"])
fopath=None

comm = MPI.COMM_WORLD
mrank = comm.Get_rank()
mranks = "["+str(mrank)+"]:"
msize = comm.Get_size()
if (mrank==0):
    print "MPI Size:", msize

out_extra=""
if len(o.oextra)>0:
    out_extra+="+"+o.oextra
if (o.realspace):
    out_extra+="+realspace"
if (o.ztrue):
    out_extra+="+ztrue"
if (o.Ngals>0):
    out_extra+="+subsamp_"+str(o.Ngals)
    
for i in range(o.Nstart,o.Nr):
    if (i%msize==mrank):
        print mranks, "Reading set ",i
        gals,inif=readColore(o.ipath+"/Set%i"%i)
        print mranks, len(gals)," galaxies read."
        if (len(gals)==0):
            print mranks, "No galaxies!"
            stop()
        # subsample if required    
        if (o.Ngals>0):
            print "Subsampling to ",o.Ngals
            # interestingly, this is superslow
            #indices=np.random.choice(xrange(len(gals)),o.Ngals, replace=False)
            # this risks repetition, but OK
            indices=np.random.randint(0,len(gals),o.Ngals)
            gals=gals[indices]
            print "Done"

        # start catalog with known dNdz and bz
        N=len(gals)
        meta={}
        for k,v in inif.items():
            meta['colore_'+k]=v
        meta['mkcat_git_version']=get_git_revision_short_hash()
        meta['timestamp']=str(datetime.datetime.now())
        meta['realspace']=o.realspace
        fields=['ra','dec','z']
        if (o.ztrue):
            fields.append('z_true')
        cat=fc.Catalog(N, fields=fields,dNdz=dNdz, bz=bz,
                        meta=meta)
        cat['ra']=gals['RA']
        cat['dec']=gals['DEC']
        if (o.realspace):
            cat['z']=gals['Z_COSMO']
        else:
            cat['z']=gals['Z_COSMO']+gals['DZ_RSD']
        if (o.ztrue):
            # make a copy of z_true various PZ algos
            # will overwrite it
            cat['z_true']=cat['z']

        ## first apply window function
        print mranks, "Creating window..."
        wfunc=getWindowFunc(o)
        print mranks, "Applying window..."
        cat.setWindow(wfunc, apply_to_data=True)

        ## next apply photoz
        pz=getPhotoZ(o)
        print mranks, "Applying photoz..."
        cat.setPhotoZ(pz,apply_to_data=True)
        ## now create full output path
        if fopath is None:
            dt=datetime.datetime.now()
            daystr="%02i%02i%02i"%(dt.year-2000,dt.month,dt.day)
            fopath=o.opath+"/"+daystr+"+"+cat.photoz.NameString()+"+"+cat.window.NameString()+out_extra
            if not os.path.exists(fopath):
                os.makedirs(fopath)
        ## write the actual catalog
        fname=fopath+'/catalog%i.h5'%(i)
        print mranks, "Writing "+fname+" ..."
        cat.writeH5(fname)

comm.Barrier()

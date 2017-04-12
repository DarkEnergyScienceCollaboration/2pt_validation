#!/usr/bin/env python
import sys
from subs import *
sys.path=["."]+sys.path
import fastcat as fc
from  fastcat.addStars import AddStars
import numpy as np
from optparse import OptionParser
import datetime
import os

if os.environ.has_key('DESC_LSS_ROOT'):
    root=os.environ['DESC_LSS_ROOT']
else:
    root="/project/projectdirs/lsst/LSSWG"
#default in path
ipath=root+"/colore_raw"
#default out path
opath=root+"/LNMocks_staging"


def main():
    o,args = setupOptions()
    initMPI(o)
    process(o)
    
def setupOptions():

    parser = OptionParser()
    parser.add_option("--params_file", dest="ipath", default=[], action="append",
                      help="Path to CoLoRe params file", type="string")
    parser.add_option("--opath", dest="opath", default=opath,
                      help="Path to output", type="string")
    parser.add_option("--oextra", dest="oextra", default="",
                      help="Extra string to be put in output path", type="string")
    parser.add_option("--Ngals", dest="Ngals", default=-1,
                      help="If non-zero, subsample to this number of gals", type="int")
    parser.add_option("--Nstars", dest="Nstars", type="int",default=0,
                      help="If Nstars>0, subsample to this number of stars. If Nstars=-1 use the full stellar catalog")
    parser.add_option("--mpi", dest="use_mpi", default=False,
                      help="If used, use mpi4py for parallelization ",action="store_true")
    ## WF and PZ options
    fc.window.registerOptions(parser)
    ## PZ options
    fc.photoz.registerOptions(parser)
    # Stars options
    AddStars.registerOptions(parser)

    ## other options
    parser.add_option("--realspace", dest="realspace", default=False,
                      help="Realspace instead of redshift space",action="store_true")
    parser.add_option("--ztrue", dest="ztrue", default=False,
                      help="Add column with true redshift ",action="store_true")

    (o, args) = parser.parse_args()
    return o,args


def getExtraString(o):
    out_extra=""
    if len(o.oextra)>0:
        out_extra+="+"+o.oextra
    if (o.realspace):
        out_extra+="+realspace"
    if (o.ztrue):
        out_extra+="+ztrue"
    if (o.Ngals>=0):
        out_extra+="+subsamp_"+str(o.Ngals)
    if(o.Nstars>0):
        out_extra+="+stars_"+str(o.Nstars)
    return out_extra


def getStars(o):
    if(o.Nstars>0):
        stars = AddStars(options=o)
    else:
        stars = None
    return stars

def initMPI(o):
    global comm, mrank, mranks, msize
    if o.use_mpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        mrank = comm.Get_rank()
        mranks = "["+str(mrank)+"]:"
        msize = comm.Get_size()
        if (mrank==0):
            print "MPI Size:", msize
    else:
        comm=None
        mrank=0
        msize=1
        mranks = ""

def process(o):
    out_extra=getExtraString(o)
    stars=getStars(o)
    do_stars = (stars is not None)
    time0 = datetime.datetime.now()
    fopath=None
    for i,param_file in enumerate(o.ipath):
        gals,inif=readColore(param_file,use_mpi=o.use_mpi)
        time1 = datetime.datetime.now()
        if mrank==0:
            print 'CoLoRe realization read. Elapsed time: ',(time1-time0).total_seconds(), ' seconds'
            dirname, _ = os.path.split(param_file)
            nzfile=inif['nz_filename']
            nzfile = os.path.join(dirname,nzfile)  
            bzfile=inif['bias_filename']
            bzfile = os.path.join(dirname,bzfile) 
            bz=np.genfromtxt(bzfile, dtype=None, names=["z","bz"])
            dNdz=np.genfromtxt(nzfile, dtype=None, names=["z","dNdz"])
            len(gals)," galaxies read."
            if (len(gals)==0):
                print "No galaxies!"
                stop()
            # subsample if required
            if (o.Ngals>=0):
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
            meta['command_line']=' '.join(sys.argv)
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
            # add star is needed
            if do_stars:
                scat=stars.generateStarCatalog(o.Nstars)
                cat.appendCatalog(scat)
            # add true z if required (apply pz will loose this)
            if (o.ztrue):
                cat['z_true']=cat['z']
            ## first apply window function
            print "Creating window..."
            wfunc=fc.window.getWindowFunc(o)
            print "Applying window..."
            cat.setWindow(wfunc, apply_to_data=True)
            ## next apply photoz
            pz=fc.photoz.getPhotoZ(o)
            print "Applying photoz..."
            cat.setPhotoZ(pz,apply_to_data=True)
            ## now create full output path
            if fopath is None:
                dt=datetime.datetime.now()
                daystr="%02i%02i%02i"%(dt.year-2000,dt.month,dt.day)
                fopath=o.opath+"/"+daystr+"+"+cat.photoz.NameString()+"+"+cat.window.NameString()+out_extra
                if not os.path.exists(fopath):
                    os.makedirs(fopath)
            ## write the actual catalog
            fname=fopath+'/fastcat_catalog%i.h5'%(i)
            print "Writing "+fname+" ..."
            cat.writeH5(fname)
        time_fin = datetime.datetime.now()
        print 'Done'
        print 'Total elapsed time: ',(time_fin-time0).total_seconds(),' seconds'


if __name__=="__main__":
    main()
    

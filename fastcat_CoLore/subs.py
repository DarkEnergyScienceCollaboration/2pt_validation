import os, sys, datetime, glob, h5py
import numpy as np
import fastcat as fc
from  fastcat.addStars import AddStars
from optparse import OptionParser


def readColoreIni(fname):
    try:
        lines=open(fname).readlines()
    except IOError:
        try:
            lines=open(fname).readlines()
        except IOError:
            print fname
            print "Could not find parameter file, giving up."
            raise IOError
    idic={}
    for line in lines:
        i=line.find('#')
        if (i>0):
            line=line[:i]
        if "= " in line:
            x,y=map(lambda x:x.strip(),line.split('= '))
            # try guessing the type
            if "." in y:
                try:
                    y=float(y)
                except:
                    pass
            if (type(y)==str):
                try:
                   y=y.split('"')[1]
                except:
                   pass
            else:
                try:
                    y=int(y)
                except:
                    pass
            idic[x]=y
            
    return idic

def readColore(params_path):
    idic=readColoreIni(params_path)
    path_out = idic['prefix_out']
    path_out = path_out+"_srcs_*.h5"
    flist=glob.glob(path_out)
    return flist,idic

def get_git_revision_short_hash():
    import subprocess
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])

def setupOptions():

    parser = OptionParser()
    if os.environ.has_key('DESC_LSS_ROOT'):
        root=os.environ['DESC_LSS_ROOT']
    else:
        root="/project/projectdirs/lsst/LSSWG"
        #default in path
    ipath=root+"/colore_raw/params.ini"
    #default out path
    opath=root+"/LNMocks_staging"
    parser.add_option("--params_file", dest="ipath", default=[], action='append',
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
    if len(o.ipath)==0:
        o.ipath.append(ipath)
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
    isc=None
    for i,param_file in enumerate(o.ipath):
        flist,inif=readColore(param_file)
        dirname, _ = os.path.split(param_file)
        nzfile=inif['nz_filename']
        nzfile = os.path.join(dirname,nzfile)
        bzfile=inif['bias_filename']
        bzfile = os.path.join(dirname,bzfile)
        bz=np.genfromtxt(bzfile, dtype=None, names=["z","bz"])
        dNdz=np.genfromtxt(nzfile, dtype=None, names=["z","dNdz"])
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
        ## create an empty catalog
        cat = fc.Catalog(0, fields=fields, dNdz=dNdz, bz=bz, meta=meta) 
        ## create window 
        wfunc=fc.window.getWindowFunc(o)
        ## next create photoz
        pz=fc.photoz.getPhotoZ(o)
        #To add the stars just once
        if do_stars & (isc is None):
            scat=stars.generateStarCatalog(o.Nstars)
            cat.appendCatalog(scat)
            isc = 0
        if fopath is None:
            dt=datetime.datetime.now()
            daystr="%02i%02i%02i"%(dt.year-2000,dt.month,dt.day)
            fopath=o.opath+"/"+daystr+"+"+pz.NameString()+"+"+wfunc.NameString()+out_extra
            if not os.path.exists(fopath):
                os.makedirs(fopath)
       
        fname=fopath+'/fastcat_catalog%i.h5'%(i)
        print "Going to write "+fname+" ..."
        if o.use_mpi:
            sizes=[0]
        for i,filename in enumerate(flist):
            if i%msize==mrank:
                print mranks, "Reading set ",i
                da=h5py.File(filename)
                data = da['sources0'].value 
                print "     ... reading : ",filename, ' with ', len(data), 'sources'
                cataux = fc.Catalog(len(data), fields=fields, dNdz=dNdz, bz=bz, meta=meta)
                cataux['ra']=data['RA']
                cataux['dec']=data['DEC']
                if (o.realspace):
                    cataux['z']=data['Z_COSMO']
                else:
                    cataux['z']=data['Z_COSMO']+data['DZ_RSD']
                # add true z if required (apply pz will loose this)
                if (o.ztrue):
                    cataux['z_true']=cataux['z']
                print "Applying window..."
                cataux.setWindow(wfunc, apply_to_data=True)
                print "Applying photoz..."
                cataux.setPhotoZ(pz,apply_to_data=True) 
                cat.appendCatalog(cataux)
                if len(cat.data)>0 & o.use_mpi:
                    sizes.append(comm.bcast(len(cat.data),root=mrank))
                t1 = datetime.datetime.now()
                print "Colore realization read. Elapsed time: ", (t1-time0).total_seconds()
         
        if o.use_mpi:    
            comm.barrier()
            if mrank==0: 
                print 'Writing... %d galaxies' % len(cat.data)
                print 'test sizes: ', sizes
                print 'len sizes: ', len(sizes) 
            with h5py.File(fname, "w",driver='mpio', comm=comm) as of:
                dset=of.create_dataset("objects", (len(cat.data),), cat.data.dtype)
                of.atomic=True
                comm.barrier()
                print 'Test sizes',sizes[0],sizes[-1]
                for i in range(len(sizes)-1):
                    if i%msize==mrank:
                        print 'sizes: ',sizes[i], sizes[i+1], i, mrank
                        #with dset.collective:
                        dset[sizes[i]:sizes[i+1]]=cat.data[sizes[i]:sizes[i+1]]
            comm.barrier()
        if mrank==0:
            if (len(cat.data)==0):
                print "No galaxies!"
                #stop()
            # subsample if required
            if (o.Ngals>=0):
                 print "Subsampling to ",o.Ngals
                 # interestingly, this is superslow
                 #indices=np.random.choice(xrange(len(gals)),o.Ngals, replace=False)
                 # this risks repetition, but OK
                 indices=np.random.randint(0,len(cat.data),o.Ngals)
                 cat.data=cat.data[indices]
                 print "Done"
        cat.setPhotoZ(pz, apply_to_data=False)
        cat.setWindow(wfunc, apply_to_data=False)
        if not o.use_mpi: 
            cat.writeH5(fname, comm)
        else:
            if type(cat.dNdz)!=type(None):
                dset=of.create_dataset("dNdz", data=cat.dNdz)
            if type(cat.bz)!=type(None):
                dset=of.create_dataset("bz", data=cat.bz)
            cat.window.writeH5(of)
            pz=of.create_dataset("photoz",data=[])
            cat.photoz.writeH5(pz)
            comm.barrier()
            of.close()
        time_fin = datetime.datetime.now()
        print 'Done'
        print 'Total elapsed time: ',(time_fin-time0).total_seconds(),' seconds'

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
                    try:
                       y=y.split('"')[1]
                    except:
                       pass
            elif type(y) is str: 
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
    flist=sorted(glob.glob(path_out))
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
    parser.add_option("--ss_frac", dest="ss_frac", default=-1,
                      help="If non-zero, subsample to this fraction of gals", type="float")
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
    if (o.ss_frac>=0):
        out_extra+="+subsamp_"+str(o.ss_frac)
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
    addFields=[]
    for i,param_file in enumerate(o.ipath):
        flist,inif=readColore(param_file)
        print flist
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
        if o.ztrue:
            addFields.append('z_true')
        ## create window 
        wfunc=fc.window.getWindowFunc(o)
        ## next create photoz
        pz=fc.photoz.getPhotoZ(o)
        if fopath is None:
            dt=datetime.datetime.now()
            daystr="%02i%02i%02i"%(dt.year-2000,dt.month,dt.day)
            fopath=o.opath+"/"+daystr+"+"+pz.NameString()+"+"+wfunc.NameString()+out_extra
            if (mrank==0):
                if not os.path.exists(fopath):
                    os.makedirs(fopath)
       
        fname=fopath+'/fastcat_catalog%i.h5'%(i)
        print "Going to write "+fname+" ..."
        Nparts=len(flist)
        for i,filename in enumerate(flist):
            if i%msize==mrank:
                print mranks, "Reading set ",i
                da=h5py.File(filename)
                data = da['sources1'].value 
                print "     ... reading : ",filename, ' with ', len(data), 'sources'
                if (o.ss_frac>=0):
                     print "Subsampling to ",o.ss_frac
                     # interestingly, this is superslow
                     #indices=np.random.choice(xrange(len(gals)),o.Ngals, replace=False)
                     # this risks repetition, but OK
                     indices=np.random.randint(0,len(data),int(len(data)*o.ss_frac))
                     data=data[indices]
                     print "Done"

                cataux = fc.Catalog(len(data), dNdz=dNdz, bz=bz, meta=meta, addFields=addFields)
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

                #Add stars
                if do_stars:
                    scat=stars.generateStarCatalog(o.Nstars/Nparts)
                    cataux.appendCatalog(scat)
                print cataux.data
                if (Nparts>0):
                    cataux.writeH5(fname, MPIComm=None,part=(i,Nparts))
                else:
                    cataux.writeH5(fname, MPIComm=None,part=None)
                    
                t1 = datetime.datetime.now()
                print "Colore realization read. Elapsed time: ", (t1-time0).total_seconds()
         

        time_fin = datetime.datetime.now()
        print 'Done'
        print 'Total elapsed time: ',(time_fin-time0).total_seconds(),' seconds'

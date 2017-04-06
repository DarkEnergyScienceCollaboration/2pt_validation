
import numpy as np
import fastcat as fc
import glob
import h5py

def readColore(params_path,use_mpi=True):
    ## first read ini file
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
            else:
                try:
                    y=int(y)
                except:
                    pass
            idic[x]=y
    data=[]
    path_out = idic['prefix_out']
    flist=glob.glob(path_out+'_srcs_*.h5")
    data=[]
    if use_mpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        mrank = comm.Get_rank()
        mranks = "["+str(mrank)+"]:"
        msize = comm.Get_size()
        if (mrank==0):
            print "MPI Size:", msize
    else:
        mrank=0
        msize=1
        mranks = ""
    for fname in flist:
        if (i%msize==mrank):
            print mranks, "Reading set ",i
            print "     ... reading : ",fname, "\r",
            da=h5py.File(fname)
            data.append(da['sources'].value)
    data=np.concatenate(data,axis=0)
    print "Read"
    return data,idic

def get_git_revision_short_hash():
    import subprocess
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])

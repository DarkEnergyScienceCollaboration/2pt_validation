
import numpy as np
import fastcat as fc

import glob 
import h5py

def readColore(path):
    ## first read ini file
    idic={}
    for line in open(path+"/params.ini").readlines():
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
    flist=glob.glob(path+"/out_*.h5")
    data=[]
    for fname in flist: 
        print "     ... reading : ",fname, "\r",
        da=h5py.File(fname)
        data.append(da['sources'].value)
    data=np.concatenate(data,axis=0)
    print "Read"
    return data,idic
        
def get_git_revision_short_hash():
    import subprocess
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])

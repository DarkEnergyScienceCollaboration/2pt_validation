import numpy as np
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
    for fname in glob.glob(path+"/out_*.h5"):
        da=h5py.File(fname)
        cdata=da['sources'].value
        if len(data)==0:
            data=cdata
        else:
            data=np.concatenate((data,cdata),axis=0)
    return data,idic
        
    
def get_git_revision_short_hash():
    import subprocess
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])

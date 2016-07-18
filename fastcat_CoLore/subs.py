
import numpy as np
import fastcat as fc

import glob
import h5py

def readColore(path):
    ## first read ini file
    idic={}
    try:
        lines=open(path+"/param.ini")
    except IOError:
        try:
            lines=open(path+"/params.ini")
        except IOError:
            print "Could not find eiter param.ini or params.ini, giving up."
            raise IOError
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

def readStars(path,zmean, zsigma):
    data=h5py.File(path)
    stars=np.array(data['data'])
    N=len(stars)
    print stars[:5]
    ## and now need to convert to same type as galaxies
    nstars=np.zeros(N,dtype=[('RA', '<f4'), ('DEC', '<f4'), ('Z_COSMO', '<f4'), ('DZ_RSD', '<f4')])
    nstars['RA']=stars['RA']
    nstars['DEC']=stars['DEC']
    stars=nstars
    if (zmean>0):
        stars['Z_COSMO']=np.random.normal(zmean, zsigma, len(stars))
    print stars[:5]

    return stars

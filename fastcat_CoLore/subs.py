
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

def getWindowFunc(o):
    if o.wftype=="radecbcut":
        wfunc= fc.window.WindowDecBcut(o.decmin, o.decmax, o.bcut)
    elif o.wftype=="healpix":
        if o.humnamap=="nodither":
            mapfn="/coaddM5Data_masked_rBand_NoDither.npz"
        elif o.humnamap=="reprandom":
            mapfn="/coaddM5Data_masked_rBand_RepulsiveRandomDitherFieldPerVisit.npz"
        else:
            print "Unknown humna type map"
            stop()
        mapfn=o.hpath+mapfn
        print "     Reading mapfn",mapfn,"..."
        hmap=np.load(mapfn)
        mask=hmap['mask']
        vals=hmap['metricValues']
        mx=vals.max()
        vals=np.exp(-o.dlogndmlim*(mx-vals))
        vals[mask]=0.0
        amask=np.where(mask==False)
        cmin,cmax,cmean=vals[amask].min(), vals[amask].max(), vals[amask].mean()
        print "Window func min, max, mean:",cmin,cmax,cmean
        info="HumnaDepthVariations map=%s dlogndmlim=%f"%(o.humnamap,o.dlogndmlim)
        shortinfo=o.humnamap+"_"+str(o.dlogndmlim)
        wfunc=fc.window.WindowHealpix(vals,info,shortinfo)
    else:
        print "Bad WF type:",o.wftype
        stop()
    return wfunc

def getPhotoZ(o):
    if o.pztype=="none":
        pz = fc.photoz.PhotoZBase()
    elif o.pztype=="gauss":
        pz = fc.photoz.PhotoZGauss(o.pz_sigma)
    elif o.pztype=="doublegauss":
        pz = fc.photoz.PhotoZDoubleGauss(o.pz_sigma,o.pz_Acat,o.pz_zcat,o.pz_sigmacat)
    elif o.pztype=="hiddenvar":
        pz = fc.photoz.PhotoZHiddenVar(o.pz_sigma)#,o.pz_zcat,o.pz_zstep)
    else:
        print "Bad PZ type:",o.pztype
        stop()
    return pz

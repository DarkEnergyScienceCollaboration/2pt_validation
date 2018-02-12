#!/usr/bin/env python
#
# Gets precision matrix from a number of files
#
import numpy as np
import sacc
import glob

globin="inputs/fastcats/GaussPZ_0.02+FullSky+run*+ztrue/twopoints*_ns2048.sacc"
fnout="inputs/fastcats/GaussPZ_0.02+FullSky+run0+ztrue_mean.sacc"

meanvecs=[]
for fn in glob.glob(globin):
    print "Loading",fn
    s=sacc.SACC.loadFromHDF(fn)
    meanvecs.append(s.mean.vector)

print "Num of files", len(meanvecs)
meanvecs=np.array(meanvecs)
tmean=meanvecs.mean(axis=0)
tvar=meanvecs.var(axis=0)
precision=1./tvar


outs=sacc.SACC(s.tracers, s.binning,sacc.MeanVec(tmean),
               sacc.Precision(matrix=np.diag(precision), mode="diagonal"),windows=s.windows)
print "Saving to ",fnout
outs.saveToHDF(fnout)


    

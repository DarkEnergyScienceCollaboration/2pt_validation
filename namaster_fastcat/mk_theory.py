import numpy as np
import pyccl as ccl
import matplotlib.pyplot as plt
import sacc
import astropy.table

def getTheories(ccl_cosmo,s,ctracers):
    theo={}
    for t1i,t2i,ells,_ in s.sortTracers(): 
        cls=ccl.angular_cl(ccl_cosmo,ctracers[t1i],ctracers[t2i],ells)
        theo[(t1i,t2i)]=cls
        theo[(t2i,t1i)]=cls
    return theo

def getTheoryVec(s, cls_theory):
    vec=np.zeros((s.size(),))
    for t1i,t2i,ells,ndx in s.sortTracers():
        vec[ndx]=cls_theory[(t1i,t2i)]
    return sacc.MeanVec(vec)
print 'Setting up cosmology'
cosmo = ccl.Cosmology(ccl.Parameters(Omega_c=0.27,Omega_b=0.045,h=0.67,sigma8=0.8,n_s=0.96,),transfer_function='eisenstein_hu',matter_power_spectrum='linear')
print 'Reading SACC file'
#SACC File with the N(z) to analyze
binning_sacc = sacc.SACC.loadFromHDF('../test/catalog0.sacc')
#Bias file (it can also be included in the SACC file in the line before)
bias_tab = astropy.table.Table.read('../test/bz_lsst.txt',format='ascii')
tracers = binning_sacc.tracers
print 'Got ',len(tracers),' tracers'
cltracers=[]
[cltracers.append(ccl.ClTracerNumberCounts(cosmo,False,False,n=(t.z,t.Nz),bias=(bias_tab['col1'],bias_tab['col2']))) for t in tracers]
print 'Cl tracers ready'
theories = getTheories(cosmo,binning_sacc,cltracers)
mean=getTheoryVec(binning_sacc,theories)
csacc=sacc.SACC(tracers,binning_sacc.binning,mean)
csacc.printInfo()
csacc.saveToHDF('theory.sacc',save_precision=False)

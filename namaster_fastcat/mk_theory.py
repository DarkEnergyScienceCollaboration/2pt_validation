import numpy as np
import pyccl as ccl
import matplotlib.pyplot as plt
import sacc
import astropy.table

hhub=0.69

def getTheories(ccl_cosmo,s,ctracers) :
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
cosmo = ccl.Cosmology(ccl.Parameters(Omega_c=0.266,Omega_b=0.049,h=hhub,sigma8=0.8,n_s=0.96,),matter_power_spectrum='linear',transfer_function='eisenstein_hu')

#Compute grid scale
zmax=2.5
ngrid=3072
a_grid=2*ccl.comoving_radial_distance(cosmo,1./(1+zmax))*(1+2./ngrid)/ngrid*hhub
print "Grid smoothing : %.3lf Mpc/h"%a_grid

print 'Reading SACC file'
#SACC File with the N(z) to analyze
binning_sacc = sacc.SACC.loadFromHDF('../test/catalog0.sacc')
#Bias file (it can also be included in the SACC file in the line before)
bias_tab = astropy.table.Table.read('../test/bz_lsst.txt',format='ascii')
tracers = binning_sacc.tracers
print 'Got ',len(tracers),' tracers'
cltracers=[ccl.ClTracer(cosmo,'nc',False,False,n=(t.z,t.Nz),bias=(bias_tab['col1'],bias_tab['col2']),r_smooth=0.5*a_grid) for t in tracers]
print 'Cl tracers ready'
theories = getTheories(cosmo,binning_sacc,cltracers)
mean=getTheoryVec(binning_sacc,theories)
csacc=sacc.SACC(tracers,binning_sacc.binning,mean)
csacc.printInfo()
csacc.saveToHDF('theory.sacc',save_precision=False)

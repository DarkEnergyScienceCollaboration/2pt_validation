#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import pyccl as ccl
import matplotlib.pyplot as plt
import sacc
import astropy.table
from optparse import OptionParser
from subs_theory import *

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

def main():
    parser = OptionParser()
    parser.add_option('--input-sacc','-i',dest="fname_in",default=None,
        help='Name of input SACC file')
    parser.add_option('--output-sacc','-o',dest="fname_out",default=None,
        help='Name of output SACC file')
    parser.add_option('--show-plot',dest="show_plot",default=False,action="store_true",
        help="Show plot comparing data and theory")
    parser.add_option('--param-file',dest="param_file",default=None,
        help="Path to CoLoRe param file")
    parser.add_option('--power-spectrum-type',dest="pk",default='halofit',type=str,
        help="Power spectrum type: [`linear`,`halofit`]")
    parser.add_option('--transfer-function',dest="tf",default='boltzmann',type=str,
        help="Type of transfer function: [`eisenstein_hu`,`bbks`,`boltzmann`,`halofit`]")
    parser.add_option("--bias-file",dest="fname_bias",default=None,
        help="Path to bias file")
    parser.add_option("--tracer-number",dest="sbin",default=4,type=int,
        help="Number of tracer to show")
    parser.add_option("--shot-noise-file",dest="fname_sn",default=None,
        help="Path to shot-noise file")
    parser.add_option("--include-rsd",dest="rsd",default=False,action="store_true",
        help="Include RSD")
    (o, args) = parser.parse_args()

    colore_dict = readColoreIni(o.param_file)
    ob = colore_dict['omega_B']
    oc = colore_dict['omega_M']-ob
    hhub = colore_dict['h']
    sigma8 = colore_dict['sigma_8']
    ns = colore_dict['ns']

    print('Reading power spectrum')
    numz=32
    numk=len(np.loadtxt('predictions/predictions_pk_srcs_pop0_z0.000.txt',unpack=True)[0])
    pk2d=np.zeros([numz,numk])
    zarr=0.05*(numz-1-np.arange(numz))
    aarr=1./(1+zarr)
    for i in np.arange(numz) :
        z=0.05*(numz-1-i)
        karr,pk,dum1,dum2=np.loadtxt('predictions/predictions_pk_srcs_pop0_z%.3lf.txt'%z,unpack=True)
        idpos=np.where(pk>0)[0];
        if len(idpos)>0 :
            idmax=idpos[-1]
            pk=np.maximum(pk,pk[idmax])
            w=1./karr**6
            pk[idmax:]=pk[idmax]*w[idmax:]/w[idmax]
        pk2d[i,:]=pk
    pk2d=pk2d.flatten()
    karr*=hhub
    pk2d/=hhub**3
    
    print('Setting up cosmology')
    cosmo = ccl.Cosmology(ccl.Parameters(Omega_c=oc,Omega_b=ob,h=hhub,sigma8=sigma8,n_s=ns,),matter_power_spectrum=o.pk,transfer_function=o.tf)
    ccl.update_matter_power(cosmo,karr,aarr,pk2d,is_linear=True)
    ccl.update_matter_power(cosmo,karr,aarr,pk2d,is_linear=False)

    #Compute grid scale
    zmax = colore_dict['z_max']
    ngrid = int(colore_dict['n_grid'])
    #rsm = colore_dict['r_smooth']/hhub 
    #shear= colore_dict['include_shear']=='true'
    #Lbox = 2*ccl.comoving_radial_distance(cosmo,1./(1+zmax))*(1+2./ngrid)
    #if shear:
    #    a_grid=Lbox/ngrid
    #else:
    #    print('No shear applied')
    #    a_grid=Lbox/ngrid
    #rsm_tot = np.sqrt(rsm**2+a_grid**2/12) 
    #print("Grid smoothing : %.3lf Mpc/h"%(rsm_tot*hhub))
    print('Reading SACC file')
    #SACC File with the N(z) to analyze
    binning_sacc = sacc.SACC.loadFromHDF(o.fname_in)
    print("Reading bias file"+o.fname_bias)
    #Bias file (it can also be included in the SACC file in the line before)
    zz,bzz=np.loadtxt(o.fname_bias,unpack=True); bzz[:]=1.
    tracers = binning_sacc.tracers
    print('Got ',len(tracers),' tracers')
    cltracers=[ccl.ClTracer(cosmo,'nc',has_rsd=o.rsd,has_magnification=False,n=(t.z,t.Nz),
                            bias=(zz,bzz)) for t in tracers]
    print('Cl tracers ready')
    theories = getTheories(cosmo,binning_sacc,cltracers)
    mean=getTheoryVec(binning_sacc,theories)
    csacc=sacc.SACC(tracers,binning_sacc.binning,mean)
    csacc.printInfo()
    csacc.saveToHDF(o.fname_out,save_precision=False)

    if o.show_plot:
        if o.sbin > len(tracers):
            print('The bin number cannot be higher than the number of tracers')
        else:
            print('Bin', o.sbin)
            plt.plot(tracers[o.sbin].z,tracers[o.sbin].Nz/np.sum(tracers[o.sbin].Nz))
            b1 = (binning_sacc.binning.binar['T1']==o.sbin) & (binning_sacc.binning.binar['T2']==o.sbin)
            xdata = binning_sacc.binning.binar['ls'][b1]
            if o.fname_sn is None:
                ydata = (binning_sacc.mean.vector[b1])*xdata*(xdata+1)
            else:
                sn_arr = np.load(o.fname_sn)
                sn = sn_arr[o.sbin,o.sbin]
                ydata = (binning_sacc.mean.vector[b1]-sn)*xdata*(xdata+1)
            b2 = (csacc.binning.binar['T1']==o.sbin) & (csacc.binning.binar['T2']==o.sbin)
            xth = csacc.binning.binar['ls'][b2]
            yth = csacc.mean.vector[b2]*xth*(xth+1)
            plt.show()
            plt.figure()
            plt.plot(xdata,ydata,label='Data')
            plt.plot(xdata,yth,label='Theory')
            plt.xlabel('$l$')
            plt.ylabel('$C_{l}l(l+1)$')
            plt.xscale('log')
            plt.ylim(-0.1,0.1)
            plt.savefig('sample_bin_%d.png'%o.sbin)
            plt.show()
    print('Done')
if __name__=="__main__":
    main()

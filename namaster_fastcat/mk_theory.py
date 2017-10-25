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
    parser.add_option('--power-spectrum-type',dest="pk",default='linear',type=str,
        help="Power spectrum type: [`linear`,`non-linear`]")
    parser.add_option('--transfer-function',dest="tf",default='eisenstein_hu',type=str,
        help="Type of transfer function: [`eisenstein_hu`,`bbks`,`boltzmann`,`halofit`]")
    parser.add_option("--bias-file",dest="fname_bias",default=None,
        help="Path to bias file")
    parser.add_option("--tracer-number",dest="sbin",default=4,type=int,
        help="Number of tracer to show")
    parser.add_option("--fsky",dest="area",default=1.,type=float,
        help="Fraction of the sky surveyed")
    (o, args) = parser.parse_args()

    colore_dict = readColoreIni(o.param_file)
    ob = colore_dict['omega_B']
    oc = colore_dict['omega_M']-ob
    hhub = colore_dict['h']
    sigma8 = colore_dict['sigma_8']
    ns = colore_dict['ns']
    print('Setting up cosmology')
    cosmo = ccl.Cosmology(ccl.Parameters(Omega_c=oc,Omega_b=ob,h=hhub,sigma8=sigma8,n_s=ns,),matter_power_spectrum=o.pk,transfer_function=o.tf)

    #Compute grid scale
    zmax = colore_dict['z_max']
    ngrid = int(colore_dict['n_grid'])
    rsm = colore_dict['r_smooth']
    shear= colore_dict['include_shear']=='true'
    if shear:
        a_grid=np.sqrt((ccl.comoving_radial_distance(cosmo,1./(1+zmax))*(1+2./ngrid)/ngrid)**2+rsm**2)
    else:
        a_grid=np.sqrt((ccl.comoving_radial_distance(cosmo,1./(1+zmax))*(1+2./ngrid)/ngrid)**2*0.5+rsm**2)

    print("Grid smoothing : %.3lf Mpc/h"%a_grid)
    print('Reading SACC file')
    #SACC File with the N(z) to analyze
    binning_sacc = sacc.SACC.loadFromHDF(o.fname_in)
    #Bias file (it can also be included in the SACC file in the line before)
    bias_tab = astropy.table.Table.read(o.fname_bias,format='ascii')
    tracers = binning_sacc.tracers
    print('Got ',len(tracers),' tracers')
    cltracers=[ccl.ClTracer(cosmo,'nc',False,False,n=(t.z,t.Nz),bias=(bias_tab['col1'],bias_tab['col2']),r_smooth=a_grid) for t in tracers]
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
            b1 = (binning_sacc.binning.binar['T1']==o.sbin) & (binning_sacc.binning.binar['T2']==o.sbin)
            xdata = binning_sacc.binning.binar['ls'][b1]
            shot_noise = 4*np.pi*o.area/np.sum(tracers[o.sbin].Nz)
            ydata = (binning_sacc.mean.vector[b1]-shot_noise)*xdata*(xdata+1)
            yth = csacc.mean.vector[b1]*xdata*(xdata+1)
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

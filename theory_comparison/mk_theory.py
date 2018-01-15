#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import pyccl as ccl
import matplotlib.pyplot as plt
import sacc
from optparse import OptionParser
import libconf as lcf
import io,os

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
        help="Power spectrum type: [`linear`,`halofit`]")
    parser.add_option('--transfer-function',dest="tf",default='boltzmann',type=str,
        help="Type of transfer function: [`eisenstein_hu`,`bbks`,`boltzmann`,`halofit`]")
    #parser.add_option("--bias-file",dest="fname_bias",default=None,
    #    help="Path to bias file")
    #parser.add_option("--shot-noise-file",dest="fname_sn",default=None, #CHECK
    #    help="Path to shot-noise file")
    parser.add_option("--include-rsd",dest="rsd",default=False,action="store_true",
        help="Include RSD")
    parser.add_option("--dz-lognormal",dest="dz_ln",default=0.05,type=float,
        help="Redshift interval to use in computation of lognormal prediction")
    (o, args) = parser.parse_args()

    #Read cosmological parameters and set cosmology object
    print("Reading CoLoRe params")
    with io.open(o.param_file) as f :
        colore_dict=lcf.load(f)
    ob=colore_dict['cosmo_par']['omega_B']
    om=colore_dict['cosmo_par']['omega_M']
    oc=om-ob
    hhub=colore_dict['cosmo_par']['h']
    ns=colore_dict['cosmo_par']['ns']
    sigma8=colore_dict['cosmo_par']['sigma_8']
    cosmo = ccl.Cosmology(ccl.Parameters(Omega_c=oc,Omega_b=ob,h=hhub,sigma8=sigma8,n_s=ns,),
                          matter_power_spectrum=o.pk,transfer_function=o.tf)

    #Calculate effective smoothing scale
    zmax=colore_dict['global']['z_max']
    ngrid=colore_dict['field_par']['n_grid']
    rsm0=colore_dict['field_par']['r_smooth']
    Lbox=2*ccl.comoving_radial_distance(cosmo,1./(1+zmax))*(1+2./ngrid)*hhub
    a_grid=Lbox/ngrid
    rsm_tot=np.sqrt(rsm0**2+a_grid**2/12.)

    #Estimate lognormal prediction
    print("Estimating lognormal prediction")
    cmd="./lnpred "
    cmd+="%lf "%om
    cmd+="%lf "%ob
    cmd+="%lf "%hhub
    cmd+="%lf "%ns
    cmd+="%lf "%sigma8
    cmd+=colore_dict['srcs1']['bias_filename']+" "
    cmd+="%lf "%rsm_tot
    cmd+="%lf "%zmax
    cmd+="%lf "%o.dz_ln
    cmd+=o.fname_out+".lnpred"
    cmd+=" > "+o.fname_out+".lnpred_log"
    os.system(cmd)

    #Read lognormal prediction
    # Note that the lognormal prediction can be negative due to numerical instabilities.
    # The lines below make sure that the power spectrum is always positive, and extrapolate
    # it beyond the range where it becomes unreliable
    numz=int(zmax/o.dz_ln)
    numk=len(np.loadtxt(o.fname_out+".lnpred_pk_z0.000.txt",unpack=True)[0])
    pk2d=np.zeros([numz,numk])
    zarr=o.dz_ln*(numz-1-np.arange(numz))
    aarr=1./(1+zarr)
    for i in np.arange(numz) :
        z=o.dz_ln*(numz-1-i)
        karr,pk,dum1,dum2=np.loadtxt(o.fname_out+".lnpred_pk_z%.3lf.txt"%z,unpack=True)
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

    #Cleanup
    os.system('rm '+o.fname_out+'.lnpred*')
    
    #Update power spectrum in cosmo to lognormal prediction (both linear and non-linear)
    ccl.update_matter_power(cosmo,karr,aarr,pk2d,is_linear=True)
    ccl.update_matter_power(cosmo,karr,aarr,pk2d,is_linear=False)

    print('Reading SACC file')
    #SACC File with the N(z) to analyze
    binning_sacc = sacc.SACC.loadFromHDF(o.fname_in)

    print("Setting up CCL tracers")
    tracers = binning_sacc.tracers
    cltracers=[ccl.ClTracer(cosmo,'nc',has_rsd=o.rsd,has_magnification=False,n=(t.z,t.Nz),
                            bias=(t.z,np.ones_like(t.z))) for t in tracers]

    print('Computing power spectra')
    theories = getTheories(cosmo,binning_sacc,cltracers)
    mean=getTheoryVec(binning_sacc,theories)
    csacc=sacc.SACC(tracers,binning_sacc.binning,mean)
    csacc.printInfo()
    csacc.saveToHDF(o.fname_out,save_precision=False)

    if o.show_plot:
        plt.figure()
        for t in tracers :
            plt.plot(t.z,t.Nz)
        plt.show()
    print('Done')
if __name__=="__main__":
    main()

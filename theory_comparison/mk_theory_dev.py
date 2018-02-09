#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from scipy.interpolate import interp1d
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
                      help="show plot comparing data and theory")
    parser.add_option('--param-file',dest="param_file",default=None,
                      help="Path to CoLoRe param file")
    parser.add_option('--power-spectrum-type',dest="pk",default='linear',type=str,
                      help="Power spectrum type: [`linear`,`halofit`]")
    parser.add_option('--transfer-function',dest="tf",default='boltzmann',type=str,
                      help="Type of transfer function: [`eisenstein_hu`,`bbks`,`boltzmann`,`halofit`]")
    parser.add_option("--include-rsd",dest="rsd",default=False,action="store_true",
                      help="Include RSD")
    parser.add_option("--dz-lognormal",dest="dz_ln",default=0.05,type=float,
                      help="Redshift interval to use in computation of lognormal prediction")
    parser.add_option("--do-lognormal",dest="do_ln",default=False,action="store_true",
                      help="Do we want to do the lognormal transformation?")
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
    rsm2=rsm_tot*rsm_tot

    #Estimate lognormal prediction
    '''
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
    cmd+=o.tf+" "
    cmd+=o.fname_out+".lnpred "
    if o.do_ln>0 :
        cmd+="1 "
    else :
        cmd+="0 "
    #cmd+="> "+o.fname_out+".lnpred_log"
    print(cmd)
    os.system(cmd)
    '''
    numz=int(zmax/o.dz_ln)
    zarr=o.dz_ln*(numz-1-np.arange(numz))
    Nk=10000
    kmin=1e-3
    kmax=50
    karr=np.array([kmin*(kmax/kmin)**(float(i)/(Nk-1)) for i in range(Nk)])
    zbias_arr,bias_arr=np.loadtxt('./inputs/BzBlue.txt',unpack=True)
    pk2d=np.zeros([numz,Nk])
    bias = interp1d(zbias_arr,bias_arr)
    for iz,z in enumerate(zarr):
        a = 1./(1+z)
        pklin = ccl.linear_matter_power(cosmo,karr,a)*(hhub**3)
        b = bias(z)
        pk = pklin*b*b*np.exp(-rsm2*karr*karr)
        r,xi = ccl.utils.pk2xi(karr,pk)
        if o.do_ln:
            xi = np.exp(xi)-1
            karr,pk = ccl.utils.xi2pk(r,xi)
            idpos=np.where(pk>0)[0];
            if len(idpos)>0 :
                idmax=idpos[-1]
                pk=np.maximum(pk,pk[idmax])
                w=1./karr**6
                pk[idmax:]=pk[idmax]*w[idmax:]/w[idmax]
            pk2d[iz,:]=pk

    pk2d=pk2d.flatten()
    karr*=hhub
    pk2d/=hhub**3
    aarr=1./(1.+zarr)
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

    lnpred_sacc=sacc.SACC.loadFromHDF("inputs/fastcats/GaussPZ_0.02+FullSky+run1+ztrue/theory_ln1_ns2048.sacc")
    
    if o.show_plot:
        plt.figure()
        # for t in tracers :
        #     plt.plot(t.z,t.Nz)
        plt.plot(mean.vector[:100],'b-')
        #plt.plot(binning_sacc.mean.vector)
        plt.plot(lnpred_sacc.mean.vector[:100],'r-')
        plt.semilogy()
        plt.show()
    print('Done')

if __name__=="__main__":
    main()

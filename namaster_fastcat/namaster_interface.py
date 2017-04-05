import numpy as np
import fastcat as fc
import sys
import os
import pymaster as nmt

fname=sys.argv[1]
nside=int(sys.argv[2])
fname_bins_z=sys.argv[3]
theta_apo=float(sys.argv[4])
templates_fname=sys.argv[5]
delta_ell=int(sys.argv[6])

#Read binning
z0_bins,zf_bins=np.loadtxt(fname_bins_z,unpack=True)
nbins=len(z0_bins)
npix=hp.nside2npix(nside)

#Read catalog
cat=fc.Catalog(read_from=fname)

#Generate bandpowers binning scheme (we're assuming all maps will use the same bandpowers!)
bpw=nmt.NmtBin(nside,nlb=delta_ell)

#Get weights, compute binary mask based on weights, and apodize it if needed
weights=cat.window
mask_noapo=np.zeros(npix);
ipix_notmasked=np.where(cat.window>0)[0]
mask_noapo[ipix_notmasked]=1.
if theta_apo>0 :
    mask_apo=nmt.mask_apodization(mask_noapo,theta_apo,apotype='Smooth')
else :
    mask_apo=mask_noapo.copy()
mask_total=mask_apo*weights

#Get contaminant templates
if templates_fname!="none" :
    templates=[[t] for t in hp.read_map(templates_fname,field=None)]
else :
    templates=[]
ntemp=len(templates)

#Generate fields
dtor=np.pi/180
fields=[]
for b in np.arange(nbins) :
    ids=np.where((cat.data['z']<zf_bins[b]) & (cat.data['z']>=z0_bins[b]))[0]
    ipix=hp.ang2pix(nside,dtor*(90-cat.data['dec'][ids]),dtor*cat.data['ra'][ids])
    mp_n=np.bincount(ipix,minlength=npix); mp_n*=mask_noapo
    mp_delta=np.zeros(npix)
    mp_delta[ipix_notmasked]=(mp_n+0.)[ip_notmasked]/weights[ip_notmasked]*(np.sum(mp_n+0.)/np.sum(weights))-1.
    if ntemp>0 :
        fields.append(nmt.NmtField(mask_total,[mp_delta],templates=templates))
    else :
        fields.append(nmt.NmtField(mask_total,[mp_delta]))
    
#Compute coupling matrix (only once assuming all maps have the same mask!)
w=nmt.NmtWorkspace()
w.compute_coupling_matri(fields[0],fields[0],bpw)

#Compute all cross-correlations
def compute_master(fa,fb,wsp,clb) :
    cl_coupled=nmt.compute_coupled_cell(fa,fb)
    cl_decoupled=wsp.decouple_cell(cl_coupled,cl_bias=clb)
    return cl_decoupled

#If attempting to deproject contaminant templates, we need an estimate of the true power spectra.
#This can be done interatively from a first guess using cl_bias=0, but I haven't coded that up yet.
#For the moment we will use cl_theory=0.
cl_theory=np.zeros(3*nside)

cls_all={}
for b1 in np.arange(nbins) :
    f1=fields[b1]
    for b2 in np.arange(nbins-b1)+b1 :
        f2=fields[b2]
        if ntemp>0 :
            cl_bias=nmt.deprojection_bias(f1,f2,w,cl_theory)
        else :
            cl_bias=None
        cls_all[(b1,b2)]=compute_master(fields[b1],fields[b2],w,clb=cl_bias)

#Missing pieces (besides all the rubbish above):
# Write cls_all into SACC format
# This should also imply adding the window functions for each bandpower using the coupling matrix encoded in w
# The effective ells for each bandpower can be obtained as bpw.get_effective_ells()
# ...


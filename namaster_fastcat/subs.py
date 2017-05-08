import numpy as np
import healpy as hp
import fastcat as fc
import sys
import pymaster as nmt
import sacc


class Mask(object) :
    """ Container for angular window function

    Attributes
    ----------
    nside : int
        HEALPix resolution parameter
    weights : np.ndarray
        Array containing fluctuations in number density of sources
    ip_nomask : np.ndarray
        Pixel ids of empty pixels
    binary : np.ndarray
        Binary mask
    mask_apo : np.ndarray
        Apodized version of the binary mask
    total : np.ndarray
        Weights map used to compute power spectra, given by mask_apo*weights
    """

    def __init__(self,cat,nside,theta_apo) :
        """ Initializes Mask object from catalog and resolution parameters
        
        Parameters
        ----------
        cat : fc.Catalog
            Catalog containing window information
        nside : int
            HEALPix resolution index
        theta_apo : float
            Apodization scale in degrees
        """
        self.weights=cat.window.map.copy()
        #Figure out which pixels are empty
        self.ip_nomask=np.where(self.weights>0)[0] 
        #Generate sky mask and apply apodization
        self.binary=np.zeros_like(self.weights)
        self.binary[self.ip_nomask]=1.
        if theta_apo>0 :
            self.mask_apo=nmt.mask_apodization(self.binary,theta_apo,apotype='C1')
        else :
            self.mask_apo=self.binary.copy()
        #Final weights are a combination of mask + depth fluctuations
        self.total=self.mask_apo*self.weights
        #Can't generate maps with pixels larger than the mask
        if nside<cat.window.nside :
            #If mask pixelization is higher, upgrade output maps
            self.nside=cat.window.nside
        else :
            #If mask pixelization is lower, upgrade mask
            self.nside=nside
            self.total  =hp.ud_grade(self.total  ,nside_out=nside)
            self.binary =hp.ud_grade(self.binary ,nside_out=nside)
            self.weights=hp.ud_grade(self.weights,nside_out=nside)
            self.ip_nomask=np.where(self.weights>0)[0] 


class Tracer(object) :
    """Object containing all the information about a given sky tracer
    
    Attributes
    ----------
    zarr,nzarr : np.ndarray
        Arrays defining the tracer's selection function
    lmax : float
        Maximum multipole for which we can trust this tracer
    nside : int
        HEALPix pixelization index
    field : NaMaster Field object containing sky fluctuations
    """

    #TODO: maybe we want to do something more general than top-hat binning in z
    def __init__(self,cat,z0,zf,lmax,mask,n_sampling=256,dz_sampling=-1,templates=None) :
        """Initialize tracer (get sample, compute N(z) and sky map
        
        Parameters
        ----------
        cat : fc.Catalog
            Catalog with object properties
        z0,zf : float
            Edges of the redshift bin defining the sample
        lmax : int
            Maximum multipole for this tracer
        mask : Mask
            Mask object map (this determines the resolution of the sky maps)
        n_sampling : int 
            Number of points where you want to sample the N(z)
        dz_sampling : float
            Redshift resolution of the sampling (only used if n_sampling=None)
        templates : np.ndarray
            Sets of maps of contaminant templatest to be subtracted
        """

        #Store lmax
        self.lmax=lmax

        #Bin by z_Mean
        #TODO: would be better to use ML, but fastcat doesn't have that out of the box
        zm=cat.photoz.getMeanRMS(cat.data)[0]
        ids=np.where((zm<zf) & (zm>=z0))[0]
        data_here=cat.data[ids]
        
        #Get interval where N(z) will be significant
        zmin_arr,zmax_arr=cat.photoz.getMinMax(data_here)
        zmin=np.amin(zmin_arr); zmax=np.amax(zmax_arr)

        #Get N(z) as stacked pdfs
        if n_sampling==None :
            dz=dz_sampling
        else :
            dz=(zmax-zmin)/n_sampling
        self.zarr,self.nzarr=cat.photoz.NofZ(data_here,zmin,zmax,dz)

        #Get sky map
        #1- Get number of objects per pixel
        dtor=np.pi/180
        self.nside=mask.nside
        npix=hp.nside2npix(self.nside)
        ipix=hp.ang2pix(self.nside,dtor*(90-data_here['dec']),dtor*data_here['ra'])
        mp_n=np.bincount(ipix,minlength=npix).astype(float)
        
        #2- Identify masked pixels
        mp_n*=mask.binary

        #3- Compute overdensity field
        #In order to transform the number counts map into an overdensity we make the following assumption:
        #  - The mean density in each pixel is n_mean_i \propto weights_i
        #  - In that case delta_i = n_i/n_mean_i-1 = norm * n_i/weights_i -1
        #    where the normalization factor is norm = <weights>/<n>
        #TODO: if parts of the map are contaminated it may be worth passing a set of pixels
        #      defining a clean region that should be used to compute the mean density.
        #      currently we use the whole (unmasked) sky.
        norm=np.sum(mp_n)/np.sum(mask.weights)
        mp_delta=np.zeros(npix)
        mp_delta[mask.ip_nomask]=norm*(mp_n[mask.ip_nomask]/mask.weights[mask.ip_nomask])-1.

        #4- Generate NaMaster fields (passing contaminants if necessary)
        if templates is None :
            self.field=nmt.NmtField(mask.total,[mp_delta],templates=templates)
        else :
            self.field=nmt.NmtField(mask.total,[mp_delta])


def process_catalog(fname_catalog,fname_bins,nside,fname_out,apodization_scale=0.,
                    fname_templates="none",bins_ell=4) :
    #Read z-binning
    print "Bins"
    z0_bins,zf_bins,lmax_bins=np.loadtxt(fname_bins,unpack=True)
    nbins=len(z0_bins)


    #Read catalog
    print "Catalog"
    cat=fc.Catalog(read_from=fname_catalog)


    #Get weights, compute binary mask based on weights, and apodize it if needed
    print "Window"
    mask=Mask(cat,nside,apodization_scale)
    nside=mask.nside


    #Get contaminant templates
    #TODO: check resolution
    if fname_templates!="none" :
        templates=[[t] for t in hp.read_map(fname_templates,field=None)]
        ntemp=len(templates)
    else :
        templates=None
        ntemp=0


    #Generate bandpowers binning scheme (we're assuming all maps will use the same bandpowers!)
    print "Bandpowers"
    bpw=nmt.NmtBin(nside,nlb=bins_ell)
    ell_eff=bpw.get_effective_ells()


    print "Subsampling"
    #Generate tracers
    dtor=np.pi/180
    npix=hp.nside2npix(nside)
    tracers=[]
    for b in np.arange(nbins) :
        tracers.append(Tracer(cat,z0_bins[b],zf_bins[b],lmax_bins[b],mask,templates=templates))


    print "Compute power spectra"
    #Compute coupling matrix
    #TODO: (only done once, assuming all maps have the same mask!)
    print "  Computing coupling matrix"
    w=nmt.NmtWorkspace()
    w.compute_coupling_matrix(tracers[0].field,tracers[0].field,bpw)

    #Compute all cross-correlations
    def compute_master(fa,fb,wsp,clb) :
        cl_coupled=nmt.compute_coupled_cell(fa,fb)
        cl_decoupled=wsp.decouple_cell(cl_coupled,cl_bias=clb)
        return cl_decoupled

    #If attempting to deproject contaminant templates, we need an estimate of the true power spectra.
    #This can be done interatively from a first guess using cl_bias=0, but I haven't coded that up yet.
    #For the moment we will use cl_guess=0.
    cl_guess=np.zeros(3*nside)

    print "  Computing power spectrum"
    cls_all={}
    for b1 in np.arange(nbins) :
        f1=tracers[b1].field
        for b2 in np.arange(b1,nbins) :
            f2=tracers[b2].field
            if ntemp>0 :
                cl_bias=nmt.deprojection_bias(f1,f2,w,cl_theory)
            else :
                cl_bias=None
            cls_all[(b1,b2)]=compute_master(f1,f2,w,clb=cl_bias)[0]


    print "Translating into SACC"
    #Transform everything into SACC format
    #1- Generate SACC tracers
    stracers=[sacc.Tracer("tr_b%d"%i,"point",t.zarr,t.nzarr,exp_sample="gals")
              for i,t in enumerate(tracers)]

    #2- Define SACC binning
    typ,ell,t1,q1,t2,q2=[],[],[],[],[],[]
    for i1 in np.arange(nbins) :
        for i2 in np.arange(i1,nbins) :
            lmax=min(tracers[i1].lmax,tracers[i2].lmax)
            for l in ell_eff[ell_eff<lmax] :
                typ.append('F')
                ell.append(l)
                t1.append(i1); t2.append(i2)
                q1.append('P'); q2.append('P')
    sbin=sacc.Binning(typ,ell,t1,q1,t2,q2)
    ssbin=sacc.SACC(stracers,sbin)

    #3- Arrange power spectra into SACC mean vector
    vec=np.zeros((ssbin.size(),))
    for t1i,t2i,ells,ndx in ssbin.sortTracers() :
        lmax=min(tracers[t1i].lmax,tracers[t2i].lmax)
        vec[ndx]=cls_all[(t1i,t2i)][np.where(ell_eff<lmax)[0]]
    svec=sacc.MeanVec(vec)

    #4- Create SACC file and write to file
    csacc=sacc.SACC(stracers,sbin,svec)
    csacc.saveToHDF(fname_out)

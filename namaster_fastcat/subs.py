import numpy as np
import healpy as hp
import fastcat as fc
import sys
import pymaster as nmt
import sacc
import matplotlib.pyplot as plt
import h5py
from optparse import OptionParser
import os
from time import time
from numpy.linalg import inv
import astropy.table
debug=False
test_photoz=False
def initMPI(o):
    global comm, mrank, mranks, msize 
    
    if o.use_mpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        mrank = comm.Get_rank()
        mranks = "["+str(mrank)+"]:"
        msize = comm.Get_size()
        if (mrank==0):
            print "MPI Size:", msize
    else:
        comm=None
        mrank=0
        msize=1
        mranks = ""

def setupOptions():
    parser = OptionParser()

    parser.add_option("--input-file", dest="fname_in",default=None,
                      help="Path to fastcat input", type="string")
    parser.add_option("--output-file", dest="fname_out", default=None,
                      help="Path to output", type="string")
    parser.add_option("--nside", dest="nside", default=2048,
                      help="Nside resolution of maps", type="int")
    parser.add_option("--nz-bins-file", dest="fname_bins_z", default=None,
                      help="Name of the binning file", type="string")
    parser.add_option("--theta-apo", dest="theta_apo", type="float",default=0.,
                      help="Apodization angle")
    # Right now the --templates option doesn't recognize None so as a temporary
    # solution I pass the "default" string "none".
    parser.add_option("--templates", dest="templates_fname", default="none",
                      type="string",help="Templates to subtract from power-spectrum")
    parser.add_option("--delta-ell", dest="delta_ell", default=50,type="int",
                      help="Width of ell binning")
    parser.add_option("--mpi", dest="use_mpi", default=False,
                      help="If used, use mpi4py for parallelization ",action="store_true")
    parser.add_option("--nmt-workspace",dest="nmt_workspace",default="none",
                      help="Path to file containing the NaMaster workspace for this window function",type="string")
    parser.add_option("--save-map",dest="save_map",default=False,action="store_true",
                      help="Save input maps to namaster")
    parser.add_option("--compute-covariance",dest="compute_covariance",default=False,action="store_true",
                      help="Compute the theoretical covariance matrix using the predicted Cls")
    parser.add_option("--subtract-shot-noise",dest="sub_sn",default=False,action="store_true",
                      help="Subtract shot-noise to NaMaster results, if not used the shot-noise will be saved in a file")
    (o, args) = parser.parse_args()
    return o,args


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
        
        if cat.window.typestr=='decbcut':
            npix=hp.nside2npix(nside)
            ind_aux=np.arange(0,npix)
            theta,phi=hp.pix2ang(nside,ind_aux)
            self.weights=cat.window(phi*180/np.pi,90-theta*180/np.pi)
        elif cat.window.typestr=='base':
            self.weights=np.ones(12*nside**2)
        else:
            self.weights=cat.window.map.copy()
        if debug:
            hp.mollview(self.weights)
            plt.show() 
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
        if (cat.window.typestr!='decbcut') & (cat.window.typestr!='base'):
            if nside<cat.window.nside:
                #If mask pixelization is higher, upgrade output maps
                self.nside=cat.window.nside
            else :
                #If mask pixelization is lower, upgrade mask
                self.nside=nside
                self.total  =hp.ud_grade(self.total  ,nside_out=nside)
                self.binary =hp.ud_grade(self.binary ,nside_out=nside)
                self.weights=hp.ud_grade(self.weights,nside_out=nside)
                self.ip_nomask=np.where(self.weights>0)[0] 
        else:
            self.nside=nside
            self.total  =hp.ud_grade(self.total  ,nside_out=nside)
            self.binary =hp.ud_grade(self.binary ,nside_out=nside)
            self.weights=hp.ud_grade(self.weights,nside_out=nside)
            self.ip_nomask=np.where(self.weights>0)[0]
        #Pixel area
        self.pixOmega=4*np.pi/hp.nside2npix(self.nside)
        if debug:
            hp.mollview(self.weights)
            plt.show()
      
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
    def __init__(self,map_n,zarr,nzarr,lmax,mask,templates=None) :
        """Initialize tracer (get sample, compute N(z) and sky map
        
        Parameters
        ----------
        map_n : np.ndarray
            Map with number of objects per pixel
        zarr : np.ndarray
            Array of redshifts
        nzarr : np.ndarray 
            N(z) for each z in zarr
        lmax : int
            Maximum multipole for this tracer
        mask : Mask
            Mask object map (this determines the resolution of the sky maps)
        templates : np.ndarray
            Sets of maps of contaminant templatest to be subtracted
        """

        #Store lmax
        self.lmax=lmax

        #Store N(z)
        self.zarr=zarr
        self.nzarr=nzarr

        #Angular resolution
        self.nside=mask.nside
        npix=hp.nside2npix(self.nside)
        npix_unmasked=len(mask.ip_nomask)+0.
        if hp.npix2nside(len(map_n))!=self.nside :
            raise ValueError("Mask and map have different resolutions")

        #Get sky map
        #1- Get number of objects per pixel
        mp_n=map_n.copy()
            
        #2- Identify masked pixels
        mp_n*=mask.binary

        #3- Compute overdensity field
        #In order to transform the number counts map into an overdensity we make the following assumption:
        #  - The mean density in each pixel is n_mean_i \propto weights_i
        #  - In that case delta_i = n_i/n_mean_i-1 = n_i/(nmean * weights_i) -1
        #    where the normalization factor is nmean = <n/weights>
        #TODO: if parts of the map are contaminated it may be worth passing a set of pixels
        #      defining a clean region that should be used to compute the mean density.
        #      currently we use the whole (unmasked) sky.
        nmean=np.mean(mp_n[mask.ip_nomask]/mask.weights[mask.ip_nomask])
        mp_delta=np.zeros(npix)
        mp_delta[mask.ip_nomask]=mp_n[mask.ip_nomask]/(nmean*mask.weights[mask.ip_nomask])-1.
        self.shotnoise=mask.pixOmega**2*np.sum(1./mask.weights[mask.ip_nomask])/(4*np.pi*nmean)
        if debug:
            print np.max(mp_delta[mask.ip_nomask]), np.min(mp_delta[mask.ip_nomask]), '; Max= ', np.max(mp_n[mask.ip_nomask]/mask.weights[mask.ip_nomask]), ' Mean= ',nmean
            hp.mollview(mp_delta)
            plt.show()
            ma_delta = hp.ma(mp_delta)
            ma_delta.mask = np.logical_not(mask)
            cl = hp.anafast(ma_delta.filled(), lmax=int(lmax))
            plt.figure()
            plt.plot(np.arange(len(cl)),cl)
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel(r'$l$')
            plt.ylabel(r'$C_{l}$')
            plt.show()
            tab = astropy.table.Table([np.arange(len(cl)),cl],names=('l','Cl'))
            tab.write('cl_anafast_%f.fits' %(np.average(zarr,weights=nzarr))) 
        #4- Generate NaMaster fields (passing contaminants if necessary)
        if templates is not None :
            self.field=nmt.NmtField(mask.total,[mp_delta],templates=templates)
        else :
            self.field=nmt.NmtField(mask.total,[mp_delta])
         
#Returns map and N(z)
def bin_catalog(cat,z0_arr,zf_arr,mask,zmin=0,zmax=4.,n_sampling=1024,dz_sampling=-1,fac_neglect=1E-5) : 
    try:
        nbins=len(z0_arr)
    except:
        z0_arr = [z0_arr]
        zf_arr = [zf_arr]
        nbins=len(z0_arr)
    #Get N(z) as stacked pdfs
    if n_sampling==None :
        dz=dz_sampling
    else :
        dz=(zmax-zmin)/n_sampling
    zarr_single=np.arange(zmin,zmax,dz)
    numz=len(zarr_single)

    nzarr_all_here=np.zeros([nbins,numz])
    npix=hp.nside2npix(mask.nside)
    maps_all_here=np.zeros([nbins,npix])
    dtor=np.pi/180
    t0 = time()

    #Compute a map with #part per pix and an N(z) for all the files corresponding to this node
    for i in np.arange(cat.Npart):
        if (i%msize==mrank):
            cat.readNextPart(i)
            zm=cat.photoz.getMeanRMS(cat.data)[0]
            for ib in np.arange(nbins) :
                z0=z0_arr[ib]; zf=zf_arr[ib];
                ids=np.where((zm<zf) & (zm>=z0))[0]
                data_here=cat.data[ids]
                try:
                    zarr,nzarr=cat.photoz.NofZ_true(data_here,zmin,zmax,dz)
                except:
                    print "NO"
                    exit(1)
                    zarr,nzarr=cat.photoz.NofZ(data_here,zmin,zmax,dz)
                ipix=hp.ang2pix(mask.nside,dtor*(90-data_here['dec']),dtor*data_here['ra'])
                mp_n=np.bincount(ipix,minlength=npix).astype(float)
                nzarr_all_here[ib,:]+=nzarr
                maps_all_here[ib,:]+=mp_n
            print 'Created map ', i, ' elapsed time, ', time()-t0, ' seconds'
    #Now reduce maps and N(z)'s for all nodes into a common one
    #Only master node will do stuff afterwards
    if msize==1 :
        maps_all=maps_all_here
        nzarr_all=nzarr_all_here
    else :
        if mrank==0 :
            maps_all=np.zeros_like(maps_all_here)
            nzarr_all=np.zeros_like(nzarr_all_here)
        else :
            maps_all=None;
            nzarr_all=None
        comm.Reduce(maps_all_here,maps_all)
        comm.Reduce(nzarr_all_here,nzarr_all)

    if mrank!=0 :
        return None,None,None

    #Find relevant redshift range for each bin
    zarr_out=[]
    nzarr_out=[]
    for n in nzarr_all :
        mx=np.amax(n)
        #Find indices that are below threshold
        ids=np.where(n<mx*fac_neglect)[0]
        if len(ids)<=0 : #If all entries are above threshold
            iz0=0; izf=len(n)-1;
        else :
            iz0=0
            if ids[0]==0 : #If first entry is already above threshold just keep going
                while ((iz0<len(ids)) and (ids[iz0+1]==ids[iz0]+1)) :
                    iz0+=1
            if iz0==len(ids)-1 :
                raise ValueError("Empty bin")

            izf=len(ids)-1
            if ids[izf]==izf : #If last entry is already above threshold just keep going
                while((izf>=0) and (ids[izf-1]==ids[izf]-1)) :
                    izf-=1
            if izf==0 :
                raise ValueError("Empty bin")

            #Save only range of redshifts with relevant N(z)
            #TODO: a different compression scheme might be more efficient
            zarr_out.append(zarr_single[iz0:izf+1])
            nzarr_out.append(n[iz0:izf+1])

    return zarr_out,nzarr_out,maps_all

def compute_covariance(w,clpred,binning,t1,t2,t3,t4,nz1,nz2,cw):
    """Routine to compute the covariance matrix using NaMaster
    needs a NaMaster workspace w, the 4 tracers considered, and an array with the predicted cls
    cl_t1t3, cl_t1t4, cl_t2t3, cl_t2t4.
    """
    
    t1t3 = np.logical_and(binning.binar['T1']==min(t1,t3),binning.binar['T2']==max(t3,t1))
    t1t4 = np.logical_and(binning.binar['T1']==min(t1,t4),binning.binar['T2']==max(t4,t1))
    t2t3 = np.logical_and(binning.binar['T1']==min(t2,t3),binning.binar['T2']==max(t2,t3))
    t2t4 = np.logical_and(binning.binar['T1']==min(t2,t4),binning.binar['T2']==max(t4,t2))

    if t1==t3:
        c13 = clpred[t1t3]+1./nz1
    else:
        c13 = clpred[t1t3]
    if t1==t4:
        c14 = clpred[t1t4]+1./nz1
    else:
        c14 = clpred[t1t4]
    if t2==t3:
        c23 = clpred[t2t3]+1./nz2
    else:
        c23 = clpred[t2t3]
    if t2==t4:
        c24 = clpred[t2t4]+1./nz2
    else:
        c24 = clpred[t2t4]  
    return nmt.gaussian_covariance(cw,c13,c14,c23,c24)

def process_catalog(o) :

    #Read z-binning
    print "Bins"
    z0_bins,zf_bins,lmax_bins=np.loadtxt(o.fname_bins_z,unpack=True)
    try:
        nbins=len(z0_bins)
    except:
        nbins=1

    cat=fc.Catalog(read_from=o.fname_in)
    

    #Get weights, compute binary mask based on weights, and apodize it if needed
    print "Window"
    mask=Mask(cat,o.nside,o.theta_apo)
    nside=mask.nside
    tot_area=4.*np.pi*np.sum(mask.weights)/len(mask.weights)

    #Get contaminant templates
    if "none" not in o.templates_fname :
        templates=[[t] for t in hp.read_map(o.templates_fname,field=None)]
        ntemp=len(templates)
    else :
        templates=None
        ntemp=0

    #Generate bandpowers binning scheme (we're assuming all maps will use the same bandpowers!)
    print "Bandpowers"
    bpw=nmt.NmtBin(nside,nlb=o.delta_ell)
    ell_eff=bpw.get_effective_ells()
    tracers=[]

    #Generate tracers
    print "Maps"
    #TODO: pass extra sampling parameters
    zs,nzs,mps=bin_catalog(cat,z0_bins,zf_bins,mask)
    if mrank!=0 :
        return
    #Make sure that we have an iterable
    try:
        len(lmax_bins)
    except:
        lmax_bins=[lmax_bins]
    for zar,nzar,mp,lmax in zip(zs,nzs,mps,lmax_bins):
        zav = np.average(zar,weights=nzar)
        print "-- z-bin: %3.2f "%zav
        tracers.append(Tracer(mp,zar,nzar,lmax,mask,templates=templates))
        if o.save_map:
            hp.write_map("map_%3.2f.fits"%zav,mp)
        cat.rewind()
        
    print "Compute power spectra"
    #Compute coupling matrix
    #TODO: (only done once, assuming all maps have the same mask!)
    print "  Computing coupling matrix"
    w=nmt.NmtWorkspace()
    if not(os.path.isfile(o.nmt_workspace)) :
        w.compute_coupling_matrix(tracers[0].field,tracers[0].field,bpw)
        if o.nmt_workspace!="none" :
            w.write_to(o.nmt_workspace)
    else :
        w.read_from(o.nmt_workspace)

    #Compute all cross-correlations
    def compute_master(fa,fb,wsp,clb=None,cln=None) :
        cl_coupled=nmt.compute_coupled_cell(fa,fb)
        cl_decoupled=wsp.decouple_cell(cl_coupled,cl_bias=clb,cl_noise=cln)
        return cl_decoupled

    #If attempting to deproject contaminant templates, we need an estimate of the true power spectra.
    #This can be done interatively from a first guess using cl_bias=0, but I haven't coded that up yet.
    #For the moment we will use cl_guess=0.
    cl_guess=np.zeros(3*nside)
    t1 = time()
    print "  Computing power spectrum"
    cls_all={}
    for b1 in np.arange(nbins) :
        f1=tracers[b1].field
        for b2 in np.arange(b1,nbins) :
            f2=tracers[b2].field
            if (b1==b2) and o.sub_sn :
                cl_noise=tracers[b1].shotnoise*np.ones([1,3*nside])
            else :
                cl_noise=None
            if ntemp>0 :
                cl_bias=nmt.deprojection_bias(f1,f2,w,cl_guess)
            else :
                cl_bias=None
            cls_all[(b1,b2)]=compute_master(f1,f2,w,clb=cl_bias,cln=cl_noise)[0]

        print 'Computed bin: ', b1, b2, ' in ', time()-t1, ' s'
        if debug:
            plt.figure()
            plt.plot(ell_eff,cls_all[(b1,b1)])
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel(r'$l$')
            plt.ylabel(r'$C_{l}$')
            plt.show()
  
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
    cw = nmt.covariance.NmtCovarianceWorkspace()
    cw.compute_coupling_coefficients(w,w)

    #4- Create SACC file and write to file
    csacc=sacc.SACC(stracers,sbin,svec)

    #5- Compute covariance if needed
    if o.compute_covariance:
        print "Computing covariance"
        sacc_th = csacc
        cov_all={}
        tcov0 = time()
        #This doesn't avoid some repetitions but it is the simplest way
        ngar=np.array([np.sum(t.nzarr)/tot_area for t in tracers])
        for i1 in np.arange(nbins):
             for i2 in np.arange(i1,nbins):
                 for i3 in np.arange(nbins):
                     for i4 in np.arange(i3,nbins): 
                         cov_all[(i1,i2,i3,i4)]=compute_covariance(w,sacc_th.mean.vector,sacc_th.binning,
                                                                   i1,i2,i3,i4,ngar[i1],ngar[i2],cw)
                         cov_all[(i2,i1,i4,i3)]=cov_all[(i1,i2,i3,i4)]
        cov=np.zeros((ssbin.size(),ssbin.size()))
        for t1i,t2i,ells,ndx in ssbin.sortTracers():
            for t3i,t4i,ells2,ndy in ssbin.sortTracers():
                lmax=min(tracers[t1i].lmax,tracers[t2i].lmax,tracers[t3i].lmax,tracers[t4i].lmax)
                cov[ndx,ndy]=cov_all[(t1i,t2i,t3i,t4i)][ell_eff<lmax,ell_eff<lmax]
        icov=inv(cov)
        precision=sacc.Precision(icov,"dense",sbin)
        csacc=sacc.SACC(stracers,sbin,svec,precision)
        print 'Computed covariance in', time()-tcov0, ' seconds'
    csacc.saveToHDF(o.fname_out)

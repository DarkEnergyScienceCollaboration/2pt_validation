import astropy.table
import glob
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt


#Routine to read astropy table and create healpix maps with some weights (to get E1, E2 maps)
def make_hp_map(nside,data,weights):
    """Routine to read astropy tables and create
    healpix maps

    Args:
    ----
        nside :int: Nside parameter of healpix, the total number of pixels
        will be 12*nside**2. Nside should be a power of 2.
        data :astropy.table.Table: table or dictionary of numpy arrays 
        containing the field RA and DEC in degrees.
        weights :np.array: Array with the same length as the table with
        the contains that you want to average out in each healpix pixel.

    Returns:
    -------
         map_gal :np.array(size=12*nside**2): Healpix output map
    """
    pix_nums = hp.ang2pix(nside,np.pi/2-data['DEC']*np.pi/180,data['RA']*np.pi/180)
    parts_per_pix = np.bincount(pix_nums)
    bin_count = np.bincount(pix_nums,weights=weights)
    w = np.zeros_like(bin_count)
    w[parts_per_pix>0] = bin_count[parts_per_pix>0]*1.0/parts_per_pix[parts_per_pix>0]
    map_gal = np.append(w,np.zeros(12*nside**2-len(w)))
    return map_gal

def make_hp_map_counts(nside,data):
    """Routine to read astropy tables and create
    healpix maps with number of counts

    Args:
    ----
        nside :int: Nside parameter of healpix, the total number of pixels
        will be 12*nside**2. Nside should be a power of 2.
        data :astropy.table.Table: table or dictionary of numpy arrays 
        containing the field RA and DEC in degrees.  

    Returns:
    -------
         map_gal :np.array(size=12*nside**2): Healpix output map
    """
    pix_nums = hp.ang2pix(nside,np.pi/2-data['DEC']*np.pi/180,data['RA']*np.pi/180)
    parts_per_pix = np.bincount(pix_nums) 
    map_gal = np.append(parts_per_pix,np.zeros(12*nside**2-len(parts_per_pix)))
    return map_gal

def create_table(filelist,outname,subsample=50,zmin=0.5,zmax=0.7):
    """ Routine to create an astropy table from the simulation filelist

    Args:
    ----
        filelist :list: list of strings containing the path to the different
        tables to read.
        outname :str: Name of the output file.
        subsample :int: Subsampling factor that you want to use to speed up
        the table reading (default: 50)
        zmin :float: Minimum redshift of the selected sample (default: 0.5).
        zmax :float: Maximum redshift of the selected sample (default: 0.7).
    
    Returns:
    -------
        tab_aux :astropy.table.Table: Table of the subsample.
    """
    for i,filename in enumerate(filelist):
        print 'Reading file', i, ' of ', len(filelist),'...'
        tab = astropy.table.Table.read(filename)
        random_sel = np.random.random_integers(0,len(tab),size=int(len(tab)/subsamp))
        if(i==0):
            tab_aux=tab[random_sel]
        else:
            tab_aux = astropy.table.vstack([tab_aux,tab[random_sel]])

    zbin = np.logical_and(tab_aux['Z_COSMO']+tab_aux['DZ_RSD']>zmin,tab_aux['Z_COSMO']+tab_aux['DZ_RSD']<zmax)

    tab_aux[zbin].write(outname)
    return tab_aux[zbin]
################################

#List of files to read
filelist = glob.glob('/global/cscratch1/sd/jsanch87/massmap_hsc_4096_srcs_*')
read_files = False
nside_out = 64

if read_files:
    subsamp = 5 # I am going to just take one every 50 galaxies
    tab_aux = create_table(filelist,'/global/cscratch1/sd/jsanch87/massmap_hsc_4096_0.5_0.7_tab.fits.gz')

    dummy_map = make_hp_map(nside_out,tab_aux[zbin],np.ones(np.count_nonzero(zbin)))
    e1_map = make_hp_map(nside_out,tab_aux[zbin],tab_aux['E1'][zbin])
    e2_map = make_hp_map(nside_out,tab_aux[zbin],tab_aux['E2'][zbin])
else:
    tab = astropy.table.Table.read('/global/cscratch1/sd/jsanch87/massmap_hsc_4096_0.5_0.7_tab.fits.gz')
    tab_aux = tab[np.logical_and(tab['Z_COSMO']>0.5,tab['Z_COSMO']<0.7)]
    tab_aux2 = tab
    dummy_map = make_hp_map_counts(nside_out,tab_aux)
    hp.mollview(dummy_map)
    plt.show()
    dummy_map2 = make_hp_map_counts(nside_out,tab_aux2)
    hp.mollview(dummy_map2)
    plt.show()
    plt.figure()
    plt.scatter(dummy_map,dummy_map2)
    plt.show()
    #plt.figure()
    #plt.hist(dummy_map,bins=15,range=(0,15),label='narrow bin',alpha=0.5)
    #plt.hist(dummy_map2,bins=15,range=(0,15),label='wide bin',alpha=0.5)
    #plt.legend(loc='best')
    #plt.yscale('log')
    #plt.show()
    e1_map = make_hp_map(nside_out,tab_aux,tab_aux['E1'])
    hp.mollview(e1_map)
    #plt.show()
    e2_map = make_hp_map(nside_out,tab_aux,tab_aux['E2'])
    hp.mollview(e2_map)
    #plt.show()
tlm, elm, blm = hp.map2alm([np.random.random(len(e1_map)),e1_map,e2_map],pol=True)
kappa = hp.alm2map(elm,nside_out,pol=False)
hp.mollview(kappa,title="Kappa reconstructed")
plt.show()
hp.write_map("test_kappa_0.5_0.7.fits",kappa)
map_orig = hp.read_map('/global/cscratch1/sd/jsanch87/massmap_hsc_4096_kappa_z002.fits')
map_orig = hp.ud_grade(map_orig,nside_out)
hp.mollview(map_orig,title="Kappa original")
plt.show()
cl_orig = hp.anafast(map_orig)
cl_rec = hp.anafast(kappa)
plt.plot(4.*cl_orig/cl_rec)
plt.ylim(0.8,1.2)
plt.xscale('log')
plt.xlabel('$l$')
plt.ylabel('$C_{l,orig}/C_{l,reco}$')
plt.show()


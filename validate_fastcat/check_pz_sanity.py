#!/usr/bin/env python
import sys
sys.path.append('./fastcat_CoLore')
from subs import *
import fastcat as fc
from optparse import OptionParser
import pylab

parser = OptionParser(usage="""

validate_fastcat/check_pz_sanity.py HDF5_fastcat_file

This code will check the sanity of PZ in a given file by computing the
cumulative PDF for PZ at the true redshift and histograming it.

Hence you need to have true redshift available (option --ztrue in mkcat.py)
""")

## PZ options

parser.add_option("--Nz", dest="Nz", default=10000,
                  help="Number of zs to consider", type="int")

parser.add_option("--zmax", dest="zmax", default=1.5,
                  help="maxz", type="float")


parser.add_option("--pztype",dest="pztype",type="string",
                  help="photo z type [none,gauss, doublegauss, hiddenvar]",
                  default="gauss")
parser.add_option("--pz_sigma", dest="pz_sigma", default=0.01,
                  help="PZ: Guass sigma for (1+z)", type="float")
parser.add_option("--pz_Acat", dest="pz_Acat", default=0.2,
                  help="PZ: A catastrophic", type="float")
parser.add_option("--pz_zcat", dest="pz_zcat", default=0.3,
                  help="PZ: z catastrophic", type="float")
parser.add_option("--pz_zstep", dest="pz_zstep", default=0.6,
                  help="PZ: z step if hiddenvar", type="float")
parser.add_option("--pz_sigmacat", dest="pz_sigmacat", default=0.05,
                  help="PZ: sigma catastrophic", type="float")

(o, args) = parser.parse_args()

print "Creating catalog..."
N=o.Nz
cat=fc.Catalog(N,fields=['z','z_true'])
cat['z']=np.linspace(0,o.zmax,N)
cat['z_true']=cat['z']
pz=getPhotoZ(o)
print "Applying photoz..."
cat.setPhotoZ(pz,apply_to_data=True)
## now get the histogram
print "Getting cumulatives..."
culs=cat.photoz.cPofZ(cat.data, cat['z_true'])
print "Plotting"
n, bins, patches = pylab.hist(culs, bins=min(N/100,100))
pylab.xlabel("cumulative probability")
pylab.ylabel("frequency")
pylab.show()

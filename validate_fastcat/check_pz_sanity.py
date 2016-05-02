#!/usr/bin/env python
import sys
import fastcat as fc
from optparse import OptionParser
import pylab

parser = OptionParser(usage="""

validate_fastcat/check_pz_sanity.py HDF5_fastcat_file

This code will check the sanity of PZ in a given file by computing the
cumulative PDF for PZ at the true redshift and histograming it.

Hence you need to have true redshift available (option --ztrue in mkcat.py)
""")

(o, args) = parser.parse_args()
if len(args) != 1:
    parser.error("incorrect number of arguments")
    sys.exit(1)

print "Reading..."
cat=fc.Catalog(read_from=args[0])
if "z_true" not in cat.data.dtype.names:
    print "Need z_true column."
    sys.exit(1)
print "Calculating cumulative prob..."
culs=cat.photoz.cPofZ(cat.data, cat['z_true'])

print "Plotting"
n, bins, patches = pylab.hist(culs, bins=1000)
pylab.xlabel("cumulative probability")
pylab.ylabel("frequency")
pylab.show()

#!/usr/bin/env python
#
# print the command line used to get catalogs
#
import h5py
import sys
for fname in sys.argv[1:]:
    print "#",fname
    try:
        print h5py.File(fname)['meta'].attrs['command_line']
    except:
        print 'echo "Fail %s "'%fname


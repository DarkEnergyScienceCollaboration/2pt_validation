#!/usr/bin/env python
import sys
from subs import *
from time import time


def main():
    t_0 = time()
    o,args = setupOptions()
    initMPI(o)
    process_catalog(o)
    print 'Total time (seconds): ', time()-t_0

if __name__=="__main__":
    main()

#TODO: add options for
# - z-binning
# - mask apodization
#TODO: extract ell-windows and write in SACC format
#TODO: subtract shot noise
#TODO: deal with recursive debiasing
#TODO: change syntax for templates_fname in subs.py to accept None

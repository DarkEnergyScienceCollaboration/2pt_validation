#!/usr/bin/env python
from subs import *
import sys

if len(sys.argv) != 8 :
    print "Input parameters : fname_in nside fname_bins_z theta_apo templates_fname delta_ell fname_out"
    exit(1)

fname_in=sys.argv[1] #
nside=int(sys.argv[2]) #
fname_bins_z=sys.argv[3] #
theta_apo=float(sys.argv[4]) #
templates_fname=sys.argv[5] #
delta_ell=int(sys.argv[6]) #
fname_out=sys.argv[7] #

process_catalog(fname_in,fname_bins_z,nside,fname_out,
                apodization_scale=theta_apo,fname_templates=templates_fname,bins_ell=delta_ell)

#TODO: add options for
# - z-binning
# - mask apodization
#TODO: extract ell-windows and write in SACC format
#TODO: subtract shot noise
#TODO: deal with recursive debiasing

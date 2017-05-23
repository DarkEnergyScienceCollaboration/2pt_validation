#!/usr/bin/env python
from subs import *
import sys
from optparse import OptionParser

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

(o, args) = parser.parse_args()

process_catalog(o.fname_in,o.fname_bins_z,o.nside,o.fname_out,
                apodization_scale=o.theta_apo,fname_templates=o.templates_fname,bins_ell=o.delta_ell)

#TODO: add options for
# - z-binning
# - mask apodization
#TODO: extract ell-windows and write in SACC format
#TODO: subtract shot noise
#TODO: deal with recursive debiasing
#TODO: change syntax for templates_fname in subs.py to accept None

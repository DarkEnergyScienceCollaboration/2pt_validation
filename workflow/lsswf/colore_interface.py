#
#

from parsl import App, DataFlowKernel, ThreadPoolExecutor
from parsl.dataflow.futures import Future
from parsl_interface import dfk
from HashedFileArray import HashedFileArray
from DirectorySetup import rootDir, maybeMakeDir
from colore_defaults import Pk_default, Nz_default, Bz_default
import hashlib
import os

class colore(dict):
    def __init__(self,params=[]):
        ## we are dictionary with configuration
        self.rootDir=maybeMakeDir("colore")
        self.commonDir=maybeMakeDir("colore/common")
        self.runs=maybeMakeDir("colore/runs")
        self.update(self.defaultConfig())
        ## syntax candy
        for name,val in params:
            sec,key=name.split('/')
            if key in self[sec]:
                self[sec][key]=val
            else:
                raise NameError(name)
        self.image="docker:slosar/desc_lss:v0.1"
        self.walltime='00:10:00'
        self.Nnodes=1
        self.partition='debug'
        self.NMPI=4
        
    def defaultConfig(self):
        return {
            "global" : {
                "prefix_out" : "this/out",
                "output_format" : "HDF5",
                "output_density" : True,
                "pk_filename" : HashedFileArray(Pk_default,self.commonDir),
                "z_min" : 0.001,
                "z_max" : 0.450,
                "seed" : 1003,
                "write_pred" : False,
                "pred_dz" : 0.1
            },
            "field_par": {
                "r_smooth" : 5.,
                "smooth_potential" : True,
                "n_grid" : 256,
                "dens_type" : 1,
                "lpt_buffer_fraction" : 0.6,
                "lpt_interp_type" : 1,
                "output_lpt" : 0
            },
            "cosmo_par": { 
                "omega_M" : 0.3,
                "omega_L" : 0.7,
                "omega_B" : 0.05,
                "h" : 0.7,
                "w" : -1.0,
                "ns" : 0.96,
                "sigma_8" : 0.803869
            },
            'srcs0':
            {
                "nz_filename" : HashedFileArray(Nz_default,self.commonDir),
                "bias_filename" : HashedFileArray(Bz_default,self.commonDir),
                "include_shear" : False
            }
        }                    

    def _generate_cfg(self):
        cfg=""
        for sec in self.keys():
            cfg+="{}:\n".format(sec)
            cfg+="{\n"
            for key,val in self[sec].items():
                if type(val) in [str,float,int]:
                    cfg+="{}= {}\n".format(key,val)
                elif type(val)==bool:
                    cfg+="{}= ".format(key)
                    if val:
                        cfg+="true\n"
                    else:
                        cfg+="false\n"
                elif isinstance(val,HashedFileArray):
                    val.genFile()
                    cfg+="{}= common/{}\n".format(key,val.filename())
                else:
                    raise NotImplemented()
            cfg+="}\n\n"
        return cfg.encode('ascii')

    ###@App('bash', dfk)
    def run(self):
        cfg=self._generate_cfg()
        m = hashlib.md5()
        m.update(cfg)
        name=m.hexdigest()
        thisdir=maybeMakeDir("colore/runs/"+name)
        open(thisdir+"/params.cfg",'wb').write(cfg)
        ## TODO
        # fiddle with DFK to set number of nodes, walltime,etc
        ##
        if os.path.exists(thisdir+"/out_srcs_0.h5"):
            return "true"
        else:
            return 'srun -n {} shifter  '.format(self.NMPI) +\
                    ' --volume="{}:/common;{}:/this" '.format(self.commonDir,thisdir) +\
                    ' /home/lss/CoLoRe/runCoLoRe /this/params.cfg >/this/output.log 2>/this/output.err'


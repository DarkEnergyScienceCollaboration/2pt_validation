#
# PARSL interface for LSS workflow
#

from parsl import App, DataFlowKernel, ThreadPoolExecutor
from parsl.dataflow.futures import Future
import os

def get_dfk():
    if 'NERSC_HOST' in os.environ:
        if os.environ['NERSC_HOST']=='cori':
            return dfk_cori()
        else:
            raise NotImplemented()
    else:
        raise NotImplemented()

def dfk_cori():
    ## temporary to test
    workers = ThreadPoolExecutor(max_workers=4)
    dfk = DataFlowKernel(executors=[workers])
    return dfk
    ## throws

    ## parsl.dataflow.error.ControllerErr: Controller init
    ## failed:Reason:IPPController failed to start: [Errno 2] No such file or directory: 'ipcontroller'
    return DataFlowKernel( config = {
        "sites" : [
            { "site" : "Cori.Remote.IPP",
              "auth" : {
                  "channel" : "ssh",
                  "hostname" : "cori.nersc.gov",
                  "username" : "yadunand",
                  "scriptDir" : "/global/homes/y/yadunand/parsl_scripts"
              },
          "execution" : {
              "executor" : "ipp",
              "provider" : "slurm",  # LIKELY SHOULD BE BOUND TO SITE
              "script_dir" : ".scripts",
              "block" : { # Definition of a block
                  "nodes" : 2,            # of nodes in that block
                  "taskBlocks" : 1,       # total tasks in a block
                  "walltime" : "00:10:00",
                  "initBlocks" : 1,
                  "minBlocks" : 0,
                  "maxBlocks" : 1,
                  "scriptDir" : ".",
                  "options" : {
                      "partition" : "debug",
                      "overrides" : '''
#SBATCH --constraint=haswell
module load python/3.5-anaconda ;
source activate /global/homes/y/yadunand/.conda/envs/parsl_env_3.5
'''
                  }
              }
          }
        }
        ],
    "globals" : {   "lazyErrors" : True },
        "controller" : { "publicIp" : '*' }
    })



### global 

dfk = get_dfk()

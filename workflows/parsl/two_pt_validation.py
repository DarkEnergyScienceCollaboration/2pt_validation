import parsl
from parsl_setup import dfk

__all__ = ['colore', 'fast_colore', 'namaster_fastcat', 'mk_theory']

def write_colore_param_file(index, output_path, pk_filename, nz_filename,
                            bias_filename):
    with open('colore_param_file_template.txt', 'r') as input_:
        param_file = 'param_colore_{}.cfg'.format(index)
        with open(param_file, 'w') as output:
            content = ''.join(input_.readlines())
            output.write(content.format(index=index, output_path=output_path,
                                        pk_filename=pk_filename,
                                        nz_filename=nz_filename,
                                        bias_filename=bias_filename))
    return param_file

@parsl.App('bash', dfk)
def colore(param_file, outputs=[], stdout=None, stderr=None):
    return '''CoLoRe {0}'''

@parsl.App('bash', dfk)
def fast_colore(param_file, out_path, pz_type='gauss', pz_sigma=0.05,
                oextra='20percent'):
    return '''module load python/2.7-anaconda
module load h5py-parallel
srun -n 32 python mkcat.py --params_file={0} --opath={1} --pz_type={pz_type} --pz_sigma={pz_sigma} --mpi --oextra={oextra}
'''

@parsl.App('bash', dfk)
def namaster_fastcat(fastcat_catalog, out_file, nz_bins_file,
                     templates='none', delta_ell=25, nside=2048):
    return '''srun -n 32 python namaster_interface.py --input-file {0} --output-file {1} --nz-bins-file {2} --templates {templates} --delta-ell {delta_ell} --nside {nside} --mpi
'''

@parsl.App('bash', dfk)
def mk_theory(input_sacc, output_sacc, param_file, bias_file, shot_noise):
    return '''mk_theory.py --input-sacc {0} -output-sacc {1} --param-file {2} --bias-file {3} --shot-noise-file {4}
'''

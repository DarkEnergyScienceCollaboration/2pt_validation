import os
from two_pt_validation import *

output_path = os.environ['OUTPUT_PATH']
pk_filename = os.environ['PK_FILENAME']
nz_filename = os.environ['NZ_FILENAME']
bias_filename = os.environ['BIAS_FILENAME']

colore_results = []
for i in range(100):
    param_file = write_colore_param_file(i, output_path, pk_filename,
                                         nz_filename, bias_filename)
    colore_results.append(colore(param_file))

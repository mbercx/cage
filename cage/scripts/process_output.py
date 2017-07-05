# Encoding: utf-8

import sys
import os

import pymatgen.io.nwchem as nwchem

output_file = sys.argv[1]
dir_name = os.path.dirname(output_file)

try:
    output = nwchem.NwOutput(output_file)
except:
    raise IOError('Could not find proper nwchem output file.')

output.to_file(os.path.join(dir_name, 'data.json'))
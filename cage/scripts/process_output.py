# Encoding: utf-8

import sys

import pymatgen.io.nwchem as nwchem

output_file = sys.argv[1]

try:
    output = nwchem.NwOutput(output_file)
except:
    raise IOError('Could not find proper nwchem output file.')

# TODO Finish this script once the NwOutput is MSONable
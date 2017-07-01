# Encoding UTF-8

import sys
import os

from pymatgen.io import nwchem

"""
Script that checks if a calculation has completed successfully from the ouput
file.
"""
# TODO Add method of extracting data more quickly

filename = sys.argv[1]

try:
    out = nwchem.NwOutput(filename)
except:
    raise IOError('Could not find proper nwchem output file.')

error = False
for data in out.data:
    if data['has_error']:
        error = True

print('File: ' + os.path.abspath(filename))
print('Calculation has error = ' + str(error))
# Encoding UTF-8

import sys
import os

from pymatgen.io import nwchem
from json import JSONDecodeError

"""
Script that checks if a calculation has completed successfully from the ouput
file.
"""
# TODO Add method of extracting data more quickly

filename = sys.argv[1]

try:
    out = nwchem.NwOutput(filename, fmt='json')
except JSONDecodeError:
    try:
        out = nwchem.NwOutput(filename)
    except:
        raise IOError('File not found.')

try:
    error = False
    for data in out.data:
        if data['has_error']:
            error = True

    print('File: ' + os.path.abspath(filename))
    if out.data[-1]['task_time'] != 0:
        print('Calculation completed in ' + str(out.data[-1]['task_time']) + 's')
    else:
        print('No timing information found. Calculation might not have '
              'completed successfully.')
    print('Calculation has error: ' + str(error))

except NameError:
    print("No data found in file!")
# Encoding utf-8
# Written for Python 3.6

import sys
import pymatgen.io.nwchem as nwchem

"""
Quick script to study the total energy change of the system during an
optimization.
"""
if sys.argv[1]:
    filename = sys.argv[1]
    print('\nAnalyzing energies in ' + filename + '...\n')
else:
    OSError('No target output file provided.')

out = nwchem.NwOutput(filename)
energy_data = out.data[0]['energies']

for i in range(len(energy_data)):
    print('Step ' + str(i + 1) + ': Energy = ' + str(energy_data[i]))
# Encoding utf-8

import sys
import pymatgen.io.nwchem as nwchem

"""
Quick script to study the total energy change of the system during an
optimization.
"""

filename = sys.argv[1]
print('\nAnalyzing energies in ' + filename + '...\n')

out = nwchem.NwOutput(filename)
energy_data = out.data[0]['energies']

for i in range(len(energy_data)):
    print('Step ' + str(i + 1) + ': Energy = ' + str(energy_data[i]))
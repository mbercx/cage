# Encoding: utf-8

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pymatgen as pmg
import pymatgen.io.nwchem as nwchem

import cage.utils as utils
from cage.core import Cage

"""
Script to extract all the necessary information from a facet energy landscape
analysis.
"""

OUTPUT_FILENAME = 'result.out'

dirname = os.path.abspath(sys.argv[1])
if not os.path.isdir(dirname):
    raise IOError("Directory " + sys.argv[1] + " not found.")

# Find all the directories in the specified directory
dirlist = [d for d in os.listdir(dirname)
           if os.path.isdir(os.path.join(dirname, d))]

# Get all of the output of the calculations in the directory
output = []
for directory in dirlist:
    try:
        data = nwchem.NwOutput(os.path.join(dirname, directory,
                                            OUTPUT_FILENAME)).data[-1]
    except FileNotFoundError:
        print('Could not find output in directory ' + directory)

    # Check if the output has an error
    if not data['has_error']:
        output.append(data)
        print('Output found in directory ' + directory)
    else:
        print('Error found in output in directory ' + directory)
        output.append(data)

# Note that now only the output is left which did not contain any errors

# TODO Find a better way: Store the initial cage in a .json This is important in case the first directory actually contains a calculation that was restarted

# Find the initial facet
cage_init = Cage.from_molecule(output[0]['molecules'][0])
facet_init = cage_init.facets[0]
for facet in cage_init.facets:
    if utils.distance(facet.center, cage_init.cart_coords[-1]) < \
            utils.distance(facet_init.center, cage_init.cart_coords[-1]):
        facet_init = facet

# Get the distance to the facet, the angle between the facet center and the Li
# position, and the energy for each result
results = []
for data in output:

    cage_final = Cage.from_molecule(data['molecules'][-1])
    facet_final = cage_final.facets[0]
    for facet in cage_final.facets:
        if utils.distance(facet.center, cage_final.cart_coords[-1]) < \
                utils.distance(facet_final.center, cage_final.cart_coords[-1]):
            facet_final = facet

    init_energy = data['energies'][0]
    final_energy = data['energies'][-1]
    final_Li_coord = cage_final.cart_coords[-1]

    dist_Li_init = utils.distance(facet_init.center, final_Li_coord)
    dist_Li_final = utils.distance(facet_final.center, final_Li_coord)

    if pmg.Element('C') in cage_init.species:

        angle_Li_init = utils.angle_between(facet_init.center, final_Li_coord)
        angle_Li_final =utils.angle_between(facet_final.center, final_Li_coord)

        results.append([dist_Li_init, angle_Li_init, final_energy])
    else:
        results.append([dist_Li_init, dist_Li_final, init_energy, final_energy])

if pmg.Element('C') in cage_init.species:

    pass

    # Find out how many theta there are
    # radius = np.array(results)
    # n_theta = 1
    # r1 = results[0][0]
    # rn = results[]
    # while

else:
    results = np.array(sorted(results, key=lambda results: results[0]))

    plt.figure(1)
    plt.subplot(211)
    plt.ylabel('Total Energy')
    plt.xlabel('Distance to initial facet')
    plt.plot(results[:,0],results[:,2])

    plt.subplot(212)
    plt.ylabel('Total Energy')
    plt.xlabel('Distance to relaxed facet')
    plt.plot(results[:,1],results[:,2])
    plt.show()


# encoding: utf-8
# Written for Python 3.6

import sys
import cage
import numpy as np
import pymatgen as pmg
import pymatgen.io.nwchem as nwchem
import os
import shutil

"""
Script to set up the calculations for all the non-equivalent facets of a cage
molecule, defined by a VASP POSCAR file provided as an argument.

This script was written quite quickly, so the code is downright dirty.

:argument filename: The filename of the VASP 
"""

filename = sys.argv[1]

# Parameters
LINE = [1, 4] # Distance from the center of the facet for the Lithium

DENSITY = 5

BASIS = {'*': "aug-pcseg-1"}

ALT_SETUP = {'dft': {'iterations': '300',
                 'xc': 'xpbe96 xpbe96',
                 'direct': '',
                 'convergence': 'energy 1e-6'}}

GEO_SETUP = {'nocenter', "units angstroms"}

# Load the POSCAR into a Cagel
mol = cage.facetsym.Cage.from_poscar(filename)

# Find the non-equivalent facets
facets = mol.find_noneq_facets()

# For each facet, set up the calculation input files
facetnumber = 1
for neq_facet in facets:

    # Define the line on which to place the Lithium
    center = neq_facet.center
    normal = neq_facet.normal / np.linalg.norm(neq_facet.normal)
    endpoints = [center + LINE[0]*normal, center + LINE[1]*normal]
    line = cage.landscape.Landscape(endpoints, DENSITY)

    # Create the list of structures with lithium at varying distances to
    # the facet.
    li = pmg.Specie('Li', 0)
    structures = []
    for point in line.points:
        configuration = mol.copy()
        configuration.append(li, point)
        structures.append(configuration)

    # Find the facet that is the farthest away from the facet being studied
    distance = 0
    for facet in mol.facets:
        newdistance = np.linalg.norm(neq_facet.center - facet.center)
        if newdistance > distance:
            far_facet = facet
            distance = newdistance

    # Find the corresponding atom numbers in string format
    site_numbers = []
    for i in range(len(mol.sites)):
        if mol.sites[i] in far_facet.sites:
            site_numbers.append(str(i+1) + ' ')

    site_numbers.append(str(len(mol.sites) + 1)) # Add the lithium site
    site_numbers = ''.join(site_numbers)

    # Calculate the constraints on the atoms, and add them to the setup
    constraints = {'fix atom': site_numbers}
    ALT_SETUP['constraints'] = constraints

    # Set up the task for the calculations
    tasks = [nwchem.NwTask(-1, None, BASIS, theory='dft', operation='optimize',
                           alternate_directives=ALT_SETUP)]

    # Set up the input files, and place the geometry files in a subdirectory
    # of the composition directory
    study = cage.study.Study(structures, tasks)
    study.set_up_input('.', sort_comp=False, geometry_options=GEO_SETUP)

    facet_dir = 'facet' + str(facetnumber)
    try:
        os.mkdir(facet_dir)
    except FileExistsError:
        pass

    for i in range(len(structures)):
        shutil.move(os.path.join('geo' + str(i + 1)),
                    os.path.join(facet_dir, 'geo' + str(i + 1)))

    facetnumber += 1

    # Write out a facet json file
    neq_facet.to(fmt='json',filename=os.path.join(facet_dir,'facet.json'))

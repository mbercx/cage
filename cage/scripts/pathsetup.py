# Encoding utf-8

import sys
import os
import io

import pymatgen as pmg
import numpy as np
import pymatgen.io.nwchem as nw

import pymatgen.io.vasp as vasp

from cage.facetsym import Cage
from cage.facetsym import Facet
from pymatgen.analysis import path_finder


"""
Script to set up the directories and input files to calculate the minimum
energy paths for all non-equivalent facet combinations in a Cage molecule.
"""

# Lithium parameters
START_DIST = 2

# Calculation parameters
# BASIS = {'*': "aug-pcseg-1"}
BASIS = {'*': "aug-cc-pVDZ"}

THEORY_SETUP = {'iterations': '300',
                'xc': 'xpbe96 xpbe96',
                'direct': '',
                'smear': '0.01',
                'convergence energy': '1e-4',
                'convergence density': '1e-2',
                'convergence gradient': '1e-2',
                'convergence damp': '70'}

GEO_SETUP = {"noautosym", "noautoz", 'nocenter', "units angstroms"}

ALT_SETUP = {"driver": {'loose': '', 'maxiter': '100'}}

OPERATION = 'optimize'

def main():

    # Read the Molecule from the input file
    filename = sys.argv[1]
    mol = Cage.from_poscar(filename)
    mol.find_surface_facets()

    # TODO Same problem as in facetsetup.py. Find a way to calculate charge automatically
    if pmg.Element('C') in [site.specie for site in mol.sites]:
        mol.set_charge_and_spin(charge=0)
    else:
        mol.set_charge_and_spin(charge=-1)

    # Find the paths, i.e. the List of facet combinations
    paths = mol.find_facet_paths()

    tasks = [nw.NwTask(mol.charge, None, BASIS, theory='dft',
                       operation=OPERATION,
                       theory_directives=THEORY_SETUP,
                       alternate_directives=ALT_SETUP),]

    path_number = 1
    for path in paths:

        # Set up the start and end molecules
        start_molecule = mol.copy()
        end_molecule = mol.copy()

        start_molecule.append(pmg.Specie('Li', 1), path[0].center +
                              START_DIST * path[0].normal)
        end_molecule.append(pmg.Specie('Li', 1), path[1].center +
                              START_DIST * path[1].normal)

        # Make the path directory
        path_dir = 'path' + str(path_number)
        try:
            os.mkdir(path_dir)
        except FileExistsError:
            pass

        # Write out the start and end molecule
        start_molecule.to(fmt='xyz',
                          filename=os.path.join(path_dir, 'start.xyz'))
        end_molecule.to(fmt='xyz',
                        filename=os.path.join(path_dir, 'end.xyz'))

        # Make the directories for the start and end molecules optimizations
        try:
            os.mkdir(os.path.join(path_dir, 'start'))
        except FileExistsError:
            pass
        try:
            os.mkdir(os.path.join(path_dir, 'end'))
        except FileExistsError:
            pass

        # Set up the input
        nw.NwInput(start_molecule, tasks, geometry_options=GEO_SETUP)\
            .write_file(os.path.join(path_dir, 'start', 'input'))
        nw.NwInput(end_molecule, tasks, geometry_options=GEO_SETUP) \
            .write_file(os.path.join(path_dir, 'end', 'input'))

        path_number += 1

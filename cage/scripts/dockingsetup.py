# encoding: utf-8
# Written for Python 3.6

import os
import sys

import numpy as np
import pymatgen as pmg
import pymatgen.io.nwchem as nwchem

import cage

"""
Script to set up the calculations for all the non-equivalent facets of a cage
molecule, defined by a VASP POSCAR file provided as an argument.

This script is written specifically to study the facets of CxByHz molecules.
In case you want to use it for anything else, it will probably need some
adjustments.

:argument filename: The filename of the VASP 
"""
# TODO Make these parameters defaults, but allow the user to change them with arguments in the CLI

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

# Input check
try:
    filename = sys.argv[1]
except IndexError:
    raise IOError('No POSCAR file provided.')


def main():

    # Load the POSCAR into a Cage
    molecule = cage.facetsym.Cage.from_poscar(filename)
    molecule.find_surface_facets()

    # Find the non-equivalent facets
    facets = molecule.find_noneq_facets()

    # For each docking point, set up the calculation input file
    dock_number = 1
    for neq_facet in facets:

        # Set up the initial lithium site
        mol = molecule.copy()
        mol.find_surface_facets()
        mol.append(pmg.Specie('Li', 1), neq_facet.center +
                   START_DIST*neq_facet.normal)

        # Set the charge for the molecule
        if pmg.Element('C') in [site.specie for site in mol.sites]:
            mol.set_charge_and_spin(charge=0)
        else:
            mol.set_charge_and_spin(charge=-1)

        # Add the constraints
        ALT_SETUP['constraints'] = find_constraints(mol, neq_facet)

        # Set up the task for the calculations
        tasks = [nwchem.NwTask(mol.charge, None, BASIS,
                               theory='dft',
                               operation=OPERATION,
                               theory_directives=THEORY_SETUP,
                               alternate_directives=ALT_SETUP)]

        dock_dir = 'dock' + str(dock_number)

        try:
            os.mkdir(dock_dir)
        except FileExistsError:
            pass

        # Set up input
        input = nwchem.NwInput(mol, tasks, geometry_options=GEO_SETUP)
        input.write_file(os.path.join(dock_dir, 'input'))

        # Write out a facet json file
        neq_facet.to(fmt='json', filename=os.path.join(dock_dir, 'facet.json'))

        dock_number += 1

###########
# METHODS #
###########


def find_constraints(mol, neq_facet):
    """
    Find the necessary constraints for the molecules, depending on the facet
    this might be fixing the whole opposite facet, or only one vertex.
    :return:
    """
    # Find the facet that is the farthest away from the facet being studied
    distance = 0
    far_facet = None
    for facet in mol.facets:
        newdistance = np.linalg.norm(neq_facet.center - facet.center)
        if newdistance > distance:
            far_facet = facet
            distance = newdistance

    # Find the corresponding atom numbers in string format
    site_numbers = []
    for i in range(len(mol.sites)):
        if mol.sites[i] in far_facet.sites:
            site_numbers.append(str(i + 1) + ' ')

    site_numbers = ''.join(site_numbers)

    # Set up the constraints on the atoms
    return {'fix atom': site_numbers}

if __name__ == '__main__':
    main()

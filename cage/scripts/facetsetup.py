# encoding: utf-8
# Written for Python 3.6

import os
import shutil
import sys

import numpy as np
import pymatgen as pmg
import pymatgen.io.nwchem as nwchem

import math
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

# !TODO Currently not functional, due to API changes

# Landscape parameters
LINE = [1, 6] # Distance from the center of the facet for the Lithium
DENSITY = 5

# Calculation parameters
#BASIS = {'*': "aug-pcseg-1"}
BASIS = {'*': "aug-cc-pVDZ"}

THEORY_SETUP = {'iterations': '300',
                'xc': 'xpbe96 xpbe96',
                'direct': '',
                'smear':'0.01',
                'convergence energy': '1e-4',
                'convergence density': '1e-2',
                'convergence gradient': '1e-2',
                'convergence damp':'70'}

GEO_SETUP = {"noautosym", "noautoz", 'nocenter', "units angstroms"}

ALT_SETUP = {}

DRIVER_SETUP = {'loose':'', 'maxiter': '100'}

OPERATION = 'energy'

# Input check
try:
    # Take the filename argument from the user
    filename = sys.argv[2]
    # Take the operation input
    OPERATION = sys.argv[1]
except IndexError:
    # Take the filename argument from the user
    try:
        filename = sys.argv[1]
    except IndexError:
        raise IOError('No POSCAR file provided.')

def main():

    # Load the POSCAR into a Cage
    mol = cage.core.Cage.from_poscar(filename)

    # Find the non-equivalent facets
    facets = mol.find_noneq_facets()

    # For each facet, set up the calculation input files
    facetnumber = 1
    for neq_facet in facets:

        # Find the molecules for each lithium distance
        molecules = set_up_molecules(mol, neq_facet)

        if OPERATION == 'optimize':
            # Add the constraints
            ALT_SETUP['constraints'] = find_constraints(mol,neq_facet)
            ALT_SETUP['driver'] = DRIVER_SETUP

        # Set up the task for the calculations
        tasks = [nwchem.NwTask(molecules[0].charge, None, BASIS,
                               theory='dft',
                               operation=OPERATION,
                               theory_directives=THEORY_SETUP,
                               alternate_directives=ALT_SETUP)]

        # Set up the input files, and place the geometry files in a subdirectory
        # of the facet directory
        study = cage.study.Study(molecules, tasks)
        study.set_up_input('.', sort_comp=False, geometry_options=GEO_SETUP)

        facet_dir = 'facet' + str(facetnumber)
        try:
            os.mkdir(facet_dir)
        except FileExistsError:
            pass

        for i in range(len(molecules)):
            shutil.move(os.path.join('geo' + str(i + 1)), facet_dir)

        facetnumber += 1

        # Write out a facet json file
        neq_facet.to(fmt='json',filename=os.path.join(facet_dir,'facet.json'))

###########
# METHODS #
###########

def set_up_molecules(mol, facet):
    """
    Set up the molecules from the lithium Landscape for a facet.
    :return:
    """
    # Define the line on which to place the Lithium
    center = facet.center
    normal = facet.normal / np.linalg.norm(facet.normal)
    endpoints = [center + LINE[0] * normal, center + LINE[1] * normal]
    line = cage.landscape.Landscape.from_vertices(endpoints, DENSITY)

    # If Carbon is in the molecule, expand the energy landscape
    # TODO: Make this less specific to carbon. This whole module will probably need an overhaul for that.
    if pmg.Element('C') in [site.specie for site in mol.sites]:
        # Find the rotation axis
        C_site = [site for site in mol.sites
                  if site.specie == pmg.Element('C')][0]
        connection_vector = C_site.coords - facet.center
        axis = np.cross(facet.normal, connection_vector)
        axis = axis/np.linalg.norm(axis)

        # Overestimate rotation angle and remove landscape points that are
        # closer to another facet.
        angle = math.pi/8
        line.extend_by_rotation(angle*axis)

    # Create the list of molecules with lithium at varying distances to
    # the facet.
    li = pmg.Specie('Li', +1)
    molecules = []
    for point in line.points:
        configuration = mol.copy()
        # If carbon is in the molecule, the charge is -1, -2 otherwise
        # TODO Find a good way to calculate the charge, if possible
        if pmg.Element('C') in [site.specie for site in configuration.sites]:
            configuration.set_charge_and_spin(charge=0)
        else:
            configuration.set_charge_and_spin(charge=-1)

        try:
            configuration.append(li, point)
            molecules.append(configuration)
        except ValueError:
            print('ValueError detected when appending the Li site. '
                  'Ignoring this point in the energy landscape.')

    return molecules

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

    # If there is carbon in the molecule, only fix one atom
    # if pmg.Element('C') in [site.specie for site in mol.sites]:
    #    site_numbers = site_numbers[0:1]
    # Left here for future reference: This does not work.

    site_numbers.append(str(len(mol.sites) + 1))  # Add the lithium site
    site_numbers = ''.join(site_numbers)

    # Calculate the constraints on the atoms
    return {'fix atom': site_numbers}

if __name__ == '__main__':
    main()
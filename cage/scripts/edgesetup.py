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
Script to set up the calculations for all the non-equivalent edge sharing paths
in a cage molecule. The script sets up the energy landscape calculations for
each non-equivalent edge sharing path.

This script is written specifically to study the facets of CxByHz molecules.
In case you want to use it for anything else, it will probably need some
adjustments.

"""
# TODO Make these parameters defaults, but allow the user to change them with arguments in the CLI

# TODO The current set up does not work very well for molucules that are not spherically shaped

# Landscape parameters
LINE = [2, 5] # Distance from the center of the molecule for the Lithium
DENSITY_R = 10 # Density of Lithium along the r coordinate
DENSITY_theta = 60 # Density of Lithium along the theta coordinate

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
    mol = cage.facetsym.Cage.from_poscar(filename)
    mol.find_surface_facets()

    # Find the non-equivalent paths
    paths = mol.find_facet_paths()

    # Find the edge paths
    edge_paths = []
    for path in paths:
        if len(set(path[0].sites) & set(path[1].sites)) == 2:
            edge_paths.append(path)

    total_mol = mol.copy()

    # For each facet, set up the calculation input files
    edge_number = 1

    for path in edge_paths[:2]+edge_paths[3:4]:#+edge_paths[2:5]:

        edge_dir = 'edge' + str(edge_number)
        try:
            os.mkdir(edge_dir)
        except FileExistsError:
            pass

        mol.to(fmt='json', filename=os.path.join(edge_dir, 'mol.json'))
        path[0].to(fmt='json', filename=os.path.join(edge_dir,
                                                     'init_facet.json'))
        path[1].to(fmt='json', filename=os.path.join(edge_dir,
                                                     'final_facet.json'))

        edge_mol = mol.copy()
        facet1 = path[0].copy()
        facet2 = path[1].copy()

        intersection = facet1.get_normal_intersection(facet2)

        # Redefine the origin for the molecule and path
        # total_mol.redefine_origin(intersection)
        # edge_mol.redefine_origin(intersection)
        # facet1.redefine_origin(intersection)
        # facet2.redefine_origin(intersection)

        landscape = set_up_landscape(facet1, facet2)

        molecules = set_up_molecules(edge_mol, landscape)

        # Set up a molecule to visualize the edge
        for point in landscape.points:
            try:
                total_mol.append(pmg.Specie('Li', 1), point,
                                validate_proximity=False)
                edge_mol.append(pmg.Specie('Li', 1), point,
                                validate_proximity=False)
            except ValueError:
                pass

        # total_mol.redefine_origin(-intersection)

        edge_file = 'edge' + str(edge_number) + '.xyz'
        edge_mol.to(fmt='xyz', filename=os.path.join(edge_dir, edge_file))

        if OPERATION == 'optimize':
            # Add the constraints
            ALT_SETUP['constraints'] = find_constraints(mol,facet1)
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
        study.set_up_input(edge_dir, sort_comp=False,
                           geometry_options=GEO_SETUP)
        edge_number += 1

    total_mol.to(fmt='xyz', filename='total_mol.xyz')

###########
# METHODS #
###########

def set_up_landscape(facet1, facet2):
    """

    :param facet1:
    :param facet2:
    :return:
    """
    line_vector = facet1.center/np.linalg.norm(facet1.center)
    lands = cage.landscape.Landscape.from_vertices(
        [line_vector * LINE[0], line_vector * LINE[1]]
    )
    axis = np.cross(facet1.normal, facet2.normal)
    angle = math.asin(np.linalg.norm(axis))
    axis = axis * angle / np.linalg.norm(axis)

    lands.extend_by_rotation(axis, DENSITY_theta)

    return lands

def set_up_molecules(mol, landscape):
    """
    Set up the molecules from the lithium Landscape for a facet.
    :return:
    """

    # Create the list of molecules with lithium at varying distances to
    # the facet.
    li = pmg.Specie('Li', +1)
    molecules = []
    for point in landscape.points:
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

def find_constraints(mol, facet):
    """
    Find the necessary constraints for the molecules, depending on the facet
    this might be fixing the whole opposite facet, or only one vertex.
    :return:
    """
    # Find the facet that is the farthest away from the facet being studied
    distance = 0
    far_facet = None
    for other_facet in mol.facets:
        newdistance = np.linalg.norm(facet.center - other_facet.center)
        if newdistance > distance:
            far_facet = other_facet
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
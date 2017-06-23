# Encoding utf-8

import pymatgen as pmg
import numpy as np
import sys

from cage.facetsym import Cage
from pymatgen.analysis import path_finder


"""
Script to set up the directories and input files to calculate the minimum
energy paths for all non-equivalent facet combinations in a Cage molecule.
"""

# Parameters

filename = sys.argv[1]

## Find all the non-equivalent facet combinations, i.e. path endpoints

mol = Cage.from_poscar(filename)
paths = mol.find_facet_paths()

# Find the docking points of the facets in the path

#TODO Load the optimal docking points, calculated using NwChem

## Set up the potential field of the molecule

# Find the distance between the origin and the furthest atom in the molecule

max_distance = 0
for site in mol.sites:
    distance = np.linalg.norm(site.coords)
    if  distance > max_distance:
        max_distance = distance

# Double that, and use it to set up the lattice

lat_const = 2*max_distance
lattice = pmg.Lattice(np.array([[lat_const, 0, 0], [0, lat_const, 0],
                                [0, 0, lat_const]]))

# Use the lattice in combination with the sites of the molecule to set up a
# Structure, for which we can calculate the potential field

struc = pmg.Structure(lattice, mol.species, mol.cart_coords,
                      coords_are_cartesian=True)

# Use the pymatgen.analysis.path_finder FreeVolumePotential to find the
# potential field of the molecule.

potential = path_finder.FreeVolumePotential(struc, [100, 100, 100])

## Find the initial set up for the NEB calculations.

# Set up the beginning and end sites of the various facets
struc.copy()
struc.append('Li', coords, coords_are_cartesian=True)


# Use the NEBPathfinder class in combination with the potential field to
# calculate the initial pathways

# Find a way to convert the pathways into the xyz format that NwChem wants

## Set up the directories and calculation input files
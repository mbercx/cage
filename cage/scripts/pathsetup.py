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

# Parameters

NIMAGES = 10  # Standard number of images along every path

# Calculation parameters
BASIS = {'*': "aug-pcseg-1"}

THEORY_SETUP = {'iterations': '300',
                'xc': 'xpbe96 xpbe96',
                'direct': '',
                'convergence energy': '1e-4',
                'convergence density': '1e-2',
                'convergence gradient': '1e-2',
                'convergence damp':'70'}

GEO_SETUP = {'nocenter', "units angstroms"}

ALT_SETUP = {'neb': {'nbeads':str(NIMAGES),
                     'xyz_path':'path.neb',
                     'loose':'',
                     'maxiter': '100'}}


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

    # Set up the structure for the initial pathfinder
    struc = set_up_structure(mol)

    # Find the paths, i.e. the List of facet combinations
    paths = mol.find_facet_paths()

    # Set up the docking points of the facets
    dock_points = extract_docking_points('../docking')

    tasks = [nw.NwTask(mol.charge, None, BASIS, theory='dft', operation='neb',
                       theory_directives=THEORY_SETUP,
                       alternate_directives=ALT_SETUP),]

    path_number = 1
    for path in paths:

        # Find the docking points for the path
        endpoints = find_path_endpoints(mol, path, dock_points)

        # Set up the start and end structures
        start_struct = struc.copy()
        start_struct.append('Li', endpoints[0], coords_are_cartesian=True)

        end_struct = struc.copy()
        end_struct.append('Li', endpoints[1], coords_are_cartesian=True)

        # Interpolate to find the images
        structures = start_struct.interpolate(end_struct, nimages=NIMAGES)

        molecules = [pmg.Molecule(struct.species, struct.cart_coords)
                     for struct in structures]

        # Move the lithium positions on an ellips (kind of)
        r1 = np.linalg.norm(endpoints[0])
        r2 = np.linalg.norm(endpoints[1])
        for m in molecules:
            lithium_coord = m.sites[-1].coords
            dist1 = np.linalg.norm(endpoints[0] - lithium_coord)
            dist2 = np.linalg.norm(endpoints[1] - lithium_coord)
            li_r = np.linalg.norm(lithium_coord)
            new_radius = r2*dist1/(dist1 + dist2) + r1*dist2/(dist1 + dist2)
            translate_vector = lithium_coord/li_r*(new_radius-li_r)
            m.translate_sites([len(m)-1,], translate_vector)

        path_mol = mol.copy()
        for lithium in [m.sites[-1].coords for m in molecules]:
            path_mol.append('Li', lithium, validate_proximity=False)

        path_dir = 'path' + str(path_number)
        try:
            os.mkdir(path_dir)
        except FileExistsError:
            pass

        path_mol.to(fmt='xyz', filename=os.path.join(path_dir, 'path.xyz'))

        start_struct.to(fmt='poscar',
                        filename=os.path.join(path_dir,
                                              'init' + str(path_number)
                                              + '.vasp'))
        end_struct.to(fmt='poscar',
                      filename=os.path.join(path_dir,
                                            'final' + str(path_number)
                                            + '.vasp'))

        nw_input = nw.NwInput(mol, tasks, geometry_options=GEO_SETUP)
        nw_input.write_file(os.path.join(path_dir, 'input'))

        plot_images(molecules, filename=os.path.join(path_dir, 'path.neb'))

        path_number += 1


def set_up_structure(mol):
    """

    :param mol:
    :return:
    """
    # Find the distance between the origin and the furthest atom in the
    # molecule

    max_distance = 0
    for site in mol.sites:
        distance = np.linalg.norm(site.coords)
        if distance > max_distance:
            max_distance = distance

    # triple that, and use it to set up the lattice

    lat_const = 3 * max_distance
    lattice = pmg.Lattice(np.array([[lat_const, 0, 0], [0, lat_const, 0],
                                    [0, 0, lat_const]]))

    # Use the lattice in combination with the sites of the molecule to set up a
    # Structure, for which we can calculate the potential field

    struc = pmg.Structure(lattice, mol.species, mol.cart_coords,
                          coords_are_cartesian=True)

    return struc

def extract_docking_points(directory):
    """
    Extract the docking points from the docking directory.
    :param directory directory:
    :return:
    """
    dir_list = [d for d in os.listdir(directory)
                if os.path.isdir(os.path.join(directory, d))]

    docking_points = []
    for dir in dir_list:
        facet = Facet.from_file(os.path.join(dir, 'facet.json'))
        output = nw.NwOutput(os.path.join(dir, 'result.out'))
        final_mol = output.data[-1]['molecules'][-1]
        Li_coord = [site.coords for site in final_mol.sites
                    if site.specie == pmg.Element('Li')][0]
        docking_points.append((facet, Li_coord))

    return docking_points

def find_path_endpoints(mol, path, dock_points):
    """

    :return:
    """
    # Find the facet list of the molecule
    facet_list = mol.set_up_facet_list()

    # Set up the beginning and end sites of the various facets
    start_coords = None
    end_coords = None
    for facet in facet_list:

        if facet['surf_facet'] == path[0]:
            for point in dock_points:
                if facet['noneq_facet'] == point[0]:
                    noneq_dock_point = point[1]
            start_coords = facet['symmop'].inverse.operate(noneq_dock_point)

        if facet['surf_facet'] == path[1]:
            for point in dock_points:
                if facet['noneq_facet'] == point[0]:
                    noneq_dock_point = point[1]
            end_coords = facet['symmop'].inverse.operate(noneq_dock_point)

    endpoints = (start_coords, end_coords)

    return endpoints

def plot_images(molecules, filename='path.neb'):
    """

    :param molecules:
    :return:
    """
    file = io.open(filename, 'w')
    for structure in molecules:
        file.write(structure.to(fmt='xyz'))
        file.write('\n')


if __name__ == '__main__':
    main()


# # Old attempts to use the pymatgen path_finder package, abandoned because the
# # found paths were not really sensible
#
# # Read the charge density from the CHGCAR
# filename = sys.argv[2]
# charge_density = vasp.Chgcar.from_file(filename)
#
# # Use the pymatgen.analysis.path_finder FreeVolumePotential to find the
# # potential field of the molecule.
# # potential = path_finder.FreeVolumePotential(struc, [50, 50, 50]).get_v()
# # Since this doesn't seem to work, fuck it
#
# # Use the pymatgen.analysis.path_finder ChgcarPotential to find the
# # potential field of the molecule.
# potential = path_finder.ChgcarPotential(charge_density, normalize=True).get_v()
#
# # Use the NEBPathfinder class in combination with the potential
# # field to calculate the initial pathway
# neb = path_finder.NEBPathfinder(start_struct, end_struct,
#                                 [len(mol.sites), ], potential,
#                                 NIMAGES)
#
# neb.plot_images('path' + str(path_number) + '.vasp')
#

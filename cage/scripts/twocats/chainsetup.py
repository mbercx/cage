# encoding: utf-8
# Written for Python 3.6

import os
import math
import cage
import argparse

import numpy as np
import pymatgen as pmg
import pymatgen.io.nwchem as nwchem

"""

This script looks for the docking points for cations on the anion, then sets up
a set of calculations to study the 2D landscape between the non-equivalent
facets in this occupied anion.

The calculation is set up in the current directory, but the user must provide
the docking directory as input (for now).
"""

# TODO Make these parameters defaults, but allow the user to change them with arguments in the CLI
# TODO The current set up does not work very well for molecules that are not spherically shaped -> improve method set_up_landscape

OUTPUT_FILE = 'result.out'

# Facetsetup parameter
IGNORE = (pmg.Element('Li'), pmg.Element('Na'), pmg.Element('H'),
          pmg.Element('I'), pmg.Element('Br'), pmg.Element('Cl'))

# Landscape parameters

CATION = "Li"  # Cation to place on the landscape
# Distance endpoints between the center of the molecule and the cation
ENDPOINT_RADII = (3, 6)
# TODO For some reason, using the density to set the number of radii did not work. However, that seems much more sensible. Fix it.
N_RADII = 30  # Number of radius points for the landscape
ANGLE_DENSITY = 50  # Density of points along the angle coordinate

# Calculation parameters
# BASIS = {"*": "aug-pcseg-1"}
BASIS = {"*": "aug-cc-pVDZ"}

THEORY_SETUP = {"iterations": "300",
                "xc": "xpbe96 xpbe96",
                "direct": "",
                "smear": "0.01",
                "convergence energy": "1e-4",
                "convergence density": "1e-2",
                "convergence gradient": "1e-2",
                "convergence damp": "70"}

GEO_SETUP = {"noautosym", "noautoz", "nocenter", "units angstroms"}

ALT_SETUP = {}

DRIVER_SETUP = {"loose": "", "maxiter": "100"}

OPERATION = "energy"

# Input parsing
parser = argparse.ArgumentParser(description=
                                 "This script looks for the docking points for"
                                 " cations on the anion, then sets up a set of"
                                 " calculations to study the 2D landscape "
                                 "between the non-equivalent facets in this "
                                 "occupied anion. "
                                 "The calculation is set up in the current "
                                 "directory, but the user must provide "
                                 "the docking directory as input (for now).")

parser.add_argument('-I', metavar='cation', type=str,
                    help="The chemical symbol of the cation under study, i.e. "
                         "'Li', 'Na', ...")
parser.add_argument('-O', metavar='operation', type=str,
                    help="The Nwchem operation to use, i.e. 'energy' for only"
                         " total energy calculations, 'optimize' in case the "
                         "structure should be optimized first.")
parser.add_argument('dir', type=str,
                    help="The docking directory.")

args = vars(parser.parse_args())

try:
    dirname = args['dir']
except IndexError:
    raise IOError("No docking directory provided.")

for argument in args.keys():

    if args[argument]:

        if argument == 'I':
            CATION = args['I']
        elif argument == 'O':
            OPERATION = args['O']

def main():

    # Get the docking directories
    dir_list = [os.path.join(dirname, dir) for dir in os.listdir(dirname)
                if os.path.isdir(os.path.join(dirname, dir))]

    dock_number = 1

    for dir in dir_list:

        # Extract the occupied cage
        try:
            out = nwchem.NwOutput(os.path.join(dir, OUTPUT_FILE))
        except:
            print('Failed to extract output from ' + os.path.abspath(dir) +
                  '. Skipping...')
            continue

        mol = out.data[-1]['molecules'][-1]

        try:
            cat_coords = [site.coords for site in mol.sites
                          if site.specie == pmg.Element(CATION)][-1]
        except IndexError:
            raise IOError("Requested cation is not found in the molecule.")

        # Set up the occupied anion
        occmol = cage.facetsym.OccupiedCage.from_molecule(mol)
        occmol.find_surface_facets(ignore=IGNORE)
        dock = occmol.find_closest_facet(cat_coords)

        occmol.add_dock(dock, cation=None)
        occmol.find_surface_facets(IGNORE)

        total_mol = occmol.copy()

        # Find the chain paths
        paths = occmol.find_noneq_chain_paths()

        dock_dir = 'dock' + str(dock_number)

        try:
            os.mkdir(dock_dir)
        except FileExistsError:
            pass

        # For each facet, set up the calculation input files
        edge_number = 1

        for path in paths:

            # Set up the edge directory
            edge_dir = "edge" + str(edge_number)
            edge_dir = os.path.join(dock_dir, edge_dir)
            try:
                os.mkdir(edge_dir)
            except FileExistsError:
                pass

            # Write out the molecule and path facets to the edge directory
            mol.to(fmt="json", filename=os.path.join(edge_dir, "mol.json"))
            path[0].to(fmt="json", filename=os.path.join(edge_dir,
                                                         "init_facet.json"))
            path[1].to(fmt="json", filename=os.path.join(edge_dir,
                                                         "final_facet.json"))

            # Get copies so the originals aren't mutated
            edge_mol = occmol.copy()
            facet1 = path[0].copy()
            facet2 = path[1].copy()

            # Set up the landscape
            landscape = set_up_edge_landscape(facet1, facet2,
                                              endpoint_radii=ENDPOINT_RADII,
                                              number_of_radii=N_RADII,
                                              angle_density=ANGLE_DENSITY)

            # Get the molecule for each landscape point
            molecules = set_up_molecules(edge_mol, landscape, CATION)

            # Set up an xyz file to visualize the edge
            for point in landscape.points:
                try:
                    total_mol.append(pmg.Specie(CATION, 1), point,
                                     validate_proximity=False)
                    edge_mol.append(pmg.Specie(CATION, 1), point,
                                    validate_proximity=False)
                except ValueError:
                    pass

            edge_mol.to(fmt="xyz", filename=os.path.join(edge_dir, "edge.xyz"))

            # In case the molecules must be optimized, add the constraints and
            # optimization setup (DRIVER)
            if OPERATION == "optimize":
                fixed_facet = occmol.find_furthest_facet(landscape.center)
                ALT_SETUP["constraints"] = find_constraints(occmol,
                                                            fixed_facet)
                ALT_SETUP["driver"] = DRIVER_SETUP

            # Set up the task for the calculations
            tasks = [nwchem.NwTask(molecules[0].charge, None, BASIS,
                                   theory="dft",
                                   operation=OPERATION,
                                   theory_directives=THEORY_SETUP,
                                   alternate_directives=ALT_SETUP)]

            # Set up the input files
            study = cage.study.Study(molecules, tasks)
            study.set_up_input(edge_dir, sort_comp=False,
                               geometry_options=GEO_SETUP)

            edge_number += 1

        # Set up an xyz file with all the paths
        total_mol.to(fmt="xyz", filename=os.path.join(dock_dir,
                                                      "total_mol.xyz"))


###########
# METHODS #
###########


def set_up_edge_landscape(facet1, facet2, endpoint_radii=(2, 5),
                          number_of_radii=None, angle_density=50):
    """
    Set up the Landscape to study the energy landscape between two facets.

    The script creates a line, whose direction is determined by the vector
    connecting the origin to the center of the first facet. The line endpoints
    are determined by the endpoint_radii parameter. The line Landscape is then
    extended by rotation along the axis defined by the cross product of the two
    facet normals.

    :param facet1:
    :param facet2:
    :param endpoint_radii:
    :param number_of_radii:
    :param angle_density:
    :return:
    """
    line_vector = facet1.center/np.linalg.norm(facet1.center)
    lands = cage.landscape.Landscape.from_vertices(
        [line_vector * endpoint_radii[0], line_vector * endpoint_radii[1]],
        num=number_of_radii
    )
    axis = np.cross(facet1.normal, facet2.normal)
    angle = math.asin(np.linalg.norm(axis))
    axis = axis * angle / np.linalg.norm(axis)

    lands.extend_by_rotation(axis, angle_density, remove_endline=True)

    return lands


def set_up_molecules(mol, landscape, cation):
    """
    Set up the List of molecules from Landscape by adding the cation on the
    various points of the landscape.

    :param mol:
    :param landscape:
    :param cation:
    :return:
    """

    # Set up the cation Species
    cation = pmg.Specie(cation, +1)

    # Set up the List of Molecules
    molecules = []
    for point in landscape.points:
        configuration = mol.copy()

        # If carbon is in the molecule, the charge is -1, -2 otherwise
        # TODO Find a good way to calculate the charge, if possible
        if pmg.Element("C") in [site.specie for site in configuration.sites]:
            configuration.set_charge_and_spin(charge=0)
        else:
            configuration.set_charge_and_spin(charge=-1)

        # Add the cation to the Molecule
        try:
            configuration.append(cation, point)
            molecules.append(configuration)
        except ValueError:
            print("ValueError detected when appending the Li site. "
                  "Ignoring this point in the energy landscape.")

    return molecules


def find_constraints(mol, facet):
    """
    Find the constraints for the calculation, by fixing a facet and the cation.

    In order to make sure the Molecule does not rotate, the script fixes an
    entire facet of the Molecule. The cation is assumed to be the last element
    in the molecule.

    :return:
    """
    # Find the corresponding atom numbers in string format
    site_numbers = []
    for i in range(len(mol.sites)):
        if mol.sites[i] in facet.sites:
            site_numbers.append(str(i + 1) + " ")

    # Add the cation site
    site_numbers.append(str(len(mol.sites) + 1))

    # Join the strings of the site numbers
    site_numbers = "".join(site_numbers)

    # Return the constraints on the atoms
    return {"fix atom": site_numbers}


if __name__ == '__main__':
    main()

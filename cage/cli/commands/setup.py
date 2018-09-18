# encoding: utf-8
# Written for Python 3.6

import os
import io
import math

from cage.utils import unit_vector
from cage.core import Cage, OccupiedCage
from cage.study import Study
from cage.landscape import Landscape

import numpy as np
import pymatgen as pmg
import pymatgen.io.nwchem as nwchem

"""
Scripts to set up calculations for studying the geometry and energy landscapes 
around cage-like molecules.

These scripts were written specifically to study CxByXz molecules. In case you 
want to use it for anything else, they might need some adjustments.

"""

##############
# PARAMETERS #
##############

# Elements to ignore for the surface facet determination
IGNORE = (pmg.Element('Li'), pmg.Element('Na'), pmg.Element('Mg'),
          pmg.Element('H'), pmg.Element('I'), pmg.Element('Br'),
          pmg.Element('Cl'), pmg.Element('F'))

# Name of the output file
OUTPUT_FILE = 'result.out'

# Calculation parameters
THEORY_SETUP = {'iterations': '300',
                'xc': 'xpbe96 xpbe96',
                'direct': '',
                'smear': '0.01',
                'convergence energy': '1e-4',
                'convergence density': '1e-2',
                'convergence gradient': '1e-2',
                'convergence damp': '70'}

BASIS = {'*': "aug-cc-pVDZ"}

GEO_SETUP = {"noautosym", "noautoz", 'nocenter', "units angstroms"}

ALT_SETUP = {}

DRIVER_SETUP = {'loose': '', 'maxiter': '100'}


def optimize(filename, charge=None):
    """
    Set up a NwChem calculation to optimize a molecule structure.

    Args:
        filename (str): Structure file of the molecule.

    """

    try:
        # Load the Cage from the file provided
        anion = Cage.from_file(filename)
    except ValueError:
        # If that fails, try the POSCAR format
        anion = Cage.from_poscar(filename)

    # Set the charge for the molecule
    if charge:
        anion.set_charge_and_spin(charge=charge)
    elif pmg.Element('C') in [site.specie for site in anion.sites]:
        anion.set_charge_and_spin(charge=-1)
    else:
        anion.set_charge_and_spin(charge=-2)

    # Set up the geometry optimization
    tasks = [nwchem.NwTask(anion.charge, None, BASIS,
                           theory='dft',
                           operation="optimize",
                           theory_directives=THEORY_SETUP,
                           alternate_directives=ALT_SETUP)]

    # Set the driver settings for the optimization
    ALT_SETUP["driver"] = DRIVER_SETUP

    # Set up input
    try:
        os.mkdir('optimize')
    except FileExistsError:
        pass

    nw_input = nwchem.NwInput(anion, tasks, geometry_options=GEO_SETUP)
    nw_input.write_file(os.path.join('optimize', 'input'))


def docksetup(filename, cation, distance, facets, verbose):
    if verbose:
        print("Loading structure file " + filename + "...")

    try:
        # Load the POSCAR into a Cage
        anion = Cage.from_poscar(filename)
    except ValueError:
        # If that fails, try other file formats supported by pymatgen
        anion = Cage.from_file(filename)

    if verbose:
        print("Setting up surface facets...")

    anion.find_surface_facets(ignore=IGNORE)

    if verbose:
        print("Found " + str(len(anion.facets)) + " facets.")

    if verbose:
        print("Studying symmetry...")

    if facets == tuple:
        # Find the non-equivalent facets, and use them for the docking sites
        facets = anion.find_noneq_facets()
    else:
        facets = [anion.facets[index] for index in facets]

    if verbose:
        print("Found " + str(len(facets)) + " non-equivalent facets.")

    docking_dir = 'docking'
    try:
        os.mkdir(docking_dir)
    except FileExistsError:
        pass

    if verbose:
        print("Setting up docking calculations...")

    # For each docking point, set up the calculation input file
    dock_number = 1
    for neq_facet in facets:

        if verbose:
            print("")
            print("DOCK " + str(dock_number))
            print("------")
            print("")
            print("Setting up cation site...")

        # Set up the initial cation site
        mol = anion.copy()
        mol.append(pmg.Specie(cation, 1), neq_facet.center +
                   distance * neq_facet.normal)
        mol.find_surface_facets(ignore=IGNORE)

        if verbose:
            print("Finding constraints for the calculation...")

        # Constrain the opposite facet --> Why?
        # far_facet = mol.find_farthest_facet(neq_facet.center)
        # ALT_SETUP['constraints'] = find_constraints(mol, far_facet.sites)

        # Set the driver settings for the optimization
        ALT_SETUP["driver"] = DRIVER_SETUP

        # Set the charge for the molecule
        if pmg.Element('C') in [site.specie for site in mol.sites]:
            mol.set_charge_and_spin(charge=0)
        else:
            mol.set_charge_and_spin(charge=-1)

        if verbose:
            print("Setting up the task...")

        # Set up the task for the calculations
        tasks = [nwchem.NwTask(mol.charge, None, BASIS,
                               theory='dft',
                               operation="optimize",
                               theory_directives=THEORY_SETUP,
                               alternate_directives=ALT_SETUP)]

        dock_dir = os.path.join(docking_dir, 'dock' + str(dock_number))

        while True:
            try:
                os.mkdir(dock_dir)
                break
            except FileExistsError:
                dock_number += 1
                dock_dir = os.path.join(docking_dir, 'dock' + str(dock_number))

        if verbose:
            print("Setting up the input file...")

        # Set up input
        nw_input = nwchem.NwInput(mol, tasks, geometry_options=GEO_SETUP)

        if verbose:
            print("Writing input file for dock " + str(dock_number) + "...")

        nw_input.write_file(os.path.join(dock_dir, 'input'))

        # Write out a facet json file
        neq_facet.to(fmt='json', filename=os.path.join(dock_dir, 'facet.json'))

        # Write a xyz file of the molecule with the docked cation
        mol.to(fmt='xyz', filename=os.path.join(dock_dir, 'dock.xyz'))

        dock_number += 1


def chainsetup(filename, cation, facets, operation, end_radii, nradii,
               adensity):
    """
    Set up calculations to study the 2D edge landscapes between a chain of
    facets.

    Args:
        filename:
        cation:
        facets:
        operation:
        end_radii:
        nradii:
        adensity:

    Returns:

    """

    # Load the Cage from the file
    try:
        # If that fails, try other file formats supported by pymatgen
        anion = Cage.from_file(filename)
    except ValueError:
        # If that fails, try the VASP POSCAR format
        anion = Cage.from_poscar(filename)

    # Center the anion around the origin
    anion.center()

    # Find the chain edges, i.e. the paths between the edge sharing facets of
    # the chain of non-equivalent facets.
    anion.find_surface_facets(ignore=IGNORE)

    if not facets == tuple:
        chosen_facets = [anion.facets[index] for index in facets]
        edges = anion.find_noneq_chain_links(chosen_facets)
    else:
        edges = anion.find_noneq_chain_links()

    total_mol = anion.copy()

    chain_dir = 'chain_' + operation
    try:
        os.mkdir(chain_dir)
    except FileExistsError:
        pass

    # For each edge, set up the calculation input files
    edge_number = 1

    for edge in edges:

        # Set up the edge directory
        edge_dir = os.path.join(chain_dir, "edge" + str(edge_number))
        try:
            os.mkdir(edge_dir)
        except FileExistsError:
            pass

        # Write out the molecule and path facets to the edge directory
        anion.to(fmt="json", filename=os.path.join(edge_dir, "molecule.json"))
        edge[0].to(fmt="json", filename=os.path.join(edge_dir,
                                                     "init_facet.json"))
        edge[1].to(fmt="json", filename=os.path.join(edge_dir,
                                                     "final_facet.json"))

        # Get copies so the originals aren't mutated
        edge_mol = anion.copy()
        facet1 = edge[0].copy()
        facet2 = edge[1].copy()

        if edge == edges[-1]:
            remove_endline = False
        else:
            remove_endline = True

        # Set up the landscape
        landscape = set_up_edge_landscape(facet1, facet2,
                                          endpoint_radii=end_radii,
                                          number_of_radii=nradii,
                                          angle_density=adensity,
                                          remove_endline=remove_endline)

        # Get the molecule for each landscape point
        molecules = set_up_molecules(edge_mol, landscape, cation)

        # Set up an xyz file to visualize the edge and total landscape
        for point in landscape.points:
            try:
                total_mol.append(pmg.Specie(cation, 1), point,
                                 validate_proximity=False)
                edge_mol.append(pmg.Specie(cation, 1), point,
                                validate_proximity=False)
            except ValueError:
                pass

        edge_mol.to(fmt="xyz", filename=os.path.join(edge_dir, "edge.xyz"))

        # In case the molecules must be optimized, add the constraints and
        # optimization setup (DRIVER)
        if operation == "optimize":
            far_facet = anion.find_farthest_facet(landscape.center)
            constraints = find_constraints(anion, far_facet.sites)
            constraints['fix atom'] += ' ' + str(len(anion.sites) + 1)
            ALT_SETUP['constraints'] = constraints
            ALT_SETUP["driver"] = DRIVER_SETUP

        # Set up the task for the calculations
        tasks = [nwchem.NwTask(molecules[0].charge, None, BASIS,
                               theory="dft",
                               operation=operation,
                               theory_directives=THEORY_SETUP,
                               alternate_directives=ALT_SETUP)]

        # Set up the input files
        study = Study(molecules, tasks)
        study.set_up_input(edge_dir, sort_comp=False,
                           geometry_options=GEO_SETUP)

        edge_number += 1

    # Set up an xyz file with all the paths
    total_mol.to(fmt="xyz", filename=os.path.join(chain_dir, "total_mol.xyz"))


def pathsetup(filename, cation, distance, facets, edges):
    # Load the Cage from the file
    try:
        # Load the POSCAR into a Cage
        anion = Cage.from_poscar(filename)
    except ValueError:
        # If that fails, try other file formats supported by pymatgen
        anion = Cage.from_file(filename)

    # TODO Find charge automatically
    if pmg.Element('C') in [site.specie for site in anion.sites]:
        anion.set_charge_and_spin(charge=0)
    else:
        anion.set_charge_and_spin(charge=-1)

    # Find the surface facets
    anion.find_surface_facets(ignore=IGNORE)

    if facets == tuple:
        # Find the paths, i.e. the List of facet combinations
        paths = anion.find_facet_links(share_edge=edges)
    else:
        chosen_facets = [anion.facets[index] for index in facets]
        paths = anion.find_noneq_chain_links(chosen_facets)

    # Set the driver settings for the optimization
    ALT_SETUP["driver"] = DRIVER_SETUP

    tasks = [nwchem.NwTask(anion.charge, None, BASIS, theory='dft',
                           operation="optimize",
                           theory_directives=THEORY_SETUP,
                           alternate_directives=ALT_SETUP)]

    paths_dir = 'paths'
    try:
        os.mkdir(paths_dir)
    except FileExistsError:
        pass

    path_number = 1
    for path in paths:

        # Set up the start and end molecules
        start_molecule = anion.copy()
        end_molecule = anion.copy()

        start_molecule.append(pmg.Specie(cation, 1), path[0].center +
                              distance * path[0].normal)
        end_molecule.append(pmg.Specie(cation, 1), path[1].center +
                            distance * path[1].normal)

        # Make the path directory
        path_dir = os.path.join(paths_dir, 'path' + str(path_number))
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
        nwchem.NwInput(start_molecule, tasks, geometry_options=GEO_SETUP) \
            .write_file(os.path.join(path_dir, 'start', 'input'))
        nwchem.NwInput(end_molecule, tasks, geometry_options=GEO_SETUP) \
            .write_file(os.path.join(path_dir, 'end', 'input'))

        path_number += 1


def nebsetup(paths_dir, nimages):
    directory = os.path.abspath(paths_dir)

    dir_list = [os.path.join(directory, d) for d in os.listdir(directory)
                if os.path.isdir(os.path.join(directory, d))]

    for directory in dir_list:

        # Get the output from the docking calculations
        try:
            start_data = nwchem.NwOutput(os.path.join(directory, 'start',
                                                      OUTPUT_FILE)).data[-1]
            end_data = nwchem.NwOutput(os.path.join(directory, 'end',
                                                    OUTPUT_FILE)).data[-1]
        except:
            print('Failed to extract start and end images in ' + directory +
                  ', ignoring directory...')
            continue

        # Set up the structures
        start_molecule = start_data['molecules'][-1]
        end_molecule = end_data['molecules'][-1]

        lattice = set_up_lattice(start_molecule)

        start_struct = pmg.Structure(lattice, start_molecule.species,
                                     start_molecule.cart_coords,
                                     coords_are_cartesian=True)
        end_struct = pmg.Structure(lattice, end_molecule.species,
                                   end_molecule.cart_coords,
                                   coords_are_cartesian=True)

        # Interpolate to find the images
        structures = start_struct.interpolate(end_struct, nimages=nimages)

        molecules = [pmg.Molecule(struct.species, struct.cart_coords)
                     for struct in structures]

        # Move the lithium positions on an ellipsoid (kind of)
        r1 = np.linalg.norm(start_molecule.sites[-1].coords)
        r2 = np.linalg.norm(end_molecule.sites[-1].coords)
        for m in molecules:
            lithium_coord = m.sites[-1].coords
            dist1 = np.linalg.norm(start_molecule.sites[-1].coords
                                   - lithium_coord)
            dist2 = np.linalg.norm(end_molecule.sites[-1].coords
                                   - lithium_coord)
            li_r = np.linalg.norm(lithium_coord)
            new_radius = r2 * dist1 / (dist1 + dist2) + r1 * dist2 / (
            dist1 + dist2)
            translate_vector = lithium_coord / li_r * (new_radius - li_r)
            m.translate_sites([len(m) - 1, ], translate_vector)

        # Make a path molecule to show the interpolation
        path_mol = start_molecule.copy()
        for m in molecules[1:]:
            for coord in [(site.coords, site.specie) for site in m.sites]:
                path_mol.append(coord[1], coord[0], validate_proximity=False)

        path_mol.to(fmt='xyz', filename=os.path.join(directory, 'path.xyz'))

        # Set up the input file
        ALT_SETUP["neb"] = DRIVER_SETUP
        ALT_SETUP["neb"].update({"nbeads": nimages,
                                 "xyz_path": "path.neb",
                                 "print_shift": 5})

        # Set the charge for the molecule
        # TODO CALCULATE THE CHARGE SOMEHOW
        if pmg.Element('C') in [site.specie for site in start_molecule.sites]:
            start_molecule.set_charge_and_spin(charge=0)
        else:
            start_molecule.set_charge_and_spin(charge=-1)

        tasks = [nwchem.NwTask(start_molecule.charge, None, BASIS,
                               theory='dft',
                               operation="neb",
                               theory_directives=THEORY_SETUP,
                               alternate_directives=ALT_SETUP), ]

        nw_input = nwchem.NwInput(start_molecule, tasks,
                                  geometry_options=GEO_SETUP)
        nw_input.write_file(os.path.join(directory, 'input'))

        plot_images(molecules, filename=os.path.join(directory, 'path.neb'))


def reference(facet_index, filename, cation, end_radius, start_radius=4.0,
              nradii=10, verbose=False):
    """
    Set up a calculation to determine a reference energy for an anion,
    based on a 1D landscape perpendicular to a chosen facet.

    Args:
        facet_index:
        filename:
        cation:
        end_radius:
        start_radius:
        nradii:
        verbose:

    Returns:

    """
    if verbose:
        print("Loading structure file " + filename + "...")

    try:
        # Load the structure file into a Cage object
        anion = OccupiedCage.from_file(filename)
    except ValueError:
        # If the conventional formats fail, try using the from_poscar method.
        anion = Cage.from_poscar(filename)

    if verbose:
        print("Setting up surface facets...")

    anion.center()
    anion.find_surface_facets(ignore=IGNORE)

    if verbose:
        print("Found " + str(len(anion.facets)) + " facets.")

    reference_dir = 'reference_fa' + str(facet_index)
    try:
        os.mkdir(reference_dir)
    except FileExistsError:
        pass

    if verbose:
        print("Setting up reference calculations...")

    reference_facet = anion.facets[facet_index]

    # Write out a facet json file
    reference_facet.to(fmt='json', filename=os.path.join(reference_dir,
                                                         'facet.json'))

    for radius in np.linspace(start_radius, end_radius, nradii):

        if verbose:
            print("Adding cation to reference facet at distance = "
                  + str(radius) + " Angstroms.")

        # Set up the cation site
        mol = anion.copy()
        mol.append(pmg.Specie(cation, 1),
                   unit_vector(reference_facet.center) * radius)

        mol.set_charge_and_spin(charge=anion.charge + 1)

        if verbose:
            print("Setting up the task...")

        # Set up the task for the calculation
        tasks = [nwchem.NwTask(mol.charge, None, BASIS,
                               theory='dft',
                               operation="energy",
                               theory_directives=THEORY_SETUP,
                               alternate_directives=ALT_SETUP)]

        if verbose:
            print("Setting up the input file...")

        # Set up input
        nw_input = nwchem.NwInput(mol, tasks, geometry_options=GEO_SETUP)

        calculation_dir = os.path.join(reference_dir, "radius="
                                       + '{0:.2f}'.format(radius))

        try:
            os.mkdir(os.path.join(calculation_dir))
        except FileExistsError:
            pass

        nw_input.write_file(os.path.join(calculation_dir, 'input'))

        # Write a xyz file of the molecule with the reference cation
        mol.to(fmt='xyz', filename=os.path.join(calculation_dir,
                                                'geometry.xyz'))


def spheresetup(filename, cation, radius, axis=None, density=20):

    anion = Cage.from_file(filename)
    anion.center()

    if not axis:
        print(anion)
        axis = input("\n Please choose the direction of the main axis (vector "
                   "format, spaces between the coordinates): ")

        axis = np.array([float(number) for number in axis.split(" ")])

    sphere_landscape = Landscape.create_sphere(
        radius=radius,
        axis=axis,
        density=density
    )

    molecule_list = set_up_molecules(anion, sphere_landscape, cation)

    # Set up the task for the calculations
    tasks = [nwchem.NwTask(molecule_list[0].charge, None, BASIS,
                           theory="dft",
                           operation="energy",
                           theory_directives=THEORY_SETUP,
                           alternate_directives=ALT_SETUP)]

    current_dir = os.getcwd()
    sphere_dir = os.path.join(current_dir, "sphere_" + str(radius))

    try:
        os.mkdir(sphere_dir)
    except FileExistsError:
        pass

    # Set up the input files
    study = Study(molecule_list, tasks)
    study.set_up_input(sphere_dir, sort_comp=False, geometry_options=GEO_SETUP)

    total_mol = anion.copy()

    # Set up an xyz file to visualize the edge and total landscape
    for point in sphere_landscape.points:
        try:
            total_mol.append(pmg.Specie(cation, 1), point,
                             validate_proximity=False)
        except ValueError:
            pass

    # Write the axis coordinates to a file
    with open(os.path.join(sphere_dir, "axis.xyz")) as file:
        file.write(axis.tostring())

    total_mol.to(fmt="xyz", filename=os.path.join(sphere_dir, "total_mol.xyz"))


def twocat_chainsetup(dock_dir, cation, operation, endradii, nradii, adensity,
                      tolerance, verbose):
    # Get the docking directories
    dir_list = [os.path.join(dock_dir, directory)
                for directory in os.listdir(dock_dir)
                if os.path.isdir(os.path.join(dock_dir, directory))]

    dock_number = 1

    chain_dir = 'twocat_chain'
    try:
        os.mkdir(chain_dir)
    except FileExistsError:
        pass

    for directory in dir_list:

        # Extract the occupied cage
        try:
            out = nwchem.NwOutput(os.path.join(directory, OUTPUT_FILE))
        except:
            print('Failed to extract output from ' + os.path.abspath(directory)
                  + '. Skipping...')
            continue

        # Check if calculation has completed successfully
        if out.data[-1]['task_time'] == 0:
            print("WARNING: No task time information found in output file: "
                  + os.path.join(directory, OUTPUT_FILE) + "\n" +
                  "The calculation might not have completed successfully.")

        for data in out.data:
            if data['has_error']:
                print("WARNING: Error found in one of the calculations of "
                      + os.path.join(directory, OUTPUT_FILE))

        # The final molecules in the optimization is used
        mol = out.data[-1]['molecules'][-1]

        try:
            cat_coords = [site.coords for site in mol.sites
                          if site.specie == pmg.Element(cation)][-1]
        except IndexError:
            raise IOError("Requested cation is not found in the molecule.")

        # Set up the occupied anion
        occmol = OccupiedCage.from_molecule(mol)
        occmol.center()
        occmol.find_surface_facets(ignore=IGNORE)
        dock = occmol.find_closest_facet(cat_coords)

        occmol.add_dock(dock, cation=None)
        occmol.find_surface_facets(ignore=IGNORE)

        total_mol = occmol.copy()

        # Find the chain paths
        paths = occmol.find_noneq_chain_links(symm_tol=tolerance,
                                              verbose=verbose)

        dock_dir = os.path.join(chain_dir, 'dock' + str(dock_number))

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
                                              endpoint_radii=endradii,
                                              number_of_radii=nradii,
                                              angle_density=adensity)

            # Get the molecule for each landscape point
            molecules = set_up_molecules(edge_mol, landscape, cation)

            # Set up an xyz file to visualize the edge
            for point in landscape.points:
                try:
                    total_mol.append(pmg.Specie(cation, 1), point,
                                     validate_proximity=False)
                    edge_mol.append(pmg.Specie(cation, 1), point,
                                    validate_proximity=False)
                except ValueError:
                    pass

            edge_mol.to(fmt="xyz", filename=os.path.join(edge_dir, "edge.xyz"))

            # In case the molecules must be optimized, add the constraints and
            # optimization setup (DRIVER)
            if operation == "optimize":
                far_facet = occmol.find_farthest_facet(landscape.center)
                ALT_SETUP["constraints"] = find_constraints(occmol,
                                                            far_facet.sites)
                ALT_SETUP['constraints'].join(
                    ' ' + str(len(mol.sites)))  # cation
                ALT_SETUP["driver"] = DRIVER_SETUP

            # Set up the task for the calculations
            tasks = [nwchem.NwTask(molecules[0].charge, None, BASIS,
                                   theory="dft",
                                   operation=operation,
                                   theory_directives=THEORY_SETUP,
                                   alternate_directives=ALT_SETUP)]

            # Set up the input files
            study = Study(molecules, tasks)
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
                          number_of_radii=None, angle_density=50,
                          remove_endline=True):
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
    line_vector1 = facet1.center / np.linalg.norm(facet1.center)
    line_vector2 = facet2.center / np.linalg.norm(facet2.center)
    lands = Landscape.from_vertices(
        [line_vector1 * endpoint_radii[0], line_vector1 * endpoint_radii[1]],
        num=number_of_radii
    )

    axis = np.cross(line_vector1, line_vector2)
    angle = math.asin(np.linalg.norm(axis))
    axis = axis * angle / np.linalg.norm(axis)

    lands.extend_by_rotation(axis, angle_density,
                             remove_endline=remove_endline)

    return lands


def set_up_molecules(molecule, landscape, cation):
    """
    Set up the List of molecules from Landscape by adding the cation on the
    various points of the landscape.

    :param molecule:
    :param landscape:
    :param cation:
    :return:
    """

    # Check if the molecule has a charge zero
    if molecule.charge == 0:
        print("WARNING: Found charge equal to zero for the anion molecule. "
              "This may be because there was no charge information on the "
              "structure file used. The charge will be set to -1 by "
              "assumption. Please provide .json structure files.")
        molecule.set_charge_and_spin(charge=-1)

    # Set up the cation Species
    cation = pmg.Specie(cation, oxidation_state=1)

    # Set up the List of Molecules
    molecules = []
    for point in landscape.points:
        configuration = molecule.copy()

        # Add the cation to the Molecule
        try:
            configuration.append(cation, point)
            configuration.set_charge_and_spin(charge=molecule.charge + 1)
            molecules.append(configuration)

        except ValueError:
            print("ValueError detected when appending the cation site. "
                  "Ignoring this point in the energy landscape.")

    return molecules


def find_constraints(mol, sites):
    """
    Find the NwChem constraints for a selection of sites on a molecule.

    """

    # Find the corresponding atom numbers in string format
    site_numbers = []
    for i in range(len(mol.sites)):
        if mol.sites[i] in sites:
            site_numbers.append(str(i + 1) + ' ')

    site_numbers = ''.join(site_numbers)

    # Set up the constraints on the atoms
    return {'fix atom': site_numbers}


def set_up_lattice(mol):
    """

    :param mol:
    :return:
    """
    # Find the distance between the origin and the farthest atom in the
    # molecule

    max_distance = 0
    for site in mol.sites:
        distance = np.linalg.norm(site.coords)
        if distance > max_distance:
            max_distance = distance

    # triple that, and use it to set up the lattice

    lat_const = int(3 * max_distance)
    lattice = pmg.Lattice(np.array([[lat_const, 0, 0], [0, lat_const, 0],
                                    [0, 0, lat_const]]))

    return lattice


def plot_images(molecules, filename='path.neb'):
    """

    :param molecules:
    :return:
    """
    file = io.open(filename, 'w')
    for structure in molecules:
        file.write(structure.to(fmt='xyz'))
        file.write('\n')

# Encoding" utf-8

import sys
import os
import io

import pymatgen as pmg
import numpy as np
import pymatgen.io.nwchem as nw

"""
Script to set up the directories and input files to calculate the minimum
energy paths for a cage molecule.
"""

# Parameters

NIMAGES = 10  # Standard number of images along every path
OUT_FILE = 'result.out'

# Calculation parameters
#BASIS = {'*': "aug-pcseg-1"}
BASIS = {'*': "aug-cc-pVDZ"}

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
OPERATION = 'neb'

def main():

    if sys.argv[1]:
        directory = os.path.abspath(sys.argv[1])
    else:
        directory = os.path.abspath('.')

    dir_list = [os.path.join(directory, d) for d in os.listdir(directory)
                if os.path.isdir(os.path.join(directory, d))]

    for dir in dir_list:

        # Get the output from the docking calculations
        try:
            start_data = nw.NwOutput(os.path.join(dir, 'start',
                                                  OUT_FILE)).data[-1]
            end_data = nw.NwOutput(os.path.join(dir, 'end',
                                                OUT_FILE)).data[-1]
        except:
            print('Failed to extract output in ' + dir +', ignoring directory...')
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
        structures = start_struct.interpolate(end_struct, nimages=NIMAGES)

        molecules = [pmg.Molecule(struct.species, struct.cart_coords)
                     for struct in structures]

        # Move the lithium positions on an ellips (kind of)
        r1 = np.linalg.norm(start_molecule.sites[-1].coords)
        r2 = np.linalg.norm(end_molecule.sites[-1].coords)
        for m in molecules:
            lithium_coord = m.sites[-1].coords
            dist1 = np.linalg.norm(start_molecule.sites[-1].coords
                                   - lithium_coord)
            dist2 = np.linalg.norm(end_molecule.sites[-1].coords
                                   - lithium_coord)
            li_r = np.linalg.norm(lithium_coord)
            new_radius = r2*dist1/(dist1 + dist2) + r1*dist2/(dist1 + dist2)
            translate_vector = lithium_coord/li_r*(new_radius-li_r)
            m.translate_sites([len(m)-1,], translate_vector)

        # Make a path molecule to show the interpolation
        path_mol = start_molecule.copy()
        for m in molecules[1:]:
            for coord in [(site.coords, site.specie) for site in m.sites]:
                path_mol.append(coord[1], coord[0], validate_proximity=False)

        path_mol.to(fmt='xyz', filename=os.path.join(dir, 'path.xyz'))

        # Set up the input file

        # Set the charge for the molecule
        # TODO CALCULATE THE CHARGE SOMEHOW
        if pmg.Element('C') in [site.specie for site in start_molecule.sites]:
            start_molecule.set_charge_and_spin(charge=0)
        else:
            start_molecule.set_charge_and_spin(charge=-1)

        tasks = [nw.NwTask(start_molecule.charge, None, BASIS, theory='dft',
                           operation=OPERATION,
                           theory_directives=THEORY_SETUP,
                           alternate_directives=ALT_SETUP), ]

        nw_input = nw.NwInput(start_molecule, tasks, geometry_options=GEO_SETUP)
        nw_input.write_file(os.path.join(dir, 'input'))

        plot_images(molecules, filename=os.path.join(dir, 'path.neb'))

def set_up_lattice(mol):
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

if __name__ == '__main__':
    main()
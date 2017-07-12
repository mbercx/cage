dir_list = [os.path.join(current_dir, d) for d in
            os.listdir(current_dir)
            if os.path.isdir(os.path.join(current_dir, d))]

if 'dock1' in dir_list and 'dock2' in dir_list:
    # Extract the final molecules from the docking directories
    pass

 # Find the docking points for the path
        path_molecules = find_path_molecules(mol, path, dock_molecules)

        # Set up the start and end structures
        start_struct = set_up_structure(path_molecules[0].copy())

        end_struct = set_up_structure(path_molecules[1].copy())

        # Interpolate to find the images
        structures = start_struct.interpolate(end_struct, nimages=NIMAGES)

        molecules = [pmg.Molecule(struct.species, struct.cart_coords)
                     for struct in structures]

        # Move the lithium positions on an ellips (kind of)
        r1 = np.linalg.norm(path_molecules[0].sites[-1].coords)
        r2 = np.linalg.norm(path_molecules[1].sites[-1].coords)
        for m in molecules:
            lithium_coord = m.sites[-1].coords
            dist1 = np.linalg.norm(path_molecules[0].sites[-1].coords
                                   - lithium_coord)
            dist2 = np.linalg.norm(path_molecules[1].sites[-1].coords
                                   - lithium_coord)
            li_r = np.linalg.norm(lithium_coord)
            new_radius = r2*dist1/(dist1 + dist2) + r1*dist2/(dist1 + dist2)
            translate_vector = lithium_coord/li_r*(new_radius-li_r)
            m.translate_sites([len(m)-1,], translate_vector)

        path_mol = mol.copy()
        for m in molecules:
            for coord in [(site.coords, site.specie) for site in m.sites]:
                path_mol.append(coord[1], coord[0], validate_proximity=False)

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

def extract_docking_molecules(directory):
    """
    Extract the docking points from the docking directory.
    :param directory directory:
    :return:
    """
    dir_list = [os.path.join(os.path.abspath(directory), d)
                for d in os.listdir(directory)
                if os.path.isdir(os.path.join(directory, d))]

    print(dir_list)

    docking_molecules = []
    for dir in dir_list:
        facet = Facet.from_file(os.path.join(dir, 'facet.json'))
        output = nw.NwOutput(os.path.join(dir, 'result.out'))
        final_mol = output.data[-1]['molecules'][-1]
        docking_molecules.append((facet, final_mol))

    return docking_molecules

def plot_images(molecules, filename='path.neb'):
    """

    :param molecules:
    :return:
    """
    file = io.open(filename, 'w')
    for structure in molecules:
        file.write(structure.to(fmt='xyz'))
        file.write('\n')
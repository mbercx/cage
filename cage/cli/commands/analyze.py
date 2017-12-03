# Encoding: utf-8

import os
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

from cage.landscape import LandscapeAnalyzer
from cage.core import Facet
from cage import utils

"""
Analysis tools for the data produced by calculation set up using the cage setup
tools. Requires expansion.

"""

START_FACET = 0  # 0 or 1 -> determines facet to start chain from

def landscape_analysis(lands_dir, cation, energy_range, interp_mesh, end_radii,
                       contour_levels, verbose):

    lands_dir = os.path.abspath(lands_dir)

    dir_list = [os.path.abspath(os.path.join(lands_dir, file))
                for file in os.listdir(lands_dir)
                if os.path.isdir(os.path.join(lands_dir, file))]

    if verbose:
        print("Extracting Landscape data from " + lands_dir + "...")

    # Extract the data from the subdirectories of the landscape directory
    chain = {}

    for directory in dir_list:

        if os.path.isfile(os.path.join(directory, 'landscape.json')):
            if os.path.isfile(os.path.join(directory, 'init_facet.json')) \
                and os.path.isfile(os.path.join(directory,
                                                'final_facet.json')):
                if verbose:
                    print("Found data in " + directory + ".")

                chain[directory] = {
                    'landscape': LandscapeAnalyzer.from_file(
                        os.path.join(directory, 'landscape.json')),
                    "link": [Facet.from_file(os.path.join(directory,
                                                          'init_facet.json')),
                             Facet.from_file(
                                 os.path.join(directory, 'final_facet.json'))]
                }
            else:

                if verbose:
                    print("No facet information found in " + directory + ".")
        else:
            print("No landscape data found in " + directory + ".")

    if verbose:
        print("Found " + str(len(chain.keys())) +
              " directories with landscape data.")
        print("")
        print("Analyzing landscape data...")

    for data in chain.keys():

        # Extract the total energies from the calculations, and assign their
        # coordinates.

        chain[data]['landscape'].analyze_cation_energies(
            facet=chain[data]['link'][0],
            cation=cation
        )

    # Check if the chain is properly constructed
    broken_chain = False
    test_chain = [chain[data]['link'] for data in chain.keys()]
    for index in range(len(test_chain)-1):
        if test_chain[index][1] != test_chain[index + 1][0]:
            broken_chain = True

    if broken_chain:
        chain_landscapes, facet_chain = reconstruct_chain(chain,
                                                          verbose=verbose)
    else:
        chain_links = [chain[data]['link'] for data in chain.keys()]
        facet_chain = [chain_links[0][0]]
        for link in chain_links:
            facet_chain.append(link[1])
        chain_landscapes = [chain[data]['landscape'] for data in chain.keys()]

    # Interpolate the landscapes to a uniform mesh
    if verbose:
        print('-----------')
        print("Finding proper radii and angles...")

    # If no radii were provided
    if end_radii == (0.0, 0.0):
        # Find the proper endradii automatically
        min_max_radius = 1e6
        max_min_radius = 0
        for landscape in chain_landscapes:
            rmax = landscape.datapoints['Distance'].max()
            if rmax < min_max_radius:
                min_max_radius = rmax
            rmin = landscape.datapoints['Distance'].min()
            if rmin > max_min_radius:
                max_min_radius = rmin
    else:
        # Use the radii provided by the user
        max_min_radius = end_radii[0]
        min_max_radius = end_radii[1]

    if verbose:
        print("")
        print('Largest minimal radius = ' + str(max_min_radius))
        print('Smallest maximal radius = ' + str(min_max_radius))

    # Adjust the angles to set up one angle coordinate for all edges

    facet_angles = [0, ]

    for index in range(1,len(facet_chain)):

        chain_landscapes[index-1].datapoints['Angle'] += facet_angles[-1]

        facet_angles.append(facet_angles[-1] +
                            utils.angle_between(facet_chain[index - 1].center,
                                                facet_chain[index].center))

    if verbose:
        print("")
        print("Facet Angles:")
        print(facet_angles)

    # facet_angles = [0, landscape_chain[0].datapoints['Angle'].max()]
    #
    # for lands in landscape_chain[1:]:
    #
    #     if verbose:
    #         print('Maximum angle = ' + str(facet_angles[-1]))
    #
    #     lands.datapoints['Angle'] += facet_angles[-1]
    #     facet_angles.append(lands.datapoints['Angle'].max())

    all_radii = []
    all_angles = []
    all_energy = []

    # TODO There needs to be a better way of interpolating this...
    # Interpolate the landscapes
    for landscape in chain_landscapes:

        data = landscape.datapoints

        data['Distance'] = np.round(data['Distance'], 5)
        data = np.sort(data, order=['Distance', 'Angle'])

        # Find the number of radii and angles
        r_init = data['Distance'][0]
        nangles = 1
        while abs(data['Distance'][nangles] - r_init) < 1e-5:
            nangles += 1
        nradii = int(len(data) / nangles)

        if verbose:
            print('')
            print('-----------')
            print('Number of Angles = ' + str(nangles))
            print('Number of Radii = ' + str(nradii))

        # Get the right format for the data
        radii = data['Distance'].reshape(nradii, nangles)
        angles = data['Angle'].reshape(nradii, nangles)
        energy = data['Energy'].reshape(nradii, nangles)

        if verbose:
            print('Shape angles: ' + str(angles.shape))
            print('Shape radii: ' + str(radii.shape))
            print('Shape energy: ' + str(energy.shape))

            print("Minimum Angle = " + str(angles.min()))
            print("Maximum Angle = " + str(angles.max()))

        if interp_mesh == (0.0, 0.0):
            new_radii = np.transpose(radii)
            new_angles = np.transpose(angles)
            new_energy = np.transpose(energy)
        else:
            new_angles, new_radii = np.mgrid[
                                    angles.min():angles.max() + interp_mesh[0]
                                    :interp_mesh[0],
                                    max_min_radius:min_max_radius
                                                   + interp_mesh[1]
                                    :interp_mesh[1]
                                    ]

            if verbose:
                print('-------------')
                print('Shape new_angles: ' + str(new_angles.shape))
                print('Shape new_radii: ' + str(new_radii.shape))
                print("New Minimum Angle = " + str(new_angles.min()))
                print("New Maximum Angle = " + str(new_angles.max()))

            tck = interpolate.bisplrep(angles, radii, energy, s=0.1)

            new_energy = interpolate.bisplev(new_angles[:, 0], new_radii[0, :],
                                             tck)

        all_radii.append(new_radii)
        all_angles.append(new_angles)
        all_energy.append(new_energy)

    total_radii = np.concatenate(tuple(all_radii))
    total_angles = np.concatenate(tuple(all_angles))
    total_energy = np.concatenate(tuple(all_energy))
    total_energy -= total_energy.min()

    # If no energy range is specified by the user, take (min, max)
    if energy_range == (0.0, 0.0):
        energy_range = (total_energy.min(), total_energy.max())

    contour_levels = np.mgrid[energy_range[0]:energy_range[1]:contour_levels]

    # Plot the landscape
    plt.figure()
    plt.pcolor(total_angles, total_radii, total_energy, vmin=energy_range[0],
               vmax=energy_range[1], cmap='viridis')
    cbar = plt.colorbar()
    cbar.set_label('Energy (eV)', size='x-large')
    cs = plt.contour(total_angles, total_radii, total_energy, colors='black',
                     levels=contour_levels, linewidths=0.6)

    for angle in facet_angles:
        plt.plot([angle, angle], [total_radii.min(), total_radii.max()],
                 color='k', linestyle='-', linewidth=1)
    xlabel = []
    for i in range(len(facet_angles)):
        xlabel.append('$\Omega_' + str(i+1) + '$')

    plt.xlabel('Angle', size='large')
    plt.ylabel('$r$ ($\mathrm{\AA}$)', size='x-large', fontname='Georgia')
    plt.xticks(facet_angles, xlabel, size='x-large')
    plt.clabel(cs, fontsize=10, inline_spacing=15, fmt='%1.1f', manual=True)
    plt.show()


###########
# METHODS #
###########

def reconstruct_chain(chain, verbose=False):

    if verbose:
        print("Analyzing facets in chain...")

    # Find all the facets in the links

    chain_facets = []
    end_facets = []

    for data in chain.keys():
        for facet in chain[data]['link']:

            if facet not in chain_facets:
                chain_facets.append(facet)

            if facet not in end_facets:
                end_facets.append(facet)
            else:
                end_facets.remove(facet)

    if verbose:
        print("Found " + str(len(chain_facets)) + " facets, of which "
              + str(len(end_facets)) + " are termination facets.")

    # Check if the chain has two termination facets.
    if len(end_facets) != 2:
        raise ValueError('Edges are not connected in a chain. Aborting...')
    # TODO Handle case of circular connection of paths

    # Choose the starting facet from the termination facets
    facet_chain = [end_facets[START_FACET], ]

    chain_links = []
    chain_landscapes = []

    if verbose:
        print("Constructing chain...")
        print('')

    while len(facet_chain) < len(chain_facets):

        last_facet = facet_chain[-1]

        if verbose:
            print('Last Facet:')
            print(str(last_facet))

        # Find the link to the last facet and is not already in the path chain
        for data in chain.keys():
            if chain[data]['link'] not in chain_links \
                and last_facet in chain[data]['link']:

                # Get the other facet in the link
                other_facet = chain[data]['link'].copy()
                other_facet.remove(last_facet)
                other_facet = other_facet[0]

                if verbose:
                    print('Other Facet:')
                    print(str(other_facet))

                facet_chain.append(other_facet)

                # Make sure that the landscape edge is set up in the right
                # direction.
                if chain[data]['link'][0] != last_facet:

                    if verbose:
                        print('Flipped path (before):')
                        print(str(chain[data]['link'][0]))
                        print(str(chain[data]['link'][1]))

                    chain[data]['link'].reverse()

                    if verbose:
                        print('Flipped path (after):')
                        print(str(chain[data]['link'][0]))
                        print(str(chain[data]['link'][1]))

                    chain[data]['landscape'].flip_coordinates('Angle')

                chain_landscapes.append(chain[data]['landscape'])
                chain_links.append(chain[data]['link'])

        if verbose:
            print('')

    if verbose:
        print('----------------')
        print('Facet Chain:')
        for facet in facet_chain:
            print(str(facet))
        print('----------------')
        print('Path Chain:')
        for link in chain_links:
            print('From')
            print(str(link[0]))
            print('To')
            print(str(link[1]))

    return (chain_landscapes, facet_chain)
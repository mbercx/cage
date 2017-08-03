# Encoding: utf-8

import os
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

from cage.landscape import LandscapeAnalyzer
from cage.core import Facet

"""
Analysis tools for the data produced by calculation set up using the cage setup
tools. Requires expansion.

"""

START_FACET = 0  # 0 or 1 -> determines facet to start chain from

#TODO This code needs serious improvement. It's damn near illegible


def landscape_analysis(lands_dir, cation, energy_range, interp_mesh,
                       contour_levels, verbose):

    lands_dir = os.path.abspath(lands_dir)

    dir_list = [os.path.abspath(os.path.join(lands_dir, file))
                for file in os.listdir(lands_dir)
                if os.path.isdir(os.path.join(lands_dir, file))]

    if verbose:
        print("Extracting Landscape data from " + lands_dir + "...")
        for directory in dir_list:
            print(directory)

    # TODO Make this part ignore directories that do not have the right files
    edges = [[LandscapeAnalyzer.from_file(os.path.join(directory,
                                                       'landscape.json')),
              [Facet.from_file(os.path.join(directory, 'init_facet.json')),
              Facet.from_file(os.path.join(directory, 'final_facet.json'))]]
             for directory in dir_list]

    if verbose:
        print("Found " + str(len(edges)) + " directories with landscape data.")

    for edge in edges:
        edge[0].analyze_cation_energies(facet=edge[1][0], cation=cation)

    # TODO Make landscapes and paths connected via dictionary

    if verbose:
        print("Analyzing facets...")

    # Find all the facets in the paths
    facets = []
    for edge in edges:
        for facet in edge[1]:
            if facet not in facets:
                facets.append(facet)

    if verbose:
        print("Found " + str(len(facets)) + " facets.")

    if verbose:
        print("Looking for end facets...")

    # Find end_facets
    end_facets = []
    for edge in edges:
        for facet in edge[1]:
            if facet not in end_facets:
                end_facets.append(facet)
            else:
                end_facets.remove(facet)

    # Check if the chain is connected
    if len(end_facets) != 2:
        raise ValueError('Edges are not connected in a chain. Aborting...')
    # TODO Handle case of circular connection of paths

    # Find the starting path
    facet_chain = [end_facets[START_FACET], ]
    path_chain = []
    landscape_chain = []

    if verbose:
        print("Constructing chain...")
        print('')

    while len(facet_chain) < len(facets):

        last_facet = facet_chain[-1]

        if verbose:
            print('Last Facet:')
            print(str(last_facet))

        # Find the path that connects to the last facet and is not already in
        # the path chain
        for edge in edges:
            if edge[1] not in path_chain and last_facet in edge[1]:
                other_facet = edge[1].copy()
                other_facet.remove(last_facet)
                other_facet = other_facet[0]

                if verbose:
                    print('Other Facet:')
                    print(str(other_facet))

                facet_chain.append(other_facet)

                # Make sure that the landscape edge is set up in the right
                # direction.
                if edge[1][0] != last_facet:

                    if verbose:
                        print('Flipped path (before):')
                        print(str(edge[1][0]))
                        print(str(edge[1][1]))

                    edge[1].reverse()

                    if verbose:
                        print('Flipped path (after):')
                        print(str(edge[1][0]))
                        print(str(edge[1][1]))

                    edge[0].flip_coordinates('Angle')

                landscape_chain.append(edge[0])
                path_chain.append(edge[1])

        if verbose:
            print('')

    if verbose:
        print('----------------')
        print('Facet Chain:')
        for facet in facet_chain:
            print(str(facet))
        print('----------------')
        print('Path Chain:')
        for edge in path_chain:
            print('From')
            print(str(edge[0]))
            print('To')
            print(str(edge[1]))

    # Interpolate the landscapes to a uniform mesh

    if verbose:
        print('-----------')
        print("Finding proper radii and angles...")
    # Find the proper radii
    min_max_radius = 1e6
    max_min_radius = 0
    for lands in landscape_chain:
        rmax = lands.datapoints['Distance'].max()
        if rmax < min_max_radius:
            min_max_radius = rmax
        rmin = lands.datapoints['Distance'].min()
        if rmin > max_min_radius:
            max_min_radius = rmin

    if verbose:
        print("")
        print('Largest minimal radius = ' + str(max_min_radius))
        print('Smallest maximal radius = ' + str(min_max_radius))

    # Adjust the angles to make one angle coordinate for all edges
    facet_angles = [0, landscape_chain[0].datapoints['Angle'].max()]

    for lands in landscape_chain[1:]:

        if verbose:
            print('Maximum angle = ' + str(facet_angles[-1]))

        lands.datapoints['Angle'] += facet_angles[-1]
        facet_angles.append(lands.datapoints['Angle'].max())

    all_radii = []
    all_angles = []
    all_energy = []

    # Interpolate the landscapes
    for lands in landscape_chain:

        data = lands.datapoints

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
        radii = data['Distance'].reshape(nradii, nangles)  # [::nradii]
        angles = data['Angle'].reshape(nradii, nangles)  # [:nangles]
        energy = data['Energy'].reshape(nradii, nangles)

        if verbose:
            print('Shape angles: ' + str(angles.shape))
            print('Shape radii: ' + str(radii.shape))
            print('Shape energy: ' + str(energy.shape))

        new_angles, new_radii = np.mgrid[
                                angles.min():angles.max():interp_mesh[0],
                                max_min_radius:min_max_radius:interp_mesh[1]
                                ]

        if verbose:
            print('-------------')
            print('Shape new_angles: ' + str(new_angles.shape))
            print('Shape new_radii: ' + str(new_radii.shape))

        tck = interpolate.bisplrep(angles, radii, energy, s=0.01)

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

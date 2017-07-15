# Encoding: utf-8

import sys
import os
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from pylab import *

from cage.landscape import LandscapeAnalyzer
from cage.facetsym import Facet

"""
Plot the landscape of connected edges. This is some damn illegible code, but
improving it would require thinking, and aint nobody got time for that.
"""

# PARAMETERS TODO DIRTY AS FUCK -> Find a better way
CATION='Na'
R_CUTOFF = 1.1
E_MAX_ADJUSTMENT = -0.5
START_FACET = 0 # 0 or 1 -> determines facet to start chain from
ANGLE_MESH = 0.003
RADIUS_MESH = 0.01
CONTOUR_LEVELS = np.mgrid[0.1:2:0.1]

# Input

dirname = sys.argv[1]
dir_list = [os.path.abspath(dir) for dir in os.listdir(dirname)
            if os.path.isdir(dir)]

edges = [[LandscapeAnalyzer.from_file(os.path.join(dir, 'landscape.json')),
          [Facet.from_file(os.path.join(dir, 'init_facet.json')),
          Facet.from_file(os.path.join(dir, 'final_facet.json'))]]
         for dir in dir_list]

for edge in edges:
    edge[0].analyze_cation_energies(facet=edge[1][0],cation=CATION)

# TODO Make landscapes and paths connected via dictionary

# Find all the facets in the paths
facets = []
for edge in edges:
    for facet in edge[1]:
        if facet not in facets:
            facets.append(facet)

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

while len(facet_chain) < len(facets):
    last_facet = facet_chain[-1]
    print('Last Facet:')
    print(str(last_facet))
    # Find the path that connects to the last facet and is not already in the
    # path chain
    for edge in edges:
        if edge[1] not in path_chain and last_facet in edge[1]:
            other_facet = edge[1].copy()
            other_facet.remove(last_facet)
            other_facet = other_facet[0]
            print('Other Facet:')
            print(str(other_facet))
            facet_chain.append(other_facet)

            if edge[1][0] != last_facet:
                print('Flipped path (before):')
                print(str(edge[1][0]))
                print(str(edge[1][1]))
                edge[1].reverse()
                print('Flipped path (after):')
                print(str(edge[1][0]))
                print(str(edge[1][1]))
                edge[0].flip_coordinates('Angle')

            landscape_chain.append(edge[0])
            path_chain.append(edge[1])

    print('')
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

# Find the proper radii
min_max_radius = 1e6
max_min_radius = 0
for landscape in landscape_chain:
    rmax = landscape.datapoints['Distance'].max()
    if rmax < min_max_radius:
        min_max_radius = rmax
    rmin = landscape.datapoints['Distance'].min()
    if rmin > max_min_radius:
        max_min_radius = rmin

max_min_radius += R_CUTOFF

print('-----------')
print('Largest minimal radius = ' + str(max_min_radius))
print('Smallest maximal radius = ' + str(min_max_radius))


# Adjust the angles to make one angle coordinate for all edges
facet_angles = [0, landscape_chain[0].datapoints['Angle'].max()]

for landscape in landscape_chain[1:]:
    print('Maximum angle = ' + str(facet_angles[-1]))
    landscape.datapoints['Angle'] += facet_angles[-1]
    facet_angles.append(landscape.datapoints['Angle'].max())

all_radii = []
all_angles = []
all_energy = []

# Interpolate the landscapes
for landscape in landscape_chain:

    data = landscape.datapoints

    data['Distance'] = np.round(data['Distance'], 5)
    data = np.sort(data, order=['Distance', 'Angle'])

    # Find the number of radii and angles
    r_init = data['Distance'][0]
    nangles = 1
    while abs(data['Distance'][nangles] - r_init) < 1e-5:
        nangles += 1
    nradii = int(len(data) / nangles)
    print('')
    print('-----------')
    print('Number of Angles = ' + str(nangles))
    print('Number of Radii = ' + str(nradii))

    # Get the right format for the data
    radii = data['Distance'].reshape(nradii, nangles)  # [::nradii]
    angles = data['Angle'].reshape(nradii, nangles)  # [:nangles]
    energy = data['Energy'].reshape(nradii, nangles)

    print('Shape angles: ' + str(angles.shape))
    print('Shape radii: ' + str(radii.shape))
    print('Shape energy: ' + str(energy.shape))

    new_angles, new_radii = np.mgrid[ angles.min():angles.max():ANGLE_MESH,
                                    max_min_radius:min_max_radius:RADIUS_MESH ]
    print('-------------')
    print('Shape new_angles: ' + str(new_angles.shape))
    print('Shape new_radii: ' + str(new_radii.shape))
    tck = interpolate.bisplrep(angles, radii, energy, s=0.01)

    new_energy = interpolate.bisplev(new_angles[:,0], new_radii[0,:], tck)

    all_radii.append(new_radii)
    all_angles.append(new_angles)
    all_energy.append(new_energy)

total_radii = np.concatenate(tuple(all_radii))
total_angles = np.concatenate(tuple(all_angles))
total_energy = np.concatenate(tuple(all_energy))
total_energy -= total_energy.min()

plt.figure()
plt.pcolor(total_angles, total_radii, total_energy, vmin=total_energy.min(),
           vmax=total_energy.mean() - E_MAX_ADJUSTMENT)
cbar = plt.colorbar()
cbar.set_label('Energy (eV)', size='large')
CS = plt.contour(total_angles, total_radii, total_energy, colors='black',
            levels=CONTOUR_LEVELS)
for angle in facet_angles:
    plt.plot([angle, angle], [total_radii.min(), total_radii.max()], color='k',
         linestyle='-', linewidth=1)
xlabel = []
for i in range(len(facet_angles)):
    xlabel.append('$\Omega_' + str(i+1) + '$')
plt.xlabel('Angle', size='large')
plt.ylabel('Radius', size='large')
plt.xticks(facet_angles, xlabel, size='large')
plt.clabel(CS, fontsize=10, inline_spacing=15, fmt='%1.1f', manual=True)
plt.show()









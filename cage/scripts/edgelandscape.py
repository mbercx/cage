# Encoding: utf-8

import sys
import os
import itertools

from cage.landscape import LandscapeAnalyzer
from cage.facetsym import Facet

"""
Plot the landscape of connected edges.
"""

dirname = sys.argv[1]
dir_list = [os.path.abspath(dir) for dir in os.listdir(dirname)
            if os.path.isdir(dir)]


# landscapes = [LandscapeAnalyzer.from_file(os.path.join(dir, 'landscape.json'))
#               for dir in dir_list]

paths = [[Facet.from_file(os.path.join(dir, 'init_facet.json')),
          Facet.from_file(os.path.join(dir, 'final_facet.json'))]
         for dir in dir_list]

# TODO Make landscapes and paths connected

# Find all the facets in the paths
facets = []
for path in paths:
    for facet in path:
        if facet not in facets:
            facets.append(facet)

# Find end_facets
end_facets = []
for path in paths:
    for facet in path:
        if facet not in end_facets:
            end_facets.append(facet)
        else:
            end_facets.remove(facet)

if len(end_facets) != 2:
    raise ValueError('Edges are not connected in a chain. Aborting...')

# Find the starting path
facet_chain = [end_facets[0], ]
path_chain = []

while len(facet_chain) < len(facets):
    for path in paths:
        last_facet = facet_chain[-1]
        if last_facet in path:
            other_facet = path.copy()
            other_facet.remove(last_facet)
            other_facet = other_facet[0]
            facet_chain.append(other_facet)
            path_chain.append(path)

        if path_chain[-1][0] != last_facet:
            path_chain[-1].reverse()

print('Facet Chain:')
for facet in facet_chain:
    print(str(facet))

print('Path Chain:')
for path in path_chain:
    print('From')
    print(str(path[0]))
    print('To')
    print(str(path[1]))

links = []
for combination in itertools.combinations(paths, 2):
        if combination[0][0] == combination[1][0]:
            links.append([combination[0], combination[1], 'inverse'])
        if combination[0][1] == combination[1][1]:
            links.append([combination[0], combination[1], 'inverse'])
        elif combination[0][1] == combination[1][0]:
            links.append([combination[0], combination[1], 'straight'])
        elif combination[0][0] == combination[1][1]:
            links.append([combination[1], combination[0], 'straight'])

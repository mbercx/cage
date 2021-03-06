import os

import cage

os.chdir('/Users/bercx1/Documents/Projects/Systems/B12')

mol = cage.core.Cage.from_poscar('CONTCAR')

# Test the JSON properties of Facet

facet = mol.facets[0]

cage.core.Facet.from_dict(facet.as_dict())

facet.to(fmt="json",filename='test.json')

facet_read = cage.core.Facet.from_file('test.json')

print(facet == facet_read)

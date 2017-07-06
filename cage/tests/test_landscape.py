
# TODO Figure out how to use the actual test files

import cage.tests as tests
import os

test_dir = os.path.dirname(tests.__file__)
test_files_dir = os.path.join(test_dir, 'testfiles')

def Test_LandscapeAnalyzer():
    from cage.facetsym import Facet
    facet = Facet.from_file(os.path.join(test_files_dir, 'facetcalc',
                                         'facet.json'))
    from cage.landscape import LandscapeAnalyzer
    ls = LandscapeAnalyzer.from_data(os.path.join(test_files_dir, 'facetcalc'))

    return ls.plot_facet_landscape(facet)
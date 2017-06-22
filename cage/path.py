# Encoding: utf-8
# Written for python 3.6

import pymatgen as pmg
from pymatgen.core.structure import SiteCollection

"""
Tools for analyzing pathways on Cage molecules.
"""

class Path(SiteCollection):
    """
    A Path is a List of Sites in a structure, plain and simple.
    """
    def __init__(self, coordinates):
        """
        Initializes an instance of a Path from a List of coordinates.
        """
        pass

    def set_up_neb(self):
        """
        Sets up the NEB path for a NwChem calculation.
        :return:
        """
        pass
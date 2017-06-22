# encoding: utf-8
# Copyright (c) Uncle Sam

import pymatgen as pmg
import pymatgen.symmetry.analyzer as syman

from itertools import combinations
import numpy as np
from math import pi

"""
Analysis tools to find the non-equivalent faces of fullerene shaped molecules. Still very much a work in progress.
"""

__author__ = "Marnik Bercx"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__status__ = "nowhere"
__date__ = "14 JUN 2017"


class Cage(pmg.Molecule):
    """
    A Cage is a molecule that is sort of shaped as a fullerene, based of the Molecule class from pymatgen.
    """

    def __init__(self, species, coords, charge=0, spin_multiplicity=None, validate_proximity=False,
                 site_properties=None):
        """
        Creates a Mutable Molecule Cage, specifically designed to work on the facets etc.

        The Cage molecule is automatically centered on its center of mass.
        :param species:
        :param coords:
        """
        super(Cage, self).__init__(species, coords)
        self._center()
        self._pointgroup = None
        self._facets = self._find_surface_facets()

    def _center(self):
        """
        Center the Cage by updating the sites.
        :return:
        """

        # Find the new Cartesian coordinates
        center = self.center_of_mass
        new_coords = np.array(self.cart_coords) - center

        # Update the sites
        sites = []
        for i in range(len(self.species)):
            prop = None
            if self.site_properties:
                prop = {k: v[i] for k, v in self.site_properties.items()}
            sites.append(pmg.Site(self.species[i], new_coords[i], properties=prop))

        self._sites = sites

    def _find_surface_facets(self):
        """
        Finds all the surface facets of the Cage object.

        Currently does not expand the facets to 4 sites in case it finds other sites which are in the plane of the site,
        as defined by 3 site points.
        :return: List of Facet objects
        """
        sites_all = self.sites

        # Find all the non Hydrogen sites
        sites_nonH = []
        for site in sites_all:
            if not site.species_string == 'H':
                sites_nonH.append(site)

        # Find all of the facets from combinations of three Non-H Sites
        facets_all = []
        for combination in combinations(sites_nonH, 3):
            facets_all.append(Facet(list(combination)))

        # Flip the normal of the facets in case it points to the origin
        for facet in facets_all:
            if facet.angle_to_normal(self.center_of_mass) < pi/2:
                facet.flip_normal()

        # Find all the facets that are "on the surface"
        facets_surf = []
        for facet in facets_all:

            # Calculate the angles between all nonH sites and the surface normal
            site_angles = []
            for site in sites_nonH:
                site_angles.append(facet.angle_to_normal(site.coords))

            # If all the angles are larger than pi/2, it's a surface site
            all_angles_smaller = True
            for angle in site_angles:
                if angle + 0.01 < pi/2:
                    all_angles_smaller = False

            # In that case, add it to the surface sites
            if all_angles_smaller:
                facets_surf.append(facet)

        return facets_surf

    @property
    def facets(self):
        """
        Surface facets of the Cage.
        :return:
        """
        return self._facets

    @classmethod
    def from_poscar(cls, filename):
        """
        Imports a Cage from a VASP POSCAR file.
        :return: Cage
        """
        # Import the structure from the POSCAR file
        structure = pmg.Structure.from_file(filename)

        # Generate the molecule object from the structure sites
        cage = cls(structure.species, structure.cart_coords)

        return cage

    def to_poscar(self, filename='POSCAR'):
        """
        Writes the Cage to a POSCAR file.
        """
        pass

    def find_PointGroup(self):
        """
        Find the Schoenflies PointGroup of the Cage molecule.
        :return: PointGroup object
        """
        return syman.PointGroupAnalyzer(self).get_pointgroup()

    def find_noneq_facets(self):
        """
        Find all of the nonequivalent facets of the Cage.

        First try: just
        :return:
        """
        # Set up the point group analyzer
        pgan = syman.PointGroupAnalyzer(self)

        # Find the full set of symmetry operations
        symmops = syman.generate_full_symmops(pgan.symmops,0.01)

        # Find all the non-equivalent facets
        facets_noneq = []
        for facet in self._facets:
            nonequivalent = True
            # Check to see if the facet is equivalent to one in the
            # nonequivalent list
            for facet_noneq in facets_noneq:
                for symm in symmops:
                    symm_facet_center = symm.operate(facet.center.tolist())
                    if np.linalg.norm(symm_facet_center - facet_noneq.center)\
                            < 1e-2:
                        nonequivalent = False
            if nonequivalent == True:
                facets_noneq.append(facet)

        return facets_noneq

class Facet(object):
    """
    Facet of a Molecule object, defined by a list of Sites
    """

    def __init__(self, sites):
        self._sites = sites
        self._center = siteCenter(tuple(self.sites))
        self._normal = self._find_normal()

    def _find_normal(self):
        """
        Finds the normal vector of the surface, pointing away from the origin
        """
        normal = np.cross(self._sites[0].coords - self._sites[1].coords,
                          self._sites[0].coords - self._sites[2].coords)

        normal /= 2  # Make the length of the normal equal to the surface area

        return normal

    @property
    def sites(self):
        """
        The sites that define the Facet.
        :return: List of Sites
        """
        return self._sites

    @property
    def center(self):
        """
        The center of the Facet.
        :return: Array of the center coordinates
        """
        return self._center

    @property
    def normal(self):
        """
        Surface normal of the Facet
        :return: Array of the normal vector
        """
        return self._normal

    def __str__(self):
        output = ['Facet with sites:']
        for site in self.sites:
            output.append(site.__str__())

        return '\n'.join(output)

    def flip_normal(self):
        """
        Flip the direction of the surface normal of the Facet.
        :return:
        """
        self._normal = -self._normal

    def angle_to_normal(self, coordinate):
        """
        Find the angle between the vector that connects the center and the coordinate and the normal.
        :param Site: 1x3 Array
        :return: angle value (in radians)
        """
        return angle_between(coordinate - self._center, self._normal)


# Utility functions that may be useful across classes

def siteCenter(sites):
    """
    Find the center of a collection of sites. Not particularly fast nor clever.
    :param sites: Tuple of Site objects
    :return: Array of the cartesian coordinates of the center of the sites
    """
    assert isinstance(sites, tuple)
    nSites = len(sites)

    # Sum the x,y, and z values of the site coordinates
    sum = np.zeros(3)
    for site in sites:
        sum[0] += site.x
        sum[1] += site.y
        sum[2] += site.z

    return sum / nSites


def schoenflies_to_HM():
    """
    Function for converting the Schoenflies point group symbol to the Hermann Mangiun one.
    :return:
    """
    pass


# Functions stolen from SO

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """
    Returns the angle in radians between vectors 'v1' and 'v2'::
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


# Other functions, might be deprecated

def molecule_from_poscar(filename):
    """
    Import a molecule object from a VASP POSCAR file and centers it at the
    center of mass.
    :param filename: Filename of the POSCAR file
    :return: Generated molecule
    """

    # Import the structure from the POSCAR file
    structure = pmg.Structure.from_file(filename)

    # Generate the molecule object from the structure sites
    mol = pmg.Molecule.from_sites(structure.sites)

    # Center molecule on center of mass
    mol = mol.get_centered_molecule()

    return mol


def find_all_symmops(mol):
    """
    Find all the symmetry operations of a molecule.
    :param mol: Molecule object
    :return: List of symmetry operations
    """

    # Find the point group of the molecule
    pgan = syman.PointGroupAnalyzer(mol)

    # Find all the symmetry operations of point group
    return syman.generate_full_symmops(pgan.symmops, 0.1)
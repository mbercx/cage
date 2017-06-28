# encoding: utf-8
# Copyright (c) Uncle Sam

import pymatgen as pmg
import pymatgen.symmetry.analyzer as syman
import numpy as np
import json

from pymatgen.core.structure import SiteCollection
from pymatgen.core.structure import Molecule
from monty.json import MSONable
from itertools import combinations
from math import pi
from monty.io import zopen

"""
Analysis tools to find the non-equivalent faces of fullerene shaped molecules.
Still very much a work in progress.
"""

__author__ = "Marnik Bercx"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__status__ = "initialization"
__date__ = "14 JUN 2017"


class Cage(Molecule):
    """
    A Cage is a molecule that is sort of shaped as a fullerene, based of the
    Molecule class from pymatgen.
    """

    def __init__(self, species, coords, charge=0, spin_multiplicity=None,
                 validate_proximity=False, site_properties=None):
        """
        Creates a Mutable Molecule Cage, specifically designed to work on the
        facets etc.

        The Cage molecule is automatically centered on its center of mass.
        :param species:
        :param coords:
        """
        super(Cage, self).__init__(species, coords, charge, spin_multiplicity,
                                   validate_proximity, site_properties)
        self._center()
        self._facets = self._find_surface_facets()
        self._pointgroup = None
        self._symmops = None

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

    def _find_surface_facets(self, ignore=(pmg.Element('H'),
                                           pmg.Element('Li'))):
        """
        Finds all the surface facets of the Cage object.

        Currently does not expand the facets to 4 sites in case it finds other
        sites which are in the plane of the site, as defined by 3 site points.

        :return: List of Facet objects
        """

        # Find all the sites which should not be ignored
        sites_valid = []
        for site in self.sites:
            if not site.specie in ignore:
                sites_valid.append(site)

        # Find all of the Facets from combinations of three valid Sites
        facets_all = []
        for combination in combinations(sites_valid, 3):
            facets_all.append(Facet(list(combination)))

        # Flip the normal of the facets in case it points to the origin
        for facet in facets_all:
            if facet.angle_to_normal(self.center_of_mass) < pi/2:
                facet.flip_normal()

        # Find all the facets that are "on the surface"
        facets_surf = []
        for facet in facets_all:

            # Calculate the angles between all valid sites and the surface normal
            site_angles = []
            for site in sites_valid:
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
    
    @property
    def symmops(self):
        """
        The symmetry operations of the Cage.
        :return: List of 
        """
        if not self._symmops:
            # Set up the point group analyzer
            pgan = syman.PointGroupAnalyzer(self)
    
            # Find the full set of symmetry operations
            self._symmops = syman.generate_full_symmops(pgan.symmops, 0.01)
        
        return self._symmops
        
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

    @classmethod
    def from_molecule(cls, mol):
        """
        Initializes a Cage from a Molecule.
        :param mol:
        :return:
        """
        assert type(mol) is Molecule
        return cls(species=mol.species, coords=mol.cart_coords,
                   charge=mol.charge, spin_multiplicity=mol.spin_multiplicity,
                   site_properties=mol.site_properties)
    
    def append(self, species, coords, validate_proximity=True,
               properties=None):
        """
        Overwrite the append method of the Molecule class, in order to
        remove the symmetry operations after the site has been appended.
        :param species:
        :param coords: 
        :param validate_proximity: 
        :param properties: 
        :return: 
        """
        super(Cage, self).append(species, coords, validate_proximity,
               properties)
        self._symmops = None

        # TODO Add parameter descriptions.

    def to_poscar(self, filename='POSCAR'):
        """
        Writes the Cage to a POSCAR file.
        """
        pass #TODO

    def find_PointGroup(self):
        """
        Find the Schoenflies PointGroup of the Cage molecule.
        :return: PointGroup object
        """
        return syman.PointGroupAnalyzer(self).get_pointgroup()

    def find_noneq_facets(self):
        """
        Find all of the nonequivalent facets of the Cage.

        :return: List of Facets
        """

        # Find all the non-equivalent facets
        facets_noneq = []
        for facet in self.facets:
            nonequivalent = True
            # Check to see if the facet is equivalent to one in the
            # nonequivalent list
            for facet_noneq in facets_noneq:
                for symm in self.symmops:
                    symm_facet_center = symm.operate(facet.center.tolist())
                    if np.linalg.norm(symm_facet_center - facet_noneq.center)\
                            < 1e-2:
                        nonequivalent = False
            if nonequivalent == True:
                facets_noneq.append(facet)

        return facets_noneq

    def find_facet_paths(self):
        """
        Find the non equivalent pathways between connected facets of the
        cage molecule. The facets can be connected by sharing an edge,
        or a vertex.
        :return:
        """

        # Find all paths, i.e. possible combinations of surface facets
        paths = list(combinations(self.facets,2))

        # Find the paths that share a vertex (this automatically finds those
        # that share an edge as well).
        vertex_sharing_paths = []
        for path in paths:
            cross_section = set(path[0].sites) & set(path[1].sites)
            if cross_section:
                vertex_sharing_paths.append(path)

        # Find the vertex sharing paths that are non equivalent
        non_eq_paths = []
        for path in vertex_sharing_paths:
            nonequivalent = True
            # Check to see if the path is equivalent with a path in the List of
            # non-equivalent paths
            for non_eq_path in non_eq_paths:
                for symm in self.symmops:
                    path_center = (path[0].center + path[1].center)/2
                    non_eq_path_center = (non_eq_path[0].center +
                                          non_eq_path[1].center)/2
                    symm_path_center = symm.operate(path_center)
                    if np.linalg.norm(symm_path_center - non_eq_path_center) \
                        < 1e-2:
                        nonequivalent = False
            if nonequivalent:
                non_eq_paths.append(path)

        return non_eq_paths


class Facet(SiteCollection, MSONable):
    """
    Facet of a Molecule object, defined by a list of Sites
    """

    #TODO Allow the facet to contain more than three sites
    #TODO Write those docstrings you lazy bastard!

    def __init__(self, sites, normal=None):
        self._sites = sites
        self._center = siteCenter(tuple(self.sites))
        if normal is not None:
            self._normal = normal # TODO Check if normal makes sense
        else:
            self._normal = self._find_normal()

    def __eq__(self, other):
        """
        Check if two Facets are the same.
        :param other:
        :return:
        """
        if (len(set(self.sites) & set(other.sites)) == 3) and \
                (self.normal == other.normal).all():
            return True
        else:
            return False



    @classmethod
    def from_str(cls, input_string, fmt="json"):
        """

        :param input_string:
        :param fmt:
        :return:
        """
        if fmt == "json":
            d = json.loads(input_string)
            return cls.from_dict(d)
        else:
            raise NotImplementedError('Only json formats have been '
                                      'implemented.')

    @classmethod
    def from_file(cls, filename):
        """

        :param filename:
        :return:
        """
        with zopen(filename) as file:
            contents = file.read()

        return cls.from_str(contents)

    @classmethod
    def from_dict(cls, d):
        """

        :param d:
        :return:
        """
        sites = []
        for site in d['sites']:
            sites.append(pmg.Site.from_dict(site))

        return cls(sites, normal=np.array(d['normal']))

    def to(self, fmt="json", filename=None):
        """

        :param fmt:
        :param filename:
        :return:
        """
        if fmt == "json":
            if filename:
                with zopen(filename, "wt", encoding='utf8') as file:
                    return json.dump(self.as_dict(), file)
            else:
                return json.dumps(self.as_dict())
        else:
            raise NotImplementedError("Currently only json format is "
                                      "supported.")

    def as_dict(self):
        """

        :return:
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "sites": []}
        for site in self:
            site_dict = site.as_dict()
            del site_dict["@module"]
            del site_dict["@class"]
            d["sites"].append(site_dict)
        d['normal'] = self.normal.tolist()
        return d


    def get_distance(self, i, j):
        """

        :param i:
        :param j:
        :return:
        """
        pass #TODO

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
    nsites = len(sites)

    # Sum the x,y, and z values of the site coordinates
    sum = np.zeros(3)
    for site in sites:
        sum[0] += site.x
        sum[1] += site.y
        sum[2] += site.z

    return sum / nsites


def schoenflies_to_HM():
    """
    Function for converting the Schoenflies point group symbol to the Hermann
    Mangiun one.
    :return:
    """
    pass  # TODO


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

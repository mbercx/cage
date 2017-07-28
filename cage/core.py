# encoding: utf-8
# Copyright (c) Uncle Sam

import json
import math

import numpy as np
import pymatgen as pmg
import cage.utils as utils
import pymatgen.symmetry.analyzer as syman

from itertools import combinations
from monty.io import zopen
from monty.json import MSONable
from pymatgen.core.structure import Molecule
from pymatgen.core.structure import SiteCollection
from pymatgen.core.operations import SymmOp

"""
Core tools of the cage package. Defines the Cage, OccupiedCage and Facet class.
"""

__author__ = "Marnik Bercx"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__status__ = "alpha"
__date__ = "14 JUN 2017"


SYMMETRY_TOLERANCE = 1e-2
# This is a tolerance value to determine the symmetry operations of the Cage.
# It is also used to determine which facets are equivalent. The standard value
# of 1E-2 is usually pretty good. In case the right non-equivalent facets are
# not found, it might be worth to try tweaking this value.


class Cage(Molecule):
    """
    A Cage is a pymatgen.Molecule-based object for molecules shaped similar
    to fullerenes.
    """

    def __init__(self, species, coords, charge=0, spin_multiplicity=None,
                 validate_proximity=False, site_properties=None):
        """
        Create a Cage instance. The Cage molecule's geometric center is
        automatically centered on the origin.

        Args:
            species (List of pymatgen.Specie): List of atomic species. Possible
                kinds of input include a list of dict of elements/species and
                occupancies, a List of elements/specie specified as actual
                Element/Specie, Strings ("Fe", "Fe2+") or atomic numbers
                (1,56).
            coords (List of (3,) numpy.ndarray): List of cartesian coordinates
                of each species.
            charge (float): Charge for the molecule. Defaults to 0.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 1 Ang apart. Defaults to False.
            site_properties (dict): Properties associated with the sites as
                a dict of sequences, e.g., {"magmom":[5,5,5,5]}. The
                sequences have to be the same length as the atomic species
                and fractional_coords. Defaults to None for no properties.

        Returns:
            (*cage.Cage*)
        """
        super(Cage, self).__init__(species, coords, charge, spin_multiplicity,
                                   validate_proximity, site_properties)
        self.center()
        self._facets = None
        self._pointgroup = None
        self._symmops = None
        self._facet_dict = None

    @property
    def facets(self):
        """
        Surface Facets of the Cage. Note that in case the surface facets have
        not been set up using *find.surface.facets()*, the property is equal to
        *None*.

        Returns:
            (List of Facets): The surface facets of the Cage, as set up using
                find_surface_facets()
        """
        return self._facets

    @property
    def pointgroup(self):
        """
        The Schoenflies PointGroup of the Cage molecule.

        Returns:
            (*pymatgen.symmetry.analyzer.PointGroup*)
        """
        if not self._pointgroup:
            self._pointgroup = syman.PointGroupAnalyzer(self).get_pointgroup()

        return self._pointgroup

    @property
    def symmops(self):
        """
        The symmetry operations of the Cage.

        Returns:
            (*List of pymatgen.Symmop*)
        """
        if not self._symmops:
            # Set up the point group analyzer
            pgan = syman.PointGroupAnalyzer(self)
    
            # Find the full set of symmetry operations
            self._symmops = syman.generate_full_symmops(pgan.symmops,
                                                        SYMMETRY_TOLERANCE)
        
        return self._symmops
        
    @classmethod
    def from_poscar(cls, filename):
        """
        Imports a Cage from a VASP POSCAR file.

        Args:
            filename (string): Filename of the POSCAR file.

        Returns:
            (*cage.Cage*)
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

        Args:
            mol (pymatgen.Molecule): The molecule from which to initialize the
                cage.

        Returns:
            (*cage.Cage*)
        """
        assert type(mol) is Molecule
        return cls(species=mol.species, coords=mol.cart_coords,
                   charge=mol.charge, spin_multiplicity=mol.spin_multiplicity,
                   site_properties=mol.site_properties)

    def to_poscar(self, filename='POSCAR'):
        """
        Writes the Cage to a POSCAR file.
        """
        pass  # TODO

    def center(self, point=None):
        """
        Center the Cage around a point by updating the sites, i.e. find the
        coordinates for the sites so that the geometric center is on the point
        provided. In case no point is provided, the molecule is centered around
        the origin.

        Args:
            point ((3,) numpy.ndarray): Point around which to center the
                molecule.
        """

        center = sum([site.coords for site in self.sites]) / len(self.sites)

        if point is not None:
            center -= point

        # Find the new coordinates
        new_coords = np.array(self.cart_coords) - center

        # Update the sites
        sites = []
        for i in range(len(self.species)):
            prop = None

            if self.site_properties:
                prop = {k: v[i] for k, v in self.site_properties.items()}

            sites.append(pmg.Site(self.species[i], new_coords[i],
                                  properties=prop))

        self._sites = sites

    def find_surface_facets(self, ignore=()):
        """
        Find all the surface facets of the Cage object. A surface facet is
        defined as a facet for which all atoms of non-ignored species are on
        one side of the surface defined by the facet.

        Currently does not expand the facets to 4 sites in case it finds other
        sites which are in the plane of the site, as defined by 3 site points.

        Args:
            ignore (Tuple of Elements/Species): The elements to ignore for the
                surface facet determination.
        """

        #TODO Expand the algorithm to more sites

        # Find all the sites which should not be ignored
        valid_sites = []
        for site in self.sites:
            if site.specie not in ignore:
                valid_sites.append(site)

        # Find all of the Facets from combinations of three valid Sites
        all_facets = []
        for combination in combinations(valid_sites, 3):
            all_facets.append(Facet(list(combination)))

        # Flip the normal of the facets in case it points to the center of mass
        # of the Cage
        for facet in all_facets:
            if facet.angle_to_normal(self.center_of_mass) < math.pi/2:
                facet.flip_normal()

        # Find all the facets that are "on the surface"
        facets_surf = []
        for facet in all_facets:

            # Calculate the angles between all valid sites and surface normal
            site_angles = [facet.angle_to_normal(site.coords)
                           for site in valid_sites]

            # If all the angles are larger than pi/2, it's a surface site
            all_angles_smaller = True
            for angle in site_angles:
                if angle + 0.01 < math.pi/2:
                    all_angles_smaller = False

            # In that case, add it to the surface sites
            if all_angles_smaller:
                facets_surf.append(facet)

        self._facets = facets_surf

    def redefine_origin(self, origin):
        """
        Change the coordinates of the Cage, in order to redefine the origin.

        Args:
            origin ((3,) numpy.ndarray): Origin coordinates.
        """
        # Find the new coordinates
        new_coords = np.array(self.cart_coords) - origin

        # Update the sites
        sites = []
        for i in range(len(self.species)):
            prop = None
            if self.site_properties:
                prop = {k: v[i] for k, v in self.site_properties.items()}
            sites.append(pmg.Site(self.species[i], new_coords[i],
                                  properties=prop))

        self._sites = sites

    def insert(self, index, species, coords, validate_proximity=True,
               properties=None):
        """
        Overwrite the append method of the Molecule class, in order to
        reset the symmetry operations and point group after the site has been
        appended.

        Args:
            index (int): Index to insert site.
            species (pymatgen.Specie): Species of inserted site.
            coords ((3,) numpy.ndarray): Coordinates of inserted site.
            validate_proximity (bool): Whether to check if inserted site is
                too close to an existing site. Defaults to True.
            properties (dict): A dictionary of properties for the Site.
        """
        super(Cage, self).append(species, coords, validate_proximity,
                                 properties)
        self._symmops = None
        self._pointgroup = None

    def append(self, species, coords, validate_proximity=True,
               properties=None):
        """
        Overwrite the append method of the Molecule class, in order to
        reset the symmetry operations and point group after the site has been
        appended.

        Args:
            species (pymatgen.Specie): Species of inserted site.
            coords ((3,) numpy.ndarray): Coordinates of inserted site.
            validate_proximity (bool): Whether to check if inserted site is
                too close to an existing site. Defaults to True.
            properties (dict): A dictionary of properties for the Site.
        """
        super(Cage, self).append(species, coords, validate_proximity,
                                 properties)
        self._symmops = None
        self._pointgroup = None

    def find_noneq_facets(self, tol=SYMMETRY_TOLERANCE):
        """
        Find the nonequivalent facets of the Cage.

        Args:
            tol (float): Tolerance for the equivalence condition, i.e. how much
                the distance between the centers is allowed to be after
                a symmetry operation.

        Returns:
            *List of Facets*: A set of non-equivalent facets of the Cage.
        """
        if not self.facets:
            print("Please set up surface facets first")
            return []

        # Find all the non-equivalent facets
        facets_noneq = []
        for facet in self.facets:
            facet_is_nonequivalent = True
            # Check to see if the facet is equivalent to one in the
            # nonequivalent list
            for facet_noneq in facets_noneq:
                for symm in self.symmops:
                    symm_facet_center = symm.operate(facet.center.tolist())
                    if np.linalg.norm(symm_facet_center - facet_noneq.center)\
                            < tol:
                        facet_is_nonequivalent = False
            if facet_is_nonequivalent:
                facets_noneq.append(facet)

        return facets_noneq

    def set_up_facet_list(self, fmt='dict', tol=SYMMETRY_TOLERANCE):
        """
        Set up a List of the surface facets, and how they relate to the
        non-equivalent facets, i.e. which non-equivalent facet they can be
        related to and using which symmetry operation:

        noneq_facet = symmop(facet)

        Args:
            fmt(str): Format of the Facet List. Can be either:

                *dict* -- A dictionary with the facets as keys and the tuple
                (noneq_facet, symmop) as the values.

                *str_array* -- A structured array of type:
                    [('surf_facet', Facet), ('noneq_facet', Facet),
                    ('symmop', SymmOp)]

            tol (float): Tolerance for the equivalence condition, i.e. how much
                the distance between the centers is allowed to be after
                a symmetry operation.

        Returns:
            (*List of cage.Facets*) -- The list of facets with their
                corresponding non-equivalent facet and symmetry operation.
                See the *fmt* argument.
        """
        if not self.facets:
            print("Please set up surface facets first")
            return []

        facet_list = []

        for facet in self.facets:
            for noneq_facet in self.find_noneq_facets():
                for symm in self.symmops:
                    symm_center = symm.operate(facet.center)
                    if np.linalg.norm(symm_center - noneq_facet.center) < tol:
                        facet_list.append((facet, noneq_facet, symm))
                        break

        if fmt == 'str_array':

            list_types = [('surf_facet', Facet), ('noneq_facet', Facet),
                          ('symmop', SymmOp)]

            facet_array = np.array(facet_list, dtype=list_types)

            if len(facet_array) == len(self.facets):
                return facet_array
            else:
                raise ValueError("Obtained array length is not equal to number"
                                 " of facets. Something must have gone wrong.")

        elif fmt == 'dict':

            facet_dict = {}

            for i in range(len(facet_list)):
                facet_dict[facet_list[i][0]] = (facet_list[i][1],
                                                facet_list[i][2])

            if sum([len(facet_list) for facet_list in facet_dict.values()])\
                    == len(self.facets):
                return facet_dict
            else:
                raise ValueError("Obtained number of facets in dict is not "
                                 "equal to number of surface facets. "
                                 "Something must have gone wrong.")

    def find_noneq_facet_chain(self, start=0, symm_tol=SYMMETRY_TOLERANCE,
                               verbose=False):
        """
        Find a chain of non equivalent facets, i.e. a collection of facets that
        are connected by edges. Automatically sorts the facets so they
        are connected by their neighbours in the list.

        Args:
            start (int): Determines from which termination facet the chain is
                constructed. This might be useful if the chain is not
                constructed as the user would like.
            symm_tol (float): Tolerance for the equivalence condition, i.e.
                how much the distance between the centers is allowed to be
                after a symmetry operation.
            verbose (bool): Print information about the analysis procedure.
                This is mainly useful when the result is not as expected.

        Returns:
            (*List of Facets*) -- A chain of facets connected by sharing an
                edge.

        """

        #TODO This code needs cleanup and serious testing

        if verbose:
            print("")
            print("Starting search for chain of non-equivalent facets in "
                  "molecule " + self.composition.__str__().replace(' ', '')
                  + "...")
            print("")
            print("Looking for non-equivalent facets...")

        facet_dict = self.set_up_facet_list('dict', tol=symm_tol)
        noneq_facets = self.find_noneq_facets(tol=symm_tol)

        if verbose:
            print("Found " + str(len(noneq_facets)) +
                  " non-equivalent facets.")
            print("")

        # Find the facets in the chain
        chain_facets = [noneq_facets[0]]
        chain_list_noneq_facets = [noneq_facets[0]]

        while len(chain_facets) < len(noneq_facets):

            new_chain_facet = False

            # Loop over the facets in the chain
            for chain_facet in chain_facets:

                # If a new facet has been appended, restart the loop
                if not new_chain_facet:

                    # Find a facet that shares an edge
                    for facet in self.facets:

                        # Check if the facet shares an edge and is not related
                        # to one of the non-equivalent facets in the chain
                        if len(set(facet.sites) & set(chain_facet.sites)) == 2\
                                and (facet_dict[facet][0] not in
                                     chain_list_noneq_facets):

                            chain_facets.append(facet)
                            chain_list_noneq_facets.append(
                                facet_dict[facet][0]
                            )
                            new_chain_facet = True
                            break

        if verbose:
            print("Found " + str(len(chain_facets)) +
                  " non-equivalent facets in chain.")
            print("")

        # Find the termination facets. These are defined as facets which only
        # have one edge with other chain facets.
        end_facets = []
        for facet in chain_facets:
            other_facets = chain_facets.copy()
            other_facets.remove(facet)
            for other_facet in other_facets:
                if len(set(facet.sites) & set(other_facet.sites)) == 2:
                    if facet in end_facets:
                        end_facets.remove(facet)
                        break
                    else:
                        end_facets.append(facet)

        if verbose:
            print("Found " + str(len(end_facets)) + " end facets in chain.")
            print("")

        if len(end_facets) == 0:
            raise ValueError("No termination facets found! Setting up chain "
                             "aborted.")

        # Sort the chain:
        try:
            facet_chain = [end_facets[start]]
        except IndexError:
            print("The requested starting index is too large. Taking the final"
                  " termination facet in the list.")
            facet_chain = [end_facets[-1]]

        other_facets = chain_facets.copy()
        other_facets.remove(facet_chain[0])

        for i in range(len(other_facets)):

            options = []

            for facet in other_facets:

                # See if the facet connects to the last facet in the chain
                if len(set(facet.sites) & set(facet_chain[-1].sites)) == 2:

                    # Check the amount of links this next facet has
                    leftover_facets = other_facets.copy()
                    leftover_facets.remove(facet)
                    number_links = 0
                    for leftover_facet in leftover_facets:
                        if len(set(facet.sites) & set(leftover_facet.sites)) \
                                                            == 2:
                            number_links += 1

                    options.append((facet, number_links))

            if len(options) == 1:
                facet_chain.append(options[0][0])
                other_facets.remove(options[0][0])
            else:
                for option in options:
                    if option[1] == 1:
                        facet_chain.append(option[0])
                        other_facets.remove(option[0])
                        break

        if len(facet_chain) < len(noneq_facets):
            print('WARNING: Could not connect all nonequivalent facets.')

        return facet_chain

    def find_facet_links(self, share_edge=False):
        """
        Find the non-equivalent links between facets of the cage
        molecule. The facets can be connected by sharing an edge, or a vertex.

        Args:
            share_edge (bool): Only return links between facets that
            share an edge.

        Returns:
            (List of (cage.Facet, cage.Facet) Tuples) - The
                non-equivalent facet links of the Cage.
        """

        # Find all links, i.e. possible combinations of surface facets
        links = list(combinations(self.facets, 2))

        # Find the links that share a vertex (this automatically finds
        # those that share an edge as well).
        vertex_sharing_links = []
        for link in links:
            cross_section = set(link[0].sites) & set(link[1].sites)
            if cross_section:

                # In case the user only wants edge-sharing paths, check that
                if share_edge:
                    if len(cross_section) == 2:
                        vertex_sharing_links.append(link)
                # Else just add the path to the list
                else:
                    vertex_sharing_links.append(link)

        # Find the vertex sharing paths that are non equivalent
        noneq_links = []
        for link in vertex_sharing_links:

            # Check to see if the path is equivalent with a path in the List of
            # non-equivalent paths
            nonequivalent = True
            for noneq_link in noneq_links:
                for symm in self.symmops:
                    link_center = (link[0].center + link[1].center)/2
                    noneq_link_center = sum((noneq_link[0].center,
                                                noneq_link[1].center))/2
                    symm_link_center = symm.operate(link_center)
                    connection_vector = symm_link_center - noneq_link_center
                    if np.linalg.norm(connection_vector) < 1e-2:
                        nonequivalent = False

            if nonequivalent:
                noneq_links.append(link)

        return noneq_links

    def find_noneq_chain_links(self, symm_tol=SYMMETRY_TOLERANCE,
                               verbose=False):
        """
        Find the links between the facets of the chain that connects a
        set of non equivalent facets.

        Args:
            symm_tol (float): Tolerance for the equivalence condition, i.e.
                how much the distance between the centers is allowed to be
                after a symmetry operation.
            verbose (bool): Print information about the analysis procedure.
                This is mainly useful when the result is not as expected.

        Returns:
            (*List of (cage.Facet, cage.Facet) Tuples*) --
                The links between the Facets in the chain of non-equivalent
                Facets.

        """

        facet_chain = self.find_noneq_facet_chain(symm_tol=symm_tol,
                                                  verbose=verbose)

        chain_links = []
        for index in range(len(facet_chain)-1):
            chain_links.append((facet_chain[index], facet_chain[index + 1]))

        return chain_links

    def find_farthest_facet(self, point):
        """
        Find the Facet of the Molecule that is the farthest away from the point
        provided.

        Args:
            point ((3,) numpy.ndarray): Point provided by user.

        Returns:
            (*cage.Facet*)
        """
        distance = 0
        furthest_facet = None

        for facet in self.facets:
            newdistance = np.linalg.norm(point - facet.center)
            if newdistance > distance:
                furthest_facet = facet
                distance = newdistance

        return furthest_facet
    
    def find_closest_facet(self, point):
        """
        Find the Facet of the Molecule that is the closest to the point 
        provided.

        Args:
            point ((3,) numpy.ndarray): Point provided by user.

        Returns:
            (*cage.Facet*)
        """
        distance = 1e6
        closest_facet = None

        for facet in self.facets:
            newdistance = np.linalg.norm(point - facet.center)
            if newdistance < distance:
                closest_facet = facet
                distance = newdistance

        return closest_facet


class OccupiedCage(Cage):
    """
    A Cage Molecule that has one or more cations docked on it.
    """

    CATIONS = (pmg.Element('Li'), pmg.Element('Na'), pmg.Element('K'),
               pmg.Element('Mg'))

    def __init__(self, species, coords, charge=0, spin_multiplicity=None,
                 validate_proximity=False, site_properties=None):
        """
        Initialize an OccupiedCage instance.

        Args:
            species (List of pymatgen.Specie): List of atomic species. Possible
                kinds of input include a list of dict of elements/species and
                occupancies, a List of elements/specie specified as actual
                Element/Specie, Strings ("Fe", "Fe2+") or atomic numbers
                (1,56).
            coords (List of (3,) numpy.ndarray): List of cartesian coordinates
                of each species.
            charge (float): Charge for the molecule. Defaults to 0.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 1 Ang apart. Defaults to False.
            site_properties (dict): Properties associated with the sites as
                a dict of sequences, e.g., {"magmom":[5,5,5,5]}. The
                sequences have to be the same length as the atomic species
                and fractional_coords. Defaults to None for no properties.

        Returns:
            (*cage.Cage*)
        """

        super(OccupiedCage, self).__init__(species, coords, charge,
                                           spin_multiplicity,
                                           validate_proximity,
                                           site_properties)
        self.center()
        self._docks = []

    @property
    def docks(self):
        return self._docks

    @property
    def facets(self):
        return self._facets

    def center(self, point=None):
        """
        Center the OccupiedCage around a point by updating the sites, i.e. find
        the coordinates for the sites so that the geometric center **of the
        anion** is moved to the point provided. In case no point is provided,
        the anion is centered around the origin.

        Args:
            point ((3,) numpy.ndarray): Point around which to center the
                molecule.
        """

        anion_coords = [site.coords for site in self.sites
                        if site.specie not in OccupiedCage.CATIONS]

        anion_center = sum(anion_coords)/len(anion_coords)

        super(OccupiedCage, self).center(anion_center)

    def add_dock(self, dock, cation='Li', docking_point=None):
        """
        Add a docking site to the OccupiedCage. If the chemical symbol of the
        cation is provided, the cation is appended to the OccupiedCage. In case
        the cation is equal to *None*, the cation is assumed to be present and

        Args:
            dock (cage.Facet): The Facet on which the cation is docked.
            cation (str): The chemical symbol of the cation element.
            docking_point ((3,) numpy.ndarray): Docking coordinates of the
                cation.
        """

        # Check if the dock is on of the facets in the OccupiedCage
        if dock not in self.facets:
            raise ValueError("Docking facet not found in the facet list of the"
                             " OccupiedCage.")

        if not cation:
            self._docks.append(dock)
        else:
            if docking_point:
                self.append(pmg.Element(cation), docking_point)
                self.set_charge_and_spin(self.charge,
                                         self.spin_multiplicity - 1)

            else:
                cation_coord = dock.center + 2 * dock.normal
                self.append(pmg.Element(cation), cation_coord)
                self.set_charge_and_spin(self.charge,
                                         self.spin_multiplicity - 1)
                self._docks.append(dock)

        self._facets = None

        # TODO Add some more checks

    @classmethod
    def from_cage_and_facets(cls, cage, facets, docking_points=(),
                             cation='Li'):
        """

        :param cage:
        :param facets:
        :param docking_points:
        :param cation:
        """
        occ_cage = cls(species=cage.species, coords=cage.cart_coords,
                       charge=cage.charge,
                       spin_multiplicity=cage.spin_multiplicity,
                       validate_proximity=True,
                       site_properties=cage.site_properties)

        # Add the docked cations to the Cage
        for index in range(len(facets)):
            try:
                occ_cage.add_dock(facets[index],
                                  docking_point=docking_points[index],
                                  cation=cation)
            except IndexError:
                occ_cage.add_dock(facets[index], cation=cation)

        return occ_cage

    @classmethod
    def from_poscar(cls, filename):
        """
        Initialize an OccupiedCage from a VASP POSCAR file.

        Args:
            filename:

        Returns:

        """
        pass #TODO

    @classmethod
    def from_file(cls, filename):
        """
        Initialize an OccupiedCage from a file.

        Args:
            filename:

        Returns:

        """
        pass #TODO

    def remove_surface_facet(self, facet):
        """

        :return:
        """
        surface_facets = self.facets
        if surface_facets:
            self._facets = surface_facets.remove(facet)
        else:
            print('Surface Facets have not been set up yet.')

    def find_surface_facets(self, ignore=None):
        """
        Find the surface facets of the OccupiedCage.
        :param ignore:
        :return:
        """
        mol = self.copy()
        anion = [site.coords for site in mol.sites
                 if site.specie not in OccupiedCage.CATIONS]

        mol.center(sum(anion)/len(anion))

        super(OccupiedCage, mol).find_surface_facets(ignore=ignore)

        surface_facets = mol.facets

        for dock in self.docks:
            surface_facets.remove(dock)

        self._facets = surface_facets


class Facet(SiteCollection, MSONable):
    """
    Facet of a Molecule object, defined by a list of Sites.
    """

    #TODO Allow the facet to contain more than three sites
    #TODO Write those docstrings you lazy bastard!

    def __init__(self, sites, normal=None):
        """
        Initialize a Facet with the provided sites.
        :param sites:
        :param normal:
        """
        self._sites = sites
        self._center = utils.site_center(tuple(self.sites))
        if normal is not None:
            self._normal = normal  # TODO Check if input normal makes sense
        else:
            self._normal = self._find_normal()

    def _find_normal(self):
        """
        Finds the normal vector of the surface, pointing away from the origin
        """
        normal = np.cross(self._sites[0].coords - self._sites[1].coords,
                          self._sites[0].coords - self._sites[2].coords)

        # Make length of the normal equal to 1
        normal = normal/np.linalg.norm(normal)

        return normal

    def __str__(self):
        output = ['Facet with sites:']
        for site in self.sites:
            output.append(site.__str__())

        return '\n'.join(output)

    def __eq__(self, other):
        """
        Check if two Facets are the same.
        :param other:
        :return:
        """
        if (len(set(self.sites) & set(other.sites)) == 3) and \
                np.allclose(self.normal, other.normal, atol=1e-3):
            return True
        else:
            return False

    def __hash__(self):
        """
        Make Facet hashable.
        :return:
        """
        return hash(str(self))

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

    def copy(self):
        return Facet(self._sites, self._normal)

    def get_normal_intersection(self, other, method='gonio'):
        """
        Find the intersection of the normal lines of the Facet and another one.

        Currently only works on an edge sharing Facet whose normal intersects
        with the normal of the facet.

        :return:
        """

        if method == 'gonio':

            edge = set(self.sites) & set(other.sites)
            if len(edge) != 2:
                raise ValueError('Provided facet does not share an edge. '
                                 'Goniometric method will not work.')

            edge_middle = sum([site.coords for site in edge]) / 2

            y1 = utils.distance(self.center, edge_middle)
            y2 = utils.distance(other.center, edge_middle)
            beta = utils.angle_between(self.normal, other.normal)
            psi = math.atan2(math.sin(beta), (y1 / y2 + math.cos(beta)))
            theta = beta - psi
            r = y1 / math.sin(theta)
            r1 = r * math.cos(theta)
            r2 = r * math.cos(psi)

            intersection1 = self.center - r1 * utils.unit_vector(self.normal)
            intersection2 = other.center - r2 * utils.unit_vector(other.normal)

            if utils.distance(intersection1, intersection2) < 1e-4:
                return intersection1
            else:
                print('Could not find perfect intersection. Returned closest '
                      'result.')
                return (intersection1 + intersection2)/2

        # Brute method. Currently abandoned
        #
        # elif method == 'brute':
        #     lines = []
        #     for facet in facets:
        #         line = Landscape.from_vertices([facet.center,
        #                                         facet.center - 5 *
        #                                         facet.normal],
        #                                        density=int(1e2))
        #         lines.append(line)
        #
        #     intersection = None
        #     distance = 1e3
        #     for point1 in lines[0].points:
        #         for point2 in lines[1].points:
        #             if utils.distance(point1, point2) < distance:
        #                 distance = utils.distance(point1, point2)
        #                 intersection = point1
        #     if distance < 1e-2:
        #         return intersection
        #     else:
        #         print('Brute force method failed.')
        #         return None

        else:
            NotImplementedError("Unknown method.")

    def redefine_origin(self, origin):
        """
        Change the coordinates of the Facet, in order to change the origin.
        :return:
        """
        # Find the new coordinates
        new_coords = np.array(self.cart_coords) - origin

        # Update the sites
        sites = []
        for i in range(len(self.species)):
            prop = None
            if self.site_properties:
                prop = {k: v[i] for k, v in self.site_properties.items()}
            sites.append(pmg.Site(self.species[i], new_coords[i],
                                  properties=prop))

        self._sites = sites

    def get_distance(self, i, j):
        """

        :param i:
        :param j:
        :return:
        """
        pass  #TODO

    def is_equivalent(self, other, symmops):
        """
        Check if a Facet is equivalent to another Facet based on a list of
        symmetry operations,
        :param other:
        :param symmops:
        :return:
        """
        is_equivalent = False
        for symm in symmops:
            symm_center = symm.operate(self.center)
            if np.linalg.norm(symm_center, other.center) < 1e-2:
                is_equivalent = True
        return is_equivalent

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

    def flip_normal(self):
        """
        Flip the direction of the surface normal of the Facet.
        :return:
        """
        self._normal = -self._normal

    def angle_to_normal(self, coordinate):
        """
        Find the angle between the vector that connects the center and the
        coordinate and the normal.
        :param coordinate:
        :return: angle value (in radians)
        """
        return utils.angle_between(coordinate - self._center, self._normal)
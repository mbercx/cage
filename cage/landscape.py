# encoding: utf-8
# Copyright (c) Uncle Sam

from monty.json import MSONable
from cage.facetsym import Cage

import math
import os

from itertools import combinations
import pymatgen as pmg
import pymatgen.io.nwchem as nw
import numpy as np
import cage.utils as utils

"""
Tools to set up and run calculations to study energy landscapes, specifically
for Cage molecules.
"""

# TODO Landscape should not be defined by vertices, but by points.
# TODO Instead, make a method 'From vertices' to initialize landscape

__author__ = "Marnik Bercx"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__status__ = "nowhere"
__date__ = "16 JUN 2017"

class Landscape(MSONable):
    """
    A line, area or volume that is to be studied.
    """
    def __init__(self, vertices, density=10):
        """
        Defines a landscape by defining the vertices (end points in the case of
        a line).

        :param vertices: A List of Array coordinates
        :param density: Float
        """
        self._vertices = vertices
        self._density = density
        self._dimension = self._find_dimension()
        self._points = self._set_up_landscape()

    def _find_dimension(self):
        """
        Derive the dimension of the Landscape from the vertices.

        WIP: Currently only designed to study lines.

        :return: Int
        """

        if len(self.vertices) == 1:
            raise IOError("Number of vertices must be at least two.")
        elif len(self.vertices) == 2:
            if np.linalg.norm(self.vertices[1] - self.vertices[0]) < 1e-3:
                raise IOError("User provided two vertices that are too close"
                              " to each other")
            else:
                return 1
        else:
            raise IOError("Higher dimensions than 2 not implemented yet.")

    def _set_up_landscape(self):
        """
        Finds all the points necessary to describe the Landscape with the
        desired density.
        :return: List of Array coordinates
        """
        if self.dimension == 1:
            vector = self.vertices[1] - self.vertices[0]
            unitvector = unit_vector(vector)
            length = np.linalg.norm(vector)
            npoints = math.ceil(length*self.density) + 1
            mesh_distance = length/(npoints)
            landscape_points = []
            for i in range(npoints):
                landscape_points.append(self.vertices[0]
                                        + i*mesh_distance*unitvector)
        else:
            raise IOError("2D and up not implemented yet.")

        return landscape_points

    @property
    def vertices(self):
        """
        The vertices of the landscape.

        :return: A List of Array coordinates
        """
        return self._vertices

    @property
    def density(self):
        """
        The density of the mesh, defined as number of points per Angstrom in
        each dimension.
        :return:
        """
        return self._density

    @property
    def dimension(self):
        """
        The dimension of the landscape.

        :return: Int
        """
        return self._dimension

    @property
    def points(self):
        """
        All of the points of the Landscape.
        :return:
        """
        return self._points

    def extend_by_rotation(self, axis, density=10):
        """
        Extends the landscape using a axis vector and turning all the vertices
        in the landscape around the origin by a value and direction determined
        by the axis vector.
        :param vector:
        :return:
        """

        # Find the rotation matrices
        rotation_matrices = []

        for i in range(int(density)):
            angle = np.linalg.norm(axis)*(i+1)/int(density)
            rotation_matrices.append(rotation_matrix(axis, angle))

        # Add all the points TODO This might be done more quickly
        points = self.points.copy()
        for matrix in rotation_matrices:
            for point in self.points:
                newpoint = point.dot(matrix)
                if not (newpoint == np.array([0, 0, 0])).all():
                    points.append(newpoint)

        self._points = points.copy()


class LandscapeAnalyzer(MSONable):
    """
    An analyzer class for interpreting data from calculations on a landscape.
    """
    def __init__(self, data, software="nwchem"):
        """
        Initialize an instance of LandscapeAnalyzer. This method will rarely
        be used directly. Usually a LandscapeAnalyzer is initialized from a
        directory or file.
        """
        self._data = data
        self._software = software

    @property
    def data(self):
        return self._data

    @property
    def software(self):
        return self._software

    @classmethod
    def from_data(cls, directory, output_file='result.out', software='nwchem'):
        """
        Looks through all subdirectories of 'directory' and looks for output
        files to extract data from to analyze.
        :return:
        """

        # TODO Make subdirectory finder recursive?

        if not os.path.isdir(directory):
            raise IOError("Directory " + directory + " not found.")

        # Find all the subdirectories in the specified directory
        dir_list = [d for d in os.listdir(directory)
                   if os.path.isdir(os.path.join(directory, d))]

        # Get all of the output of the calculations in the directory
        output = []
        if software == 'nwchem':
            for dir in dir_list:
                try:
                    out = nw.NwOutput(os.path.join(directory, dir,
                                                      output_file))
                except FileNotFoundError:
                    print('Did not find output in directory ' + dir)

                # Check if the output has an error
                if out.data[-1]['has_error']:
                    print('Error found in output in directory ' + dir)
                else:
                    output.append(out)
                    print('Grabbed output in directory ' + dir)

            data = [out.data[-1] for out in output]

        else:
            raise NotImplementedError("Only NwChem is currently supported.")

        return LandscapeAnalyzer(data)

    def plot_facet_landscape(self, coordinates="polar"):
        """
        Plot the landscape of a certain property on a facet.
        :return:
        """
        # TODO This method is horrendously specific, but I need to write it quick so I can show some results. Improve this later.

        # Find the initial facet
        cage_init = Cage.from_molecule(self.data[0]['molecules'][0])
        facet_init = cage_init.facets[0]
        for facet in cage_init.facets:
            if utils.distance(facet.center, cage_init.cart_coords[-1]) < \
                    utils.distance(facet_init.center,
                                   cage_init.cart_coords[-1]):
                facet_init = facet

        datapoints = []

        for data in self.data:

            if coordinates == "polar":

                Li_coord = [site.coords for site in data["molecules"][0].sites
                            if site.specie == pmg.Element('Li')][0]

                cage_init = Cage.from_molecule(data["molecules"][0])
                cage_final = Cage.from_molecule(data["molecules"][-1])

                r = np.linalg.norm(Li_coord)
                theta = angle_between(facet_init.center, Li_coord)
                if theta > math.pi/8:
                    theta = math.pi - theta
                coordinate = [r, theta]

            energy_initial = data['energies'][0]
            energy_final = data['energies'][-1]
            coordinate.append(energy_initial)

            datapoints.append(coordinate)

        return datapoints

        # TODO Finish, in case it makes sense


    def as_dict(self):
        """
        Return a dictionary representing the LandscapeAnalyzer instance.
        :return:
        """
        pass

    @classmethod
    def from_dict(cls, d):
        """
        Initialize the LandscapeAnalyzer from a dictionary.
        :param d:
        :return:
        """
        pass

    @classmethod
    def from_string(cls):
        """
        Initialize a LandscapeAnalyzer from a string.
        :return:
        """
        pass

    @classmethod
    def from_file(cls, fmt='json'):
        """
        Initialize an instance of LandscapeAnalyzer from a file.
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

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

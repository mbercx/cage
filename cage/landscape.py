# encoding: utf-8
# Copyright (c) Uncle Sam

from monty.json import MSONable

import math

from itertools import combinations

import numpy as np

"""
Tools to set up and run calculations to study energy landscapes, specifically
for Cage molecules.
"""

__author__ = "Marnik Bercx"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__status__ = "nowhere"
__date__ = "16 JUN 2017"

class Landscape(MSONable):
    """
    A line, area or volume that needs to be studied.
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
            raise IOError("Number of vertices must be at least one.")
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



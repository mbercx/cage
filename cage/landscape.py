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

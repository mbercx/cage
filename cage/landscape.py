# encoding: utf-8
# Copyright (c) Marnik Bercx

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
__status__ = "alpha"
__date__ = "16 JUN 2017"

class Landscape(MSONable):
    """
    A discrete mesh representing a line, area or volume.
    """
    def __init__(self, points):
        """
        Initializes a Landscape from the points that it consists of.

        :param points: (List of numpy.array)
        """
        self._points = points

    def __add__(self, other):
        """
        Add two Landscapes to each other into a new Landscape.

        :param other: (Landscape)
        :return: (Landscape)
        """
        points = self.points + other.points

        return Landscape(points)

    @property
    def points(self):
        """
        All of the points of the Landscape.

        :return:
        """
        return self._points

    @classmethod
    def from_vertices(cls, vertices, density=10):
        """
        Defines a landscape by defining the vertices (end points in the case of
        a line).

        :param vertices: (List of numpy.array)
        :param density: (float)
        """
        # Calculate the points of the Landscape depending on the number of
        # vertices
        if len(vertices) == 1:
            raise IOError("Number of vertices must be at least two.")
        elif len(vertices) == 2:
            if np.linalg.norm(vertices[1] - vertices[0]) < 1e-3:
                raise IOError("Vertices are too close to each other.")
            else:
                vector = vertices[1] - vertices[0]
                length = np.linalg.norm(vector)
                unitvector = vector / length
                npoints = int(math.ceil(length * density) + 1)
                mesh_distance = length / npoints
                points = []
                for i in range(npoints):
                    points.append(vertices[0] + i * mesh_distance * unitvector)
        else:
            raise IOError("Higher dimensions than 2 not implemented yet.")

        return Landscape(points)

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

    def as_dict(self):
        """
        A JSON serializable dictionary of the Landscape.
        :return: (dict)
        """
        dict = {"points": self.points}

        return dict


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

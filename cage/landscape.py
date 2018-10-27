# encoding: utf-8
# Copyright (c) Marnik Bercx

from monty.json import MSONable
from cage.core import Cage
from matplotlib.pyplot import subplots

import json
import math
import os

import pdb

import matplotlib.pyplot as plt
import pymatgen as pmg
import pymatgen.io.nwchem as nw
import numpy as np
import cage.utils as utils

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
    A selection of points representing a line, area or volume.
    """
    def __init__(self, points):
        """
        Initializes a Landscape from the points that it consists of.

        Args:
            points (list): List of 3x1 numpy.array coordinates of the points
            in the landscape.
        """

        if type(points) is list:
            self._points = points
        else:
            raise TypeError("Provided points are not formatted in a List.")

    def __add__(self, other):
        """
        Add two Landscapes to each other into a new Landscape.

        Args:
            other (cage.landscape.Landscape): Landscape to which the current
            Landscape should be added.

        Returns:
            Sum of the two landscapes, i.e. a landscape that consists of the
            points of the two landscapes combined.

        """

        points = self.points + other.points

        return Landscape(points)

    @property
    def points(self):
        """
        All of the points of the Landscape.

        Returns:
            (list) List of 3x1 numpy.array coordinates of the points
            in the landscape.

        """
        return self._points

    @property
    def center(self):
        """
        Center of the points of the Landscape.

        Returns:
            (3x1 numpy.array) Coordinates of the mathematical center of the
            points of the landscape.
        """
        return sum(self.points)/len(self.points)

    def change_center(self, center):
        """
        Change the center of the landscape. This simply translates all the
        coordinates so the new mathematical center is the one provided by
        the user.

        Args:
            center (3x1 numpy.array): New center of the landscape.

        """
        new_points = []

        for point in self.points:
            new_points.append(point - self.center + center)

        self._points = new_points

    def extend_by_rotation(self, axis, density=10.0, remove_endline=False,
                           distance_tol=1e-3):
        """
        Extends the landscape using an axis vector and turning all the vertices
        in the landscape around the origin by a value and direction determined
        by the axis vector.

        Args:
            axis (numpy.ndarray): Axis around which the landscape is
            rotated. The length of the vector determines the total rotation
            angle in radians.
            density (float): Density of grid points along the rotation
            angle. Defined in #grid points per radians.
            remove_endline (bool): Include the final rotation angle grid
            points.
            distance_tol (float): Minimum distance between two points in
            the landscape.

        Returns:

        """
        # TODO Extend this method
        # so it also allows rotations around other points than origin

        # Find the rotation matrices
        rotation_matrices = []
        total_angle = np.linalg.norm(axis)
        npoints = int(total_angle * density)

        if remove_endline:
            for i in range(npoints-1):
                angle = (i+1)/npoints*total_angle
                rotation_matrices.append(rotation_matrix(axis, angle))
        else:
            for i in range(npoints):
                angle = (i+1)/npoints*total_angle
                rotation_matrices.append(rotation_matrix(axis, angle))

        # Add all the points TODO This might be done more quickly
        points = self.points.copy()
        for matrix in rotation_matrices:
            for point in self.points:
                newpoint = point.dot(matrix)
                distance = np.linalg.norm(point - newpoint)
                if distance > distance_tol:
                    points.append(newpoint)

        self._points = points.copy()

    def extend_by_translation(self, vector, density=10):
        """
        Extends the Landscape by translating the points along a certain vector.

        Args:
            vector (3x1 numpy.array): Translation vector of the extension.
            density (float): Density of the grid.

        Returns:

        """
        pass #TODO

    def extend_from_point(self, point, extension, density=10):
        """
        Extends the Landscape by scaling the points from a specific point or
        origin, i.e. by homothetic transformations.

        Args:
            point:
            extension:
            density:

        Returns:

        """
        pass #TODO

    def as_dict(self):
        """
        A JSON serializable dictionary of the Landscape.

        Returns:

        """
        d = {"points": self.points}

        return d

    @classmethod
    def from_file(cls, filename, fmt="json"):
        """
        Load a Landscape from a file.

        Args:
            filename:
            fmt:

        Returns:

        """
        pass

    @classmethod
    def from_vertices(cls, vertices, num=10):
        """
        Define a landscape by providing the vertices (end points in the case of
        a line).

        Args:
            vertices:
            num:

        Returns:

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
                npoints = num  # int(length * num + 1)
                mesh_distance = length / npoints
                points = []
                for i in range(npoints):
                    points.append(vertices[0] + i * mesh_distance * unitvector)
        else:
            raise NotImplementedError("Higher dimensions than 1 not implemented yet.")

        return Landscape(points)

    @classmethod
    def create_sphere(cls, radius, center=np.array([0, 0, 0]),
                      axis=np.array([0, 0, 1]), density=15):
        """
        Set up a spherical landscape with a specified radius.

        Args:
            radius (float): Radius of the spherical landscape.
            center (3x1 numpy.array): Coordinates of the center of the
            spherical landscape.
            axis (3x1 numpy.array): Rotational axis which is used to
            construct the landscape.
            density (float): Grid density of the landscape, provided in
            number of grid points per radians.

        Returns:
            (Landscape): Spherical Landscape object as specified by the user.

        """
        starting_point = radius * unit_vector(axis)

        sphere_landscape = Landscape([starting_point, ])

        if angle_between(np.array([1, 0, 0]), axis) < 1e-2:
            phi_axis = unit_vector(
                np.cross(np.array([0, 1, 0]), starting_point)
            )
        else:
            phi_axis = unit_vector(
                np.cross(np.array([1, 0, 0]), starting_point)
            )

        sphere_landscape.extend_by_rotation(
            axis=phi_axis*math.pi,
            density=density)

        sphere_landscape.extend_by_rotation(
            axis=unit_vector(axis) * 2 * math.pi,
            density=density
        )

        sphere_landscape.change_center(center)

        return sphere_landscape


class LandscapeAnalyzer(MSONable):
    """
    An analyzer class for interpreting data from calculations on a landscape.
    """
    def __init__(self, data, datapoints=(), software="nwchem"):
        """
        Initialize an instance of LandscapeAnalyzer. This method is rarely
        used directly. Usually a LandscapeAnalyzer is initialized from a
        directory or file.
        """
        self._data = data
        self._software = software
        self._datapoints = datapoints

    @property
    def data(self):
        return self._data

    @property
    def software(self):
        return self._software

    @property
    def datapoints(self):
        return self._datapoints

    @classmethod
    def from_data(cls, data_dir, output_file='result.out', software='nwchem'):
        """
        Looks through all subdirectories of the provided directory and
        extracts all output.

        Args:
            directory:
            output_file:
            software:

        Returns:

        """

        # TODO Make subdirectory finder recursive?
        # Currently the method only checks in the immediate subdirectories
        # of the directory provided by the user. It might be useful to
        # extend this to subdirectories further down the tree. However,
        # maybe it's best to stick to calculating landscape properties in
        # the same directory.

        # Check if the directory exists
        if not os.path.isdir(data_dir):
            raise IOError("Directory " + data_dir + " not found.")

        # Find all the subdirectories in the specified directory
        dir_list = [d for d in os.listdir(data_dir)
                    if os.path.isdir(os.path.join(data_dir, d))]

        # Get all of the output of the calculations in the directory

        output = []

        if software == 'nwchem':
            for directory in dir_list:
                try:
                    out = nw.NwOutput(
                        os.path.join(data_dir, directory, output_file)
                    )

                except FileNotFoundError:
                    print('Did not find output file in directory ' + directory)
                except IndexError:
                    print('Did not find output in ' + output_file +
                          ' in directory ' + directory)
                else:
                    # Check if the output has an error
                    if out.data[-1]['has_error']:
                        print(
                            'Error found in output in directory ' + directory)
                    else:
                        if out.data[-1]['task_time'] == 0:
                            print('No timing data found for final task. '
                                  'Calculation might not have completed '
                                  'successfully.')

                        output.append(out)
                        print('Grabbed output in directory ' + directory)

            # TODO This currently only considers the final NwChem task.
            # Not a very general approach, but maybe the most sensible?
            data = [out.data[-1] for out in output]

        else:
            raise NotImplementedError("Only NwChem is currently supported.")

        return LandscapeAnalyzer(data, software=software)

    def analyze_cation_energies(self, coordinates, reference=None,
                                cation="Li"):
        """
        Extract the total energies for all the calculations, and connect
        them with the proper coordinates. Currently supports the following
        landscapes for analysis:

        Interfacet wedges: 2D landscapes that connect two Facets via a
        wedge. Polar coordinates are used to plot this landscape, and the
        program needs a reference axis to determine an appropriate angle.

        Spherical landscapes:

        :return:
        """

        datapoints = []
        dtype = None

        if coordinates == "polar":

            facet = reference

            dtype = [('Distance', float), ('Angle', float), ('Energy', float)]

            # If a facet is not provided, try to find the closest one to the
            # first ion coordinate data point. This might not always work.
            if reference is None:

                print(
                    "No Facet was provided. Since the facet is important for "
                    "defining the coordinates of the landscape, the program "
                    "will automatically determine the closest facet to the "
                    "cation in the first datapoint."
                )

                cage_init = Cage.from_molecule(self.data[0]['molecules'][0])
                init_cation_coords = [
                    site.coords for site in cage_init.sites
                    if site.specie == pmg.Element(cation)][-1]
                facet_init = cage_init.facets[0]
                for cage_facet in cage_init.facets:
                    if utils.distance(cage_facet.center, init_cation_coords) \
                            < utils.distance(facet_init.center,
                                             init_cation_coords):
                        facet = cage_facet

            for data in self.data:

                # Extract the cartesian coordinates of the cation
                cage = Cage.from_molecule(data["molecules"][0])
                cation_coords = [site.coords for site in cage.sites
                                 if site.specie == pmg.Element(cation)]

                # Check to see how many cations are found in the structure
                if len(cation_coords) == 0:
                    # If no cation is found, raise an error
                    raise ValueError("Requested cation not found in molecule.")
                elif len(cation_coords) == 1:
                    cation_coords = cation_coords[0]
                else:
                    # Take the last cation, this one is usually the one that
                    # is part of the landscape
                    cation_coords = cation_coords[-1]

                # Determine the polar coordinates with the facet center as a
                # reference axis for the angle
                r = np.linalg.norm(cation_coords - cage.anion_center)
                theta = angle_between(facet.center, cation_coords)

                if theta > math.pi / 2:
                    print("GOTCHA")
                    theta = math.pi - theta

                coordinate = [r, theta]

                energy_final = data['energies'][-1]
                coordinate.append(energy_final)

                datapoints.append(coordinate)

        if coordinates == "spherical":

            axis = reference

            dtype = [('Theta', float), ('Phi', float), ('Energy', float)]

            if reference is None:

                raise IOError("No reference axis was provided for the "
                              "analysis of the landscape energies. The "
                              "program currently has no recourse for "
                              "determining this axis automatically.")

            # Find a suitable perpendicular axis
            cage_init = Cage.from_molecule(self.data[0]['molecules'][0])

            i = 0
            phi_axis = None
            while i < len(cage_init) and phi_axis is None:
                v = cage_init.sites[i].coords
                if axis.dot(cage_init.sites[i].coords) > 1e-2:
                    phi_axis = unit_vector(
                        perpendicular_part(v, axis)
                    )
                i += 1

            for data in self.data:

                # Extract the cartesian coordinates of the cation
                cage = Cage.from_molecule(data["molecules"][0])
                cation_coords = [site.coords for site in cage.sites
                                 if site.specie == pmg.Element(cation)]

                # Check to see how many cations are found in the structure
                if len(cation_coords) == 0:
                    # If no cation is found, raise an error
                    raise ValueError("Requested cation not found in molecule.")
                elif len(cation_coords) == 1:
                    cation_coords = cation_coords[0]
                else:
                    # Take the last cation, this one is usually the one that
                    # is part of the landscape
                    cation_coords = cation_coords[-1]

                # Determine the spherical coordinates with the reference axis
                # for the angles
                theta = angle_between(axis, cation_coords)
                phi = angle_between(
                    perpendicular_part(cation_coords, axis), phi_axis
                )
                if angle_between(np.cross(phi_axis, axis), cation_coords)\
                        > \
                        math.pi/2:
                    phi = 2*math.pi - phi

                coordinate = [theta, phi]

                energy_final = data['energies'][-1]
                coordinate.append(energy_final)

                datapoints.append(coordinate)
                print("Appended datapoint:\n"
                      "E = " + str(energy_final) + "\n"
                      "coordinate = " + str(theta) + "," + str(phi))

        data_tuples = [tuple(point) for point in datapoints]

        darray = np.array(data_tuples, dtype=dtype)

        self._datapoints = darray

    def flip_coordinates(self, coord_name):
        """
        Flip the coordinate values for a chosen coordinate.
        :param coord_name:
        :return:
        """
        if len(self._datapoints) == 0:
            raise ValueError('No processed data present.')

        data = self._datapoints
        data[coord_name] = data[coord_name].max() - data[coord_name]
        self._datapoints = data

    def plot_energies(self, dimension, coordinates="polar", style='trisurf'):
        """
        Plot the energy landscape from the datapoints.

        Args:
            dimension:
            coordinates:
            style:

        Returns:

        """
        if len(self._datapoints) == 0:
            self.analyze_cation_energies(coordinates=coordinates)

        data = self.datapoints

        data['Distance'] = np.round(data['Distance'], 5)
        data = np.sort(data, order=['Distance', 'Angle'])

        if dimension == 1:
            plt.figure()
            plt.xlabel('Distance (Angstrom)')
            plt.ylabel('Energy (eV)')
            plt.plot(data['Distance'], data['Energy'])
            plt.show()

        if dimension == 2:
            if style == 'trisurf':
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.plot_trisurf(data['Distance'], data['Angle'],
                                data['Energy'], linewidth=0.2,
                                antialiased=True)
                plt.show()

            elif style == 'pcolor':
                # Find the number of radii and angles
                r_init = data['Distance'][0]
                nangles = 1
                while abs(data['Distance'][nangles] - r_init) < 1e-5:
                    nangles += 1
                nradii = int(len(data) / nangles)

                # Get the right format for the data
                radii = data['Distance'].reshape(nradii, nangles)  # [::nradii]
                angles = data['Angle'].reshape(nradii, nangles)  # [:nangles]
                energy = data['Energy'].reshape(nradii, nangles)

                # Plot
                fig, ax = subplots()
                p = ax.pcolor(angles, radii, energy, vmin=energy.min(),
                              vmax=energy.mean())
                cb = fig.colorbar(p)
                plt.xlabel('Angle (radians)')
                plt.ylabel('Distance (Angstrom)')
                plt.show()

    def as_dict(self):
        """
        Return a dictionary representing the LandscapeAnalyzer instance.
        :return:
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__}

        data_list = []
        for chunk in self.data:
            data_dict = chunk.copy()

            if data_dict['structures']:
                struc_dicts = []
                for structure in data_dict['structures']:
                    struc_dicts.append(structure.as_dict())
                data_dict['structures'] = struc_dicts

            if data_dict['molecules']:
                mol_dicts = []
                for molecule in data_dict['molecules']:
                    mol_dicts.append(molecule.as_dict())
                data_dict['molecules'] = mol_dicts

            data_list.append(data_dict)

        d["data"] = data_list
        d["datapoints"] = self.datapoints

        # TODO MAKE THIS WORK datapoints after analysis
        return d

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
    def from_file(cls, filename, fmt='json'):
        """
        Initialize an instance of LandscapeAnalyzer from a file.
        :return:
        """
        data = []

        if fmt == "json":

            with open(filename, "r") as f:
                file_data = json.load(f)

            for chunk in file_data["data"]:
                data_dict = chunk.copy()

                if data_dict['structures']:
                    structures = []
                    for struc_dict in data_dict['structures']:
                        structures.append(pmg.Structure.from_dict(struc_dict))
                    data_dict['structures'] = structures

                if data_dict['molecules']:
                    molecules = []
                    for mol_dict in data_dict['molecules']:
                        molecules.append(pmg.Molecule.from_dict(mol_dict))
                    data_dict['molecules'] = molecules

                data.append(data_dict)

        else:
            raise NotImplementedError("Only json format is currently "
                                      "supported.")

        return LandscapeAnalyzer(data, datapoints=file_data["datapoints"])

    def to(self, filename, fmt="json"):
        """
        Write the landscape to a file for quick reading.
        :return:
        """
        if fmt == "json":
            with open(filename, "w") as f:
                json.dump(self.as_dict(), f)
            return
        else:
            raise NotImplementedError("Only json format is currently "
                                      "supported.")
        # TODO Add support for .yaml

def perpendicular_part(v1, v2):
    """ Returns the projection of v1 on the plane perpendicular to v2. """
    return v1 - v1.dot(v2) * v2 / np.linalg.norm(v2)**2


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
    Return the rotation matrix associated with clockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

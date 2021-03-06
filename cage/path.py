# Encoding: utf-8
# Written for python 3.6

import pymatgen as pmg
import numpy as np
import re

import matplotlib.pyplot as plt
import scipy.interpolate as inter

from pymatgen.core.structure import Molecule
from monty.io import zopen
from pymatgen.core.units import Energy

"""
Tools for analyzing pathways on Cage molecules.

"""


class Path(object):
    """
    A Path is defined by a List of SiteCollections.

    """

    # Allowed distance between SiteCollection sites for checking if two paths can
    # be added.
    DISTANCE_TOL = 0.05

    def __init__(self, site_collections):
        """
        Initializes an instance of a Path from a List of SiteCollections.

        Args:
            site_collections (list):
        """
        self._site_collections = site_collections
        self._energies = None

    def __add__(self, other):
        """
        Add one Path to another. This will simply make a new Path which
        contains the SiteCollections of both paths.

        Will only work of the final SiteCollection of the first Path matches
        the initial SiteCollection of the second Path.

        The energies of the SiteCollections are only included if both Paths
        have them defined.

        Args:
            other (Path): Second path which to connect to this one.

        Returns:
            Path: Sum of the two paths.
        """
        d = [np.linalg.norm(site1.coords - site2.coords) for site1, site2 in
             zip(self.site_collections[-1].sites, other.site_collections[0].sites)]

        if max(d) < Path.DISTANCE_TOL:
            path = self.__class__(
                self.site_collections[:-1] + other.site_collections)
            if self.energies is None or other.energies is None:
                return path
            else:
                path.set_energies(self.energies[:-1] + other.energies)
                return path
        else:
            raise ValueError('Path endpoints do not match.')

    @property
    def site_collections(self):
        """
        The List of SiteCollections that comprise the Path.

        """
        return self._site_collections

    @property
    def energies(self):
        """
        The energies that correspond to each SiteCollection.

        """
        return self._energies

    @property
    def barrier(self):
        """
        Size of the energy barrier along the path. This is defined as the
        largest energy minus the initial energy. Expressed in eV.

        Returns:
            float: Size of the energy barrier in eV.

        """
        barrier_energy = Energy(max(self.energies) - self.energies[0], "Ha")
        return barrier_energy.to("eV")

    def set_energies(self, energies):
        """

        :param energies:
        :return:
        """
        if len(energies) == len(self.site_collections):
            self._energies = energies
        else:
            raise ValueError("Number of provided energies does not correspond "
                             "to the number of SiteCollections.")

    def flip(self):
        """
        Flip the direction of the path.

        Returns:
            None
        """
        self._site_collections.reverse()
        self._energies.reverse()

    def set_up_neb(self):
        """
        Sets up the NEB path for a NwChem calculation.

        """
        pass  # TODO

    @classmethod
    def from_file(cls, filename, fmt="xyz"):
        """
        Read the path from the NwChem input.neb.xyz file.

        Args:
            filename (str): Path to the file.
            fmt (str): Format of the file. Currently only support .xyz.

        Returns:
            Path
        """
        # TODO Change units to eV immediately
        # TODO determine format from the filename, as in Structure.from_file()
        if fmt == "xyz":
            with zopen(filename) as f:
                data = f.read()

            lines = data.split("\n")
            lines = [line for line in lines if line is not '']

            number_sites = int(lines[0].strip())
            number_molecules = int(len(lines) / number_sites)

            # Find the molecules
            mol_lines = [lines[i * (number_sites + 2):(i + 1) * (number_sites + 2)]
                         for i in range(number_molecules)]
            mol_strings = [''.join([line + '\n' for line in mol])
                           for mol in mol_lines]

            molecules = [Molecule.from_str(mol_string, fmt='xyz')
                         for mol_string in mol_strings]

            # Find the energies
            energies = []
            for comment in lines[1::number_sites + 2]:
                if re.search("energy", comment):
                    energies.append(float(comment.split()[-1]))

            path = cls(molecules)
            path.set_energies(energies)

            return path
        else:
            raise NotImplementedError("Requested format is not supported")

    def plot_energies(self, energy_range=None, interpolation='cubic spline'):
        """

        :return:
        """

        # Convert the energies from Hartree to eV
        energies = np.array([Energy(energy, "Ha").to("eV")
                             for energy in self.energies])

        energies = (energies - energies.min()) * 1000

        if energy_range is None:
            pass

        images = np.linspace(0, len(self.energies) - 1, num=len(self.energies))

        plt.figure()

        if interpolation == 'cubic spline':
            tck = inter.splrep(images, energies, s=0.01)
            images_inter = np.mgrid[images.min():images.max():0.01]
            energies_inter = inter.splev(images_inter, tck, der=0)
            plt.plot(images_inter, energies_inter, color='#143264')

        elif interpolation == 'cubic':
            energies_inter = inter.interp1d(images, energies, kind='cubic')
            images_inter = np.linspace(0, images.max(), num=100, endpoint=True)
            plt.plot(images_inter, energies_inter(images_inter), color='#143264')

        plt.plot(images, energies, 'o', color='#143264')
        plt.xticks([], [])
        plt.xlabel('Diffusion Pathway', size='x-large')
        plt.ylabel('Energy (meV)', size='x-large')
        plt.show()

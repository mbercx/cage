# Encoding: utf-8
# Written for python 3.6

import pymatgen as pmg
import numpy as np
import re

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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
    def __init__(self, site_collections):
        """
        Initializes an instance of a Path from a List of SiteCollections.
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
        :param other:
        :return:
        """
        # TODO Give an error if paths dont match
        path = self.__class__(self.site_collections + other.site_collections)
        if self.energies == None or other.energies == None:
            return path
        else:
            path.set_energies(self.energies + other.energies)
            return path

    @property
    def site_collections(self):
        """

        :return:
        """
        return self._site_collections

    @property
    def energies(self):
        """

        :return:
        """
        return self._energies

    def set_energies(self, energies):
        """

        :param energies:
        :return:
        """
        if len(energies) == len(self.site_collections):
            self._energies = energies
        else:
            raise ValueError("Number of provided energies does not correspond"
                             "to the number of SiteCollections.")

    def set_up_neb(self):
        """
        Sets up the NEB path for a NwChem calculation.
        :return:
        """
        pass

    @classmethod
    def from_file(cls, filename, fmt="xyz"):
        """
        Read the path from a file.
        :param filename:
        :return:
        """
        if fmt == "xyz":
            with zopen(filename) as f:
                data = f.read()

            lines = data.split("\n")
            lines = [line for line in lines if line is not '']

            number_sites = int(lines[0].strip())
            number_molecules = int(len(lines)/number_sites)

            # Find the molecules
            mol_lines = [lines[i*(number_sites + 2):(i+1)*(number_sites + 2)]
                         for i in range(number_molecules)]
            mol_strings = [''.join([line+'\n' for line in mol])
                           for mol in mol_lines]

            molecules = [Molecule.from_str(mol_string, fmt='xyz')
                         for mol_string in mol_strings]

            # Find the energies
            energies = []
            for comment in lines[1::number_sites+2]:
                if re.search("energy", comment):
                    energies.append(float(comment.split()[-1]))

            path = cls(molecules)
            path.set_energies(energies)

            return path

    def plot_energies(self):
        """

        :return:
        """

        energies = np.array([Energy(energy, "Ha").to("eV")
                             for energy in self.energies])
        energies = (energies - energies.min())*1000
        images = np.linspace(0, len(self.energies)-1, num=len(self.energies))

        e_inter = interp1d(images, energies, kind='cubic')

        im_inter = np.linspace(0, images.max(), num=100, endpoint=True)

        plt.figure()
        plt.plot(im_inter, e_inter(im_inter), 'black')
        plt.plot(images, energies, 'bo')
        plt.xlabel('Image Number')
        plt.ylabel('Energy (meV)')
        plt.show()










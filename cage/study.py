# coding: utf-8
# Copyright (c) Uncle Sam

import pymatgen as pmg
import pymatgen.io.nwchem as nwchem

from monty.json import MSONable

import os

"""
Defines a Study class, that contains all the information for a study of a
selection of structures. It allows the user to easily set up a lot of
calculations in a systematic way.
"""

__author__ = "Marnik Bercx"
__version__ = "0.1"
__maintainer__ = ""
__email__ = "marnik.bercx@uantwerpen.be"
__status__ = "Nowhere"
__date__ = "16 JUN 2017"

class Study(MSONable):
    """
    A Study is a List of Tasks which need to be performed on a List of
    Structures.
    """

    def __init__(self, structures, tasks, software='nwchem'):
        """
        Initialise the Study.

        :param systems: List of IStructures to study
        :param tasks: List of NwTasks or whatever we will use for VASP
        :param options: Dictionary of options to apply to calculations
        :param software: String that describes the code to use for the calculations
        """
        self._structures = structures
        self._tasks = tasks
        self._software = software
        self._comp_dict = classify_by_composition(structures)

    @property
    def structures(self):
        return self._structures

    @property
    def tasks(self):
        return self._tasks

    @property
    def software(self):
        return self._software

    def set_up_input(self, directory, sort_comp=False, **kwargs):
        """
        Set up all of the input files for the calculations and systems.

        :param directory: The full path to the directory where the calculation
        should be set up.
        :return:
        """
        #TODO Let the script set up the directory in case it does not exist

        # Get the absolute path to the directory
        abs_dir = os.path.abspath(directory)

        ##############
        #   NWCHEM   #
        ##############

        if self.software == 'nwchem':

            # Set up the directory tree and input files
            for comp in self._comp_dict.keys():

                if sort_comp:
                    comp_dir = str(comp).replace(' ','')

                    try:
                        os.mkdir(os.path.join(abs_dir,comp_dir))
                    except FileExistsError:
                        print(comp_dir + ' exists, skipping...')
                else:
                    comp_dir = ''

                geo_number = 1

                for structure in self._comp_dict[comp]:

                    geo_dir = 'geo' + str(geo_number)

                    try:
                        os.mkdir(os.path.join(abs_dir, comp_dir, geo_dir))
                    except FileExistsError:
                        print(geo_dir + ' exists, skipping...')

                    # Set up the input file
                    nwinput = nwchem.NwInput(structure, self.tasks,
                                             directives=self._options,
                                             **kwargs)
                    nwinput.write_file(os.path.join(abs_dir,
                                                    comp_dir, geo_dir,'input'))

                    geo_number += 1

        ############
        #   VASP   #
        ############

        #TODO Add functionality for VASP studies

        elif self.software == 'vasp':
            raise IOError("Currently, vasp inputs are not implemented yet.")

    def add_structure(self, structure):
        """
        Add a structure to the list of structures.
        """
        self._structures.append(structure)
        self._comp_dict = classify_by_composition(self.structures)

    def add_task(self):
        """
        Add a task to the list of tasks.
        :return:
        """
        pass

    def add_restraint(self):
        """
        Add a restraint to the list of restraints on the calculations.
        """
        pass

    def change_software(self):
        """
        Change the software being used for the calculations.
        :return:
        """
        pass


def classify_by_composition(structures):
    """
    Classify the different IStructures by composition.
    :param structures: List of structures
    :return: Dictionary of composition strings with a List of the corresponding
    IStructures.
    """
    comp_dic = {}
    for structure in structures:
        if comp_dic.get(str(structure.composition),False):
            comp_dic[str(structure.composition)].append(structure)
        else:
            comp_dic[str(structure.composition)] = [structure,]
    return comp_dic
# Encoding = utf-8

import numpy as np
import pymatgen.io.nwchem as nw

def write_init_final(filename):
    """
    Quickly write out the initial and the final configuration of a nwchem
    optimization.
    :param filename:
    :return:
    """
    data = nw.NwOutput(filename).data[-1]
    data['molecules'][0].to(fmt='xyz',filename='initial_mol.xyz')
    data['molecules'][-1].to(fmt='xyz', filename='final_mol.xyz')

def distance(coord1, coord2):
    """
    Calculate the distance between two coordinates, defined by arrays.
    :param coord1:
    :param coord2:
    :return:
    """
    return np.linalg.norm(coord1 - coord2)
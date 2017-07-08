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

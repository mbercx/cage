# Encoding = utf-8

import numpy as np

def distance(coord1, coord2):
    """
    Calculate the distance between two coordinates, defined by arrays.
    :param v1:
    :param v2:
    :return:
    """
    return np.linalg.norm(coord1 - coord2)
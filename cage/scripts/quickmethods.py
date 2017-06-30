# Encoding = utf-8

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
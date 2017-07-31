
import pymatgen.io.nwchem as nwchem

"""
Small utility scripts that support the use of the cage package.

"""

def geo(output_file):
    """
    Quickly write out the initial and the final configuration of a nwchem
    optimization to xyz files.

    Args:
        output_file (str): Output file of the nwchem calculation.
    """
    data = nwchem.NwOutput(output_file).data[-1]
    data['molecules'][0].to(fmt='xyz',filename='initial_mol.xyz')
    data['molecules'][-1].to(fmt='xyz', filename='final_mol.xyz')
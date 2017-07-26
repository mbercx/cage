"""
cage

Usage:
    cage -h | --help
    cage dock
    cage chain

Options:
    -h --help

"""

from cage import __version__ as VERSION
import click

@click.group()
def main():
    """
    Python tools to set up NwChem calculations for investigating the local
    cation-anion interaction around cage-like molecules.
    """
    pass

main.command

@main.command()
@click.argument('filename')
@click.option('--cation', '-C', default='Li')
@click.option('--distance', '-d', default=2)
def dock(filename, cation, distance):
    """ Set up the docking sites """
    from cage.cli.commands.setup import docksetup

    docksetup(filename, cation, distance)

@main.command()
@click.argument('filename')
@click.option('--cation', '-C', default='Li')
@click.option('--operation', '-O', default='energy')
@click.option('--endradii', '-e', default=(3, 6))
@click.option('--nradii', default=30)
@click.option('--adensity', default=50)
def chain(filename, cation, operation, endradii, nradii, adensity):
    """ Set up the non-equivalent paths between facets. """
    from cage.cli.commands.setup import chainsetup

    chainsetup(filename, cation, operation, endradii, nradii, adensity)

@main.group()
def twocat():
    """
    Set up calculations for anions with two cations. These require the results
    from docking calculations for single cations on the anion.
    """
    pass

@twocat.command()
@click.argument('dock_dir')
@click.option('--cation', '-C', default='Li')
@click.option('--operation', '-O', default='energy')
@click.option('--endradii', '-e', default=(3, 6))
@click.option('--nradii', default=30)
@click.option('--adensity', default=50)
def chain(dock_dir, cation, operation, endradii, nradii, adensity):
    """
    Similar to the single cation case, this command sets up a 2D landscape
    between the normals of a chain of non-equivalent facets.
    """
    from cage.cli.commands.setup import twocat_chainsetup

    twocat_chainsetup(dock_dir, cation, operation, endradii, nradii, adensity)
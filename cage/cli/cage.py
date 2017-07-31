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
    Python tools for studying cage-like anions for solid electrolytes in
    batteries.
    """
    pass

@main.group()
def setup():
    """
    Python tools to set up NwChem calculations for investigating the local
    cation-anion interaction around cage-like molecules.
    """
    pass

@setup.command()
@click.argument('filename')
def optimize(filename):
    """ Set up the initial anion optimization. """
    from cage.cli.commands.setup import optimize

    optimize(filename)

@setup.command()
@click.argument('filename')
@click.option('--cation', '-C', default='Li')
@click.option('--distance', '-d', default=2)
def dock(filename, cation, distance):
    """
    Set up the docking sites. It is recommended to use an anion which has
    first been optimize using 'cage setup optimize'.
    """
    from cage.cli.commands.setup import docksetup

    docksetup(filename, cation, distance)

@setup.command()
@click.argument('filename')
@click.option('--cation', '-C', default='Li')
@click.option('--operation', '-O', default='energy')
@click.option('--endradii', '-e', default=(3, 6))
@click.option('--nradii', default=30)
@click.option('--adensity', default=50)
def chain(filename, cation, operation, endradii, nradii, adensity):
    """ Set up a 2D landscape along the chain of non-equivalent facets. """
    from cage.cli.commands.setup import chainsetup

    chainsetup(filename, cation, operation, endradii, nradii, adensity)

@setup.command()
@click.argument('filename')
@click.option('--cation', '-C', default='Li')
@click.option('--distance', '-d', default=2)
@click.option('--edges', is_flag=True)
def path(filename, cation, distance, edges):
    """ Set up the paths between facets that share a vertex. """
    from cage.cli.commands.setup import pathsetup

    pathsetup(filename, cation, distance, edges)

@setup.command()
@click.argument('paths_dir')
@click.option('--nimages', '-n', default=10)
def neb(paths_dir, nimages):
    """ Set up the paths between facets that share a vertex. """
    from cage.cli.commands.setup import nebsetup

    nebsetup(paths_dir, nimages)

@setup.group()
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
@click.option('--tolerance', default=1e-2)
@click.option('--verbose', '-v', is_flag=True)
def chain(dock_dir, cation, operation, endradii, nradii, adensity, tolerance,
          verbose):
    """
    Similar to the single cation case, this command sets up a 2D landscape
    between the normals of a chain of non-equivalent facets.
    """
    from cage.cli.commands.setup import twocat_chainsetup

    twocat_chainsetup(dock_dir, cation, operation, endradii, nradii, adensity,
                      tolerance, verbose)

@main.group()
def analyze():
    """
    Python tools to analyze NwChem calculations for investigating the local
    cation-anion interaction around cage-like molecules.
    """
    pass

@analyze.command()
def landscape():
    """
    Analyze the landscape
    """

@main.group()
def util():
    """
    A set of utility scripts for the cage package.
    """
    pass

@util.command()
def geo():
    """ Write the initial and final geometry."""
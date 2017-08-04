
import click

"""
Command Line Interface for the cage package.

"""


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
    Set up Nwchem calculations for the anion.

    Tools for investigating the local cation-anion interaction around cage-like
    molecules.
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
@click.option('--distance', '-d', default=2.0)
@click.option('--verbose', '-v', is_flag=True)
def dock(filename, cation, distance, verbose):
    """
    Set up the docking sites.

    It is recommended to use an anion which has first been optimize using
    'cage setup optimize'.
    """
    from cage.cli.commands.setup import docksetup

    docksetup(filename, cation, distance, verbose)


@setup.command()
@click.argument('filename')
@click.option('--cation', '-C', default='Li')
@click.option('--facets', '-f', type=str ,  default='tuple')
@click.option('--operation', '-O', default='energy')
@click.option('--endradii', '-e', default=(3.0, 6.0))
@click.option('--nradii', default=30)
@click.option('--adensity', default=50)
def chain(filename, cation, facets, operation, endradii, nradii, adensity):
    """ Set up a 2D landscape along the chain of non-equivalent facets. """
    from cage.cli.commands.setup import chainsetup

    facets = [int(number) for number in facets.split()]

    chainsetup(filename=filename,
               facets=facets,
               cation=cation,
               operation=operation,
               endradii=endradii,
               nradii=nradii,
               adensity=adensity)


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
    """
    Set up the nudged elastic band calculation.

    """
    from cage.cli.commands.setup import nebsetup

    nebsetup(paths_dir, nimages)


@setup.group()
def twocat():
    """
    Set up calculations for anions with two cations.

    These require the results from docking calculations for single cations on
    the anion.
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
    Analyze the results from NwChem calculations.

    Allows a quick analysis of the output of several calculations to quickly
    visualize results.
    """
    pass


@analyze.command()
@click.argument('lands_dir')
@click.option('--cation', '-C', default='Li')
@click.option('--energy_range', '-E', default=(0.0, 0.0))
@click.option('--interp_mesh', '-I', default=(0.003, 0.01))
@click.option('--contour_levels', '-l', default=0.1)
@click.option('--verbose', '-v', is_flag=True)
def landscape(lands_dir, cation, energy_range, interp_mesh, contour_levels,
              verbose):
    """
    Analyze the landscape data.
    """
    from cage.cli.commands.analyze import landscape_analysis

    landscape_analysis(lands_dir=lands_dir,
                       cation=cation,
                       energy_range=energy_range,
                       interp_mesh=interp_mesh,
                       contour_levels=contour_levels,
                       verbose=verbose)


@main.group()
def util():
    """
    A set of utility scripts for the cage package.
    """
    pass


@util.command()
@click.argument('output_file')
def geo(output_file):
    """ Write the initial and final geometry of a nwchem optimization. """
    from cage.cli.commands.util import geo

    geo(output_file=output_file)


@util.command()
@click.argument('output_file')
def check(output_file):
    """  """
    from cage.cli.commands.util import check_calculation

    check_calculation(output_file=output_file)


@util.command()
@click.argument('output_file')
def process(output_file):
    """  """
    from cage.cli.commands.util import process_output

    process_output(output_file=output_file)


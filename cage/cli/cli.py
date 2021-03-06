import click

"""
Command Line Interface for the cage package.

"""

# This is only used to make '-h' a shorter way to access the CLI help
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    """
    Python tools for studying cage-like anions for solid electrolytes in
    batteries.
    """
    pass


# region * Setup
@main.group(context_settings=CONTEXT_SETTINGS,
            short_help="Set up Nwchem calculations for the anion.")
def setup():
    """
    Tools for setting up calculations for investigating the local cation-anion
    interaction around cage-like molecules.
    """
    pass


@setup.command(context_settings=CONTEXT_SETTINGS,
               short_help="Initial anion geometric optimization.")
@click.argument('filename')
@click.option('--charge', '-c', default=0)
def optimize(filename, charge):
    """ Set up the initial anion optimization. """
    from cage.cli.commands.setup import optimize

    optimize(filename, charge)


@setup.command(context_settings=CONTEXT_SETTINGS,
               short_help="Set up docking sites.")
@click.argument('filename')
@click.option('--cation', '-C', default='Li',
              help="The cation to be placed on the dock, provided as a string "
                   "of the chemical symbol, e.g. 'Li' or 'Na'.")
@click.option('--distance', '-d', default=2.0)
@click.option('--facets', '-f', type=str, default='tuple')
@click.option('--verbose', '-v', is_flag=True)
def dock(filename, cation, distance, facets, verbose):
    """
    Set up a geometry optimization of a cation docked on a facet in NwChem for
    all non-equivalent facets of the anion or a list of chosen facets.

    It is recommended to use an anion which has first been optimized using
    'cage setup optimize'.
    """
    from cage.cli.commands.setup import docksetup

    if facets == 'tuple':
        facets = tuple
    else:
        facets = [int(number) for number in facets.split()]

    docksetup(filename, cation, distance, facets, verbose)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.argument('filename')
@click.option('--cation', '-C', default='Li')
@click.option('--facets', '-f', type=str, default='tuple')
@click.option('--operation', '-O', default='energy')
@click.option('--end_radii', '-R', default=(3.0, 6.0))
@click.option('--nradii', default=30)
@click.option('--adensity', default=50)
def chain(filename, cation, facets, operation, end_radii, nradii, adensity):
    """ Set up a 2D landscape along a chain of facets. """
    from cage.cli.commands.setup import chainsetup

    if facets == 'tuple':
        facets = tuple
    else:
        facets = [int(number) for number in facets.split()]

    chainsetup(filename=filename,
               facets=facets,
               cation=cation,
               operation=operation,
               end_radii=end_radii,
               nradii=nradii,
               adensity=adensity)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.argument('filename')
@click.option('--cation', '-C', default='Li')
@click.option('--distance', '-d', default=2)
@click.option("--facets", "-f", type=str, default="tuple")
@click.option('--edges', is_flag=True)
def path(filename, cation, distance, facets, edges):
    """ Set up the paths between facets that share a vertex. """
    from cage.cli.commands.setup import pathsetup

    if facets == 'tuple':
        facets = tuple
    else:
        facets = [int(number) for number in facets.split()]

    pathsetup(filename, cation, distance, facets, edges)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.argument('paths_dir')
@click.option('--nimages', '-n', default=10)
def neb(paths_dir, nimages):
    """ Set up the nudged elastic band calculation. """
    from cage.cli.commands.setup import nebsetup

    nebsetup(paths_dir, nimages)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.argument("facet_index")
@click.argument("filename")
@click.option("--cation", "-C", default="Li",
              help="Cation for which to calculate the reference energy.")
@click.option("--distance", "-d", default=8.0,
              help="Maximum distance at which to place the cation for the "
                   "reference energy point.")
@click.option("--verbose", "-v", is_flag=True)
def reference(facet_index, filename, cation, distance, verbose):
    """Set up a calculation to determine a reference energy."""

    from cage.cli.commands.setup import reference

    # Convert string input to integer
    facet_index = int(facet_index)

    reference(facet_index=facet_index,
              filename=filename,
              cation=cation,
              end_radius=distance,
              verbose=verbose)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.argument('filename')
@click.option('--cation', '-C', default='Li')
@click.option('--radius', '-R', default=6.0)
@click.option('--density', default=20)
def sphere(filename, cation, radius, density):
    """ Set up a spherical landscape calculation. """
    from cage.cli.commands.setup import spheresetup

    spheresetup(filename=filename,
                cation=cation,
                radius=radius,
                density=density)


# endregion

# region * Setup - Twocat
@setup.group(context_settings=CONTEXT_SETTINGS)
def twocat():
    """
    Set up calculations for anions with two cations.

    These require the results from docking calculations for single cations on
    the anion.
    """
    pass


@twocat.command(context_settings=CONTEXT_SETTINGS)
@click.argument('dock_dir')
@click.option('--cation', '-C', default='Li')
@click.option('--operation', '-O', default='energy')
@click.option('--endradii', '-R', default=(3, 6))
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


# endregion

# region * Analyze
@main.group(context_settings=CONTEXT_SETTINGS,
            short_help="Analyze the results from NwChem calculations.")
def analyze():
    """
    Scripts to help analyze the output of several calculations to quickly
    visualize results.
    """
    pass


@analyze.command(context_settings=CONTEXT_SETTINGS)
@click.argument('lands_dir')
@click.option('--cation', '-C', default='Li')
@click.option('--energy_range', '-E', default=(0.0, 0.0))
@click.option('--interp_mesh', '-I', default=(0.03, 0.01))
@click.option('--end_radii', '-R', default=(0.0, 0.0))
@click.option('--contour_levels', '-l', default=0.1)
@click.option('--verbose', '-v', is_flag=True)
@click.option("--coulomb_charge", "-c", default=0)
@click.option("--reference_energy", "-r", default=0.0)
@click.option("--interp_method", "-m", default="griddata")
def landscape(lands_dir, cation, energy_range, interp_mesh, end_radii,
              contour_levels, verbose, coulomb_charge, reference_energy,
              interp_method):
    """ Analyze the landscape data. """
    from cage.cli.commands.analyze import landscape_analysis

    if reference_energy == 0.0:
        reference_energy = None

    landscape_analysis(directory=lands_dir,
                       cation=cation,
                       energy_range=energy_range,
                       interp_mesh=interp_mesh,
                       end_radii=end_radii,
                       contour_levels=contour_levels,
                       verbose=verbose,
                       coulomb_charge=coulomb_charge,
                       reference_energy=reference_energy,
                       interp_method=interp_method,
                       set_contour_levels_manually=True)


@analyze.command(context_settings=CONTEXT_SETTINGS)
@click.argument('lands_dir')
@click.option('--cation', '-C', default='Li')
@click.option('--interp_mesh', '-I', default=(0.03, 0.01))
@click.option('--end_radii', '-R', default=(0.0, 0.0))
@click.option('--verbose', '-v', is_flag=True)
@click.option("--coulomb_charge", "-c", default=0)
@click.option("--reference_energy", "-r", default=0.0)
def barrier(lands_dir, cation, interp_mesh, end_radii,
            verbose, coulomb_charge, reference_energy):
    """ Analyze the barriers in landscape data. """
    from cage.cli.commands.analyze import barrier_analysis

    if reference_energy == 0.0:
        reference_energy = None

    barrier_analysis(lands_dir=lands_dir,
                     cation=cation,
                     interp_mesh=interp_mesh,
                     end_radii=end_radii,
                     verbose=verbose,
                     coulomb=coulomb_charge,
                     reference_energy=reference_energy)


@analyze.command(context_settings=CONTEXT_SETTINGS)
@click.argument("reference_dir")
@click.option("--coulomb_charge", "-c", default=0)
def reference(reference_dir, coulomb_charge):
    from cage.cli.commands.analyze import reference

    reference(reference_dir, coulomb_charge)


@analyze.command(context_settings=CONTEXT_SETTINGS)
@click.argument('lands_dir')
@click.option('--cation', '-C', default='Li')
@click.option('--energy_range', '-E', default=(0.0, 0.0))
@click.option('--interp_mesh', '-I', default=(0.01, 0.01))
@click.option('--contour_levels', '-l', default=0.1)
@click.option("--reference_energy", "-r", default=0.0)
@click.option("--interp_method", "-m", default="griddata")
def sphere(lands_dir, cation, energy_range, interp_mesh, contour_levels,
           reference_energy, interp_method):
    """ Analyze the landscape data. """
    from cage.cli.commands.analyze import sphere_analysis

    if reference_energy == 0.0:
        reference_energy = None

    sphere_analysis(directory=lands_dir,
                    cation=cation,
                    interp_mesh=interp_mesh,
                    energy_range=energy_range,
                    contour_levels=contour_levels,
                    reference_energy=reference_energy,
                    interp_method=interp_method)


# endregion

# region * Util
@main.group(context_settings=CONTEXT_SETTINGS)
def util():
    """
    A set of utility scripts for the cage package.
    """
    pass


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument('output_file')
def geo(output_file):
    """ Write the initial and final geometry of a nwchem optimization. """

    from cage.cli.commands.util import geo

    geo(output_file=output_file)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument('output_file')
def energy(output_file):
    from cage.cli.commands.util import energy

    energy(output_file=output_file)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument('output_file')
def check(output_file):
    """ Check the output of calculations. """

    from cage.cli.commands.util import check_calculation

    check_calculation(output=output_file)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument('output_file')
def process(output_file):
    """ Process the output of calculations. """

    from cage.cli.commands.util import process_output

    process_output(output=output_file)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument('directory')
def gather(directory):
    """ Gather the results of a landscape calculation. """

    from cage.cli.commands.util import gather_landscape

    gather_landscape(directory=directory)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument('filename')
def visualize(filename):
    """ Visualize the facets of a molecule. """

    from cage.cli.commands.util import visualize_facets

    visualize_facets(filename=filename)


# endregion

# region * Workflow
@main.group(context_settings=CONTEXT_SETTINGS)
def workflow():
    """
    Workflow setup scripts.
    """
    pass


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument('filename')
@click.option('--cation', '-C', default='Li')
@click.option('--facets', '-f', type=str, default='tuple')
@click.option('--operation', '-O', default='energy')
@click.option('--end_radii', '-R', default=(3.0, 6.0))
@click.option('--nradii', default=30)
@click.option('--adensity', default=50)
def landscape(filename, cation, facets, operation, end_radii, nradii,
              adensity):
    """ Set up a 2D landscape along a chain of facets. """
    from cage.workflow import landscape_workflow

    if facets == 'tuple':
        facets = tuple
    else:
        facets = [int(number) for number in facets.split()]

    landscape_workflow(filename=filename,
                       facets=facets,
                       cation=cation,
                       operation=operation,
                       end_radii=end_radii,
                       nradii=nradii,
                       adensity=adensity)


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument('filename')
@click.option('--cation', '-C', default='Li')
@click.option("--radius", '-R', default=6.0)
@click.option('--density', default=20)
def sphere(filename, cation, radius, density):
    """ Set up a 2D spherical landscape. """
    from cage.workflow import sphere_workflow

    sphere_workflow(filename=filename,
                    cation=cation,
                    radius=radius,
                    density=density)


@workflow.command(context_settings=CONTEXT_SETTINGS,
                  short_help="Anion geometric optimization.")
@click.argument('filename')
@click.option('--charge', '-c', default=0)
def optimize(filename, charge):
    """ Set up the initial anion optimization. """
    from cage.workflow import optimize_workflow

    optimize_workflow(filename, charge)


@workflow.command(context_settings=CONTEXT_SETTINGS)
def test():
    """
    Testing if I can get workflows to work on quartz...

    Returns:

    """
    from cage.workflow import test_workflow

    test_workflow()

# endregion workflow
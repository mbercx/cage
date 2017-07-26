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
    """ Main CLI entry point."""
    pass

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
@click.option('--endpoints', '-e', default=(3, 6))
@click.option('--nradii', default=30)
@click.option('--adensity', default=50)
def chain(filename, cation, operation, endpoints, nradii, adensity):
    """ Set up the non-equivalent paths between facets. """
    from cage.cli.commands.setup import chainsetup

    chainsetup(filename, cation, operation, endpoints, nradii, adensity)
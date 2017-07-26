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
def dock(filename):
    """ Set up the docking sites """
    from cage.cli.commands.dock import docksetup

    docksetup(filename)


@main.command()
def paths():
    """ Set up the non-equivalent paths between facets. """
    print('Implement me, foolish human!')
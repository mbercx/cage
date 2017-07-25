"""
cage

Usage:
    cage -h | --help
    cage dock
    cage chain

Options:
    -h --help

"""

from docopt import docopt
from . import __version__ as VERSION

def main():
    """ Main CLI entry point."""

    import cage.commands
    options = docopt(__doc__, version=VERSION)
    print(options)
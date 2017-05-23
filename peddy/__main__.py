import sys

from .cli import peddy as cli

"""
peddy.__main__
~~~~~~~~~~~~~~~~~~~~~
The main entry point for the command line interface.
Invoke as ``python -m peddy``.
"""

if __name__ == "__main__":
    # exit using whatever exit code the CLI returned
    sys.exit(cli())

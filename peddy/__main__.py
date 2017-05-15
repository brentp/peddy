import sys

from peddy.cli import peddy as cli

"""
peddy.__main__
~~~~~~~~~~~~~~~~~~~~~
The main entry point for the command line interface.
Invoke as ``peddy`` (if installed)
or ``python -m peddy`` (no install required).
"""

if __name__ == "__main__":
    # exit using whatever exit code the CLI returned
    sys.exit(cli())
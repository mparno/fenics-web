# This script copies the Python docstrings module from DOLFIN to
# ../source/programmer-reference/python/docstrings such that building the docs
# only needs to look for the local copy.
# The script requires path/to/dolfin/site-packages/dolfin to be in PYTHON_PATH.

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2010-08-31"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-08-31

import os
import shutil

try:
    # Import the docstrings module.
    import docstrings
    # Get absolute path to docstrings module.
    src = os.path.abspath(docstrings.__path__[0])
    # Create destination name and delete iff it exists.
    dst = os.path.join(os.pardir, "source", "programmers-reference", "python", "docstrings")
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    # Dump module in destination directory.
    print "\n  Creating local copy of the Python docstrings module from DOLFIN"
    print "  (destination directory: %s" % dst
    shutil.copytree(src, dst)
    print "  Done!\n"
except Exception as what:
    print "  Could not import the Python docstrings module from DOLFIN"
    print "  (error: %s),\n update your PYTHONPATH variable?" % what
    raise ImportError


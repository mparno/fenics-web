"""
Integrate generated documentation in .rst files with fenics-web
"""

__authors__ = "Marie E. Rognes <meg@simula.no>"
__copyright__ = "Copyright (C) 2011 " + __authors__
__license__  = "GNU LGPL version 3 or any later version"

# Last changed: 2011-05-13

import os, shutil

def init_sphinx_project_for_dolfin(project_dir, web_dir):

    # Do sphinx-quickstart by hand.
    try:
        os.makedirs(os.path.join(project_dir, "source"))
    except:
        pass
    try:
        os.makedirs(os.path.join(project_dir, "build"))
    except:
        pass

    shutil.copy(os.path.join(web_dir, "Makefile"), project_dir)
    shutil.copy(os.path.join(web_dir, "source", "conf.py"),
                os.path.join(project_dir, "source"))
    shutil.copy(os.path.join(web_dir, "source", "index.rst"),
                os.path.join(project_dir, "source"))
    shutil.copy(os.path.join(web_dir, "source", "mathjax.py"),
                os.path.join(project_dir, "source"))
    try:
        shutil.copytree(os.path.join(web_dir, "source", "_themes"),
                        os.path.join(project_dir, "source", "_themes"))
    except:
        print "Couldn't update themes. Ignoring themes updates"

if __name__ == "__main__":

    version = "0.N.M"

    # The path to the main fenics web directory
    web_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               os.path.pardir, os.path.pardir)

    # The path to directory where sphinx-project should be made. (Tune
    # at will)
    project_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    os.path.pardir, os.path.pardir,
                                    "documentation", "dolfin",
                                    "dolfin-%s" % version)
    cpp_project_dir = os.path.join(project_dir, "cpp")
    init_sphinx_project_for_dolfin(cpp_project_dir, web_dir)

    python_project_dir = os.path.join(project_dir, "python")
    init_sphinx_project_for_dolfin(python_project_dir, web_dir)

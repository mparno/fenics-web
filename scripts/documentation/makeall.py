from generaterest import *
from integraterest import *

import os

def generate_doc(version, project_dir):

    # Get input directory: FIXME
    input_dir = "/home/meg/local/fenics/src/dolfin"

    # Generate complete dolfin documentation
    print "Generating documentation in %s" % project_dir
    generate_dolfin_doc(version, input_dir, project_dir)

    return project_dir

def initialize_sphinx(project_dir):

    # The path to the main fenics web directory
    web_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               os.path.pardir, os.path.pardir)

    # Init sphinx repositories for cpp and python
    cpp_project_dir = os.path.join(project_dir, "cpp")
    init_sphinx_project_for_dolfin(cpp_project_dir, web_dir)

    python_project_dir = os.path.join(project_dir, "python")
    init_sphinx_project_for_dolfin(python_project_dir, web_dir)

    return (cpp_project_dir, python_project_dir)

def generate_html(project_dirs):

    # Go to these repositories and do make html
    for project_dir in project_dirs:
        os.chdir(project_dir)
        os.system("make html")

    return

def main(version):

    # Set output directory: everything else will be generated relative
    # to this directory. (Tune at will)
    project_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               os.path.pardir, os.path.pardir,
                               "documentation", "dolfin",
                               "dolfin-%s" % version)

    # Generate .rst documentation
    generate_doc(version, project_dir)

    # Initialize sphinc projects
    project_dirs = initialize_sphinx(project_dir)

    # Call sphinx (make html)
    generate_html(project_dirs)


if __name__ == "__main__":

    version = "0.N.M"
    main(version)

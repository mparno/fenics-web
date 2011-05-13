__authors__ = "Anders Logg <logg@simula.no>, Kristian Oelgaard <k.b.oelgaard@gmail.com> and Marie E. Rognes <meg@simula.no>"
__copyright__ = "Copyright (C) 2010 " + __authors__
__license__  = "GNU LGPL version 3 or any later version"

# Last changed: 2011-05-13

__all__ = ["generate_dolfin_doc"]

import os

from generatecpprest import generate_cpp_api_documentation
from generatepythonrest import generate_python_api_documentation

def generate_demo_documentation(language, version):

    """
    Generate/Copy demo documentation

    Arguments:

      * language
        language should be either "cpp" or "python"

      * version
        version should be the version number, example: 0.9.11 or dev
    """

    print "Copying %s demo documentation for DOLFIN-%s" % (language, version)

    return

def generate_dolfin_doc(version, input_dir, output_dir):

    # Make output directory (or use current if existing)
    try:
        os.makedirs(output_dir)
    except:
        pass

    # Generate .rst for C++ documentation
    api_output_dir = os.path.join(output_dir, "cpp", "source",
                                  "programmers-reference")
    generate_cpp_api_documentation(version, input_dir, api_output_dir)
    generate_demo_documentation("cpp", version)

    # Generate .rst for Python documentation
    api_output_dir = os.path.join(output_dir, "python", "source",
                                  "programmers-reference")
    generate_python_api_documentation(version, api_output_dir)
    generate_demo_documentation("python", version)

if __name__ == "__main__":

    version = "0.N.M"

    # Set output directory: everything else will be generated relative
    # to this directory. (Tune at will)
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              os.path.pardir, os.path.pardir,
                              "documentation", "dolfin",
                              "dolfin-%s" % version)

    # Get input directory FIXME
    input_dir = "/home/meg/local/fenics/src/dolfin"

    # Generate complete dolfin documentation
    print "Generating documentation in %s" % output_dir
    generate_dolfin_doc(version, input_dir, output_dir)



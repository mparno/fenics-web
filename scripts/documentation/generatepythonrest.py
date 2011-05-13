"""
Generate .rst files for Python DOLFIN Programmer's Reference in a very
basic manner using Sphinx automodule for each .py file
"""

__authors__ = "Marie E. Rognes <meg@simula.no>"
__copyright__ = "Copyright (C) 2010 " + __authors__
__license__  = "GNU LGPL version 3 or any later version"

# Last changed: 2011-05-13

import os

# Import some utilities
from dolfin_utils.documentation import indent
import indextemplates as templates

def generate_python_api_documentation(version, output_dir):
    """
    Generate .rst files for Python DOLFIN Programmer's Reference

    Arguments:
      * version
        version should be the version number, example 0.9.11 or dev
    """

    print "Generating .rst for Python API documentation of DOLFIN-%s" % version

    try:
        import dolfin
    except:
        print "Couldn't import dolfin"

    input_dir = os.path.dirname(dolfin.__file__)
    print "input_dir = ", input_dir
    print "output_dir = ", output_dir
    try:
        os.makedirs(output_dir)
    except:
        print "Dirs already existing %s" % output_dir

    for root, dirs, files in os.walk(input_dir):
        for directory in dirs:
            directory_name = os.path.join(output_dir, directory)
            try:
                os.makedirs(directory_name)
            except:
                print "Dirs already existing %s" % directory_name

        for filename in files:
            (basename, ext) = os.path.splitext(filename)
            if ext == ".py" and basename != "__init__":

                abs_filename = os.path.join(root, filename)
                rel_filename =  os.path.relpath(abs_filename, input_dir)

                # Make rst file for this file
                py_filename = os.path.join(output_dir, rel_filename)
                rst_basename = os.path.splitext(py_filename)[0]
                print "rst_basename = ", rst_basename
                rst_filename = rst_basename + ".rst"
                #print rst_filename
                print

                file = open(rst_filename, "w")

                foo = os.path.relpath(rst_basename, output_dir)
                foo = foo.replace("/", ".")
                module_name = "dolfin." + foo
                print "module_name = ", module_name
                text = templates.automodule_template % module_name
                file.write(text)
                file.close()
        print

    return



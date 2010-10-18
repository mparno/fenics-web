# This script generates the reST input for the C++ Programmer's
# Reference from the DOLFIN C++ header files found in DOLFIN_DIR.
#
# The script looks for comments starting with /// and appearing
# before class and function declarations. The comments are copied
# verbatim to the generate reST files.
#
# Ths script also replaces _Foo_ with links to class pages.

__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2010-08-26"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2010.

# Last changed: 2010-10-06

import os, sys

# Set output directory
output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),\
                          os.pardir, "source", "programmers-reference", "cpp")

# Set directory for DOLFIN source code
if not "DOLFIN_DIR" in os.environ:
    raise RuntimeError, "You need to set the DOLFIN_DIR environment variable."
# Make sure we have an absolute path.
dolfin_dir = os.path.abspath(os.environ["DOLFIN_DIR"])

# Add path to dolfin_utils and import the documentation extractor.
#sys.path.append(os.path.abspath(os.path.join(dolfin_dir, "site-packages")))
sys.path = [os.path.abspath(os.path.join(dolfin_dir, "site-packages"))] + sys.path
from dolfin_utils.documentation import extract_doc_representation, indent, add_links

def write_documentation(documentation, header, module, classnames):
    "Write documentation for given header in given module"

    # For quick testing
#    if not header == "Mesh.h":
#        return

    print "Writing documentation for %s..." % header

    # Create containing directory
    directory = os.path.join(output_dir, module)
    try:
        os.makedirs(directory)
    except:
        pass

    # Set location of documentation
    prefix = header.split(".")[0]
    rstfile = prefix + ".rst"
    outfile = os.path.join(directory, rstfile)

    # Check that it is ok to overwrite the output file
    try:
        f = open(outfile, "r")
        if "DONT_TOUCH" in f.read():
            print "Not touching file %s" % outfile
            f.close()
            return
        f.close()
    except:
        pass

    # Open output file
    output = ""
    output += ".. Documentation for the header file dolfin/%s/%s\n" % (module, header)
    output += "\n"
    output += ".. _programmers_reference_cpp_%s_%s:\n" % (module, prefix.lower())
    output += "\n"
    output += header + "\n"
    output += len(header)*"=" + "\n"
    output += "\n"

    # Write a note that this file was generated automatically
    output += ".. note::\n"
    output += "\n"
    output += indent("""\
The documentation on this page was automatically extracted from
the DOLFIN C++ code and may need to be edited or expanded.""", 4)
    output += "\n"
    output += "\n"

    # Write reST for all functions
    for (classname, parent, comment, function_documentation) in documentation:

        # Document class
        if classname is not None:

            # Write class documentation
            output += ".. cpp:class:: %s\n" % classname
            output += "\n"

            # Write header if any
            if parent is not None:
                output += indent("*Parent class*\n", 4)
                output += "\n"
                output += indent("* :cpp:class:`%s`\n" % parent, 8)
                output += "\n"

            # Write class documentation
            if comment is not None:
                comment = add_links(comment, classnames, ":cpp:class:")
                output += indent(comment, 4)
                output += "\n"
                output += "\n"

        # Iterate over class functions
        for (signature, comment) in function_documentation:

            # Handle indentation of signature
            if "\n" in signature:
                lines = signature.split("\n")
                signature = "\n".join([lines[0]] + [(len(".. cpp.function::") + 1)*" " + l for l in lines[1:]])

            # Add crosslinks in comment
            comment = add_links(comment, classnames, ":cpp:class:")

            # Write function documentation
            output += indent(".. cpp:function:: %s\n" % signature, 4)
            output += "\n"
            output += indent(comment, 8)
            output += "\n"
            output += "\n"

    # Write file
    f = open(outfile, "w")
    f.write(output)
    f.close()

def generate_index(module, headers):
    "Generate index file for module"

    print "Generating index file for module '%s'...'" % module

    # Sort headers
    headers.sort()

    # Set heading
    heading = "DOLFIN ``%s`` library" % module
    stars = len(heading)*"*"

    # Set contents
    contents = "\n".join(h.split(".h")[0] for h in headers)
    contents = indent(contents, 4)

    # Write top of file
    f = open(os.path.join(output_dir, module, "index.rst"), "w")
    f.write("""\
.. Index file for the %s directory

.. _programmers_reference_cpp_%s_index:

%s
%s
%s

.. toctree::
    :maxdepth: 2

%s
""" % (module, module, stars, heading, stars, contents))

    # Close file
    f.close()

# Get representation and write documentation.
documentation, classnames = extract_doc_representation(dolfin_dir)
for module in documentation:
#    if not module == "quadrature":
#        continue
    headers = []
    for (header, doc) in documentation[module]:
#        if not header == "MeshEntity.h":
#            continue
        write_documentation(doc, header, module, classnames)
        headers.append(header)

    # Generate index
    generate_index(module, headers)


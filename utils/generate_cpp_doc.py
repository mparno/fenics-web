__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2010-08-26"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-08-26

import os

def generate_documentation(filename, module, dolfin_dir):
    "Generate documentation for given filename in given module"

    if not filename == "DirichletBC.h":
        return

    # Extract documentation and sort alphabetically
    documentation = extract_documentation(filename, module, dolfin_dir)
    documentation.sort()

    # Format documentation
    format_documentation(documentation, filename, module)

def extract_documentation(filename, module, dolfin_dir):
    "Extract documentation from given filename in given module"

    print "Generating documentation for %s..." % filename

    # List of (signature, comment)
    documentation = []

    # Comment and signature
    comment = None
    signature = None

    # Indentation of signatures
    indent = 0

    # Iterate over each line
    f = open(os.path.join(dolfin_dir, "dolfin", module, filename))
    for line in f:

        # Look for "///"
        if "///" in line:
            c = line.split("///")[1].strip()

            # Found start of new comment
            if comment is None:
                comment = c

            # Continuing comment on next line
            else:
                comment += "\n" + c

        elif comment is not None:
            s = line.strip()

            # Found start of new signature
            if signature is None:
                signature = s
                indent = (len(s.split("(")[0]) + 1)*" "

            # Continuing signature on next line
            else:
                signature += "\n" + indent + s

            # Signature ends when we find ")"
            if ")" in s:

                # Get function name
                #function = signature.split("(")[0].split(" ")[-1]

                # Store documentation
                documentation.append((signature, comment))

                # Reset comment and signature
                comment = None
                signature = None

    f.close()

    return documentation

def format_documentation(documentation, filename, module):
    "Format documentation for given filename and module"

    # Write reST for all functions
    for (signature, comment) in documentation:

        print indent(".. cpp:function:: %s" % signature, 4)
        print ""
        print indent(comment, 8)
        print ""

def indent(string, num_spaces):
    "Indent given text block given number of spaces"
    return "\n".join(num_spaces*" " + l for l in string.split("\n"))

# Set directory for DOLFIN source code
if not "DOLFIN_DIR" in os.environ:
    raise RuntimeError, "You need to set the DOLFIN_DIR environment variable."
dolfin_dir = os.environ["DOLFIN_DIR"]

# Extract modules from dolfin.h
modules = []
f = open(os.path.join(dolfin_dir, "dolfin", "dolfin.h"))
for line in f:
    if line.startswith("#include <dolfin/"):
        module = line.split("/")[1]
        modules += [module]
f.close()

# Iterate over modules
for module in modules:

    # Extract header files from dolfin_foo.h
    f = open(os.path.join(dolfin_dir, "dolfin", module, "dolfin_%s.h" % module))
    for line in f:

        # Generate documentation for header file
        if line.startswith("#include <dolfin/"):
            filename = line.split("/")[2].split(">")[0]
            generate_documentation(filename, module, dolfin_dir)

    f.close()

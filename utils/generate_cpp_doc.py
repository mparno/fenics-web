__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2010-08-26"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-08-26

import os

def generate_documentation(filename, module, dolfin_dir):
    "Generate documentation for given filename in given module"

    # Extract documentation and sort alphabetically
    documentation = extract_documentation(filename, module, dolfin_dir)

    # Write documentation
    write_documentation(documentation, filename, module)

def extract_documentation(filename, module, dolfin_dir):
    "Extract documentation from given filename in given module"

    print "Generating documentation for %s..." % filename

    # List of classes with documentation
    documentation = []

    # Class name and parent class name
    classname = None
    parent = None

    # Comment and signature
    comment = None
    signature = None

    # Indentation of signatures
    indent = 0

    # Iterate over each line
    f = open(os.path.join(dolfin_dir, "dolfin", module, filename))
    for line in f:

        # Check for comment
        if "///" in line:
            c = line.split("///")[1].strip()

            # Found start of new comment
            if comment is None:
                comment = c

            # Continuing comment on next line
            else:
                comment += "\n" + c

        # Check for class
        elif " class " in line and not ";" in line and not "//" in line:

            # Get class name and parent
            classname = line.split("class")[1].split(":")[0].strip()
            if "public" in line:
                parent = line.split("public")[1].strip()
            else:
                parent = None

            # Store documentation
            documentation.append((classname, parent, comment, []))
            classname = None
            parent = None
            comment = None

        # Check for signature
        elif comment is not None:
            s = line.strip()

            # Found start of new signature
            if signature is None:
                signature = s
                indent = (len(s.split("(")[0]) + 1)*" "

            # Continuing signature on next line
            else:
                signature += "\n" + indent + s

            # Signature ends when we find ";" or "{"
            if ";" in s or "{" in s:

                # Strip out last part
                signature = signature.split(";")[0]
                signature = signature.split("{")[0]
                signature = signature.strip()

                # Get function name
                #function = signature.split("(")[0].split(" ")[-1]

                # Store documentation
                if len(documentation) > 0:
                    documentation[-1][-1].append((signature, comment))
                else:
                    documentation = [(None, None, None, [(signature, comment)])]

                # Reset comment and signature
                comment = None
                signature = None

    # Close file
    f.close()

    # Sort documentation alphabetically within each class
    for (classname, parent, comment, function_documentation) in documentation:
        function_documentation.sort()

    return documentation

def write_documentation(documentation, filename, module):
    "Write documentation for given filename and module"

    # Create containing directory
    directory = os.path.join("source", "programmers-reference", "test", "cpp", module)
    try:
        os.makedirs(directory)
    except:
        pass

    # Set location of documentation
    rstfile = filename.split(".")[0] + ".rst"
    outfile = os.path.join(directory, rstfile)

    # Check that it is ok to overwrite the output file
    try:
        f = open(outfile, "r")
        if "DONT_TOUCH" in f.read():
            print "Not touching file %s" % outfile
            f.close()
            return
    except:
        pass

    # Open output file
    f = open(outfile, "w")

    # Write top of file
    f.write(".. Documentation for the header file dolfin/%s/%s\n" % (module, filename))
    f.write("\n")
    f.write(".. _programmers_reference_cpp_%s_Mesh:\n" % module)
    f.write("\n")
    f.write(filename + "\n")
    f.write(len(filename)*"=" + "\n")
    f.write("\n")

    # Write a note that this file was generated automatically
    f.write(".. note::\n")
    f.write("\n")
    f.write(indent("""\
The documentation on this page was automatically extracted from
the DOLFIN C++ code and needs to be edited and expanded.""", 4))
    f.write("\n")
    f.write("\n")

    # Write reST for all functions
    for (classname, parent, comment, function_documentation) in documentation:

        # Document class
        if classname is not None:

            # Write class documentation
            f.write(".. cpp:class:: %s\n" % classname)
            f.write("\n")

            # Write header if any
            if parent is not None:
                f.write(indent("*Parent class*\n", 4))
                f.write("\n")
                f.write(indent("* :cpp:class:`%s`\n" % parent, 8))
                f.write("\n")

            # Write class documentation
            if comment is not None:
                f.write(indent(comment, 8))
                f.write("\n")
                f.write("\n")

        # Iterate over class functions
        for (signature, comment) in function_documentation:

            # Handle indentation of signature
            if "\n" in signature:
                lines = signature.split("\n")
                signature = "\n".join([lines[0]] + [(len(".. cpp.function::") + 1)*" " + l for l in lines[1:]])

            # Write function documentation
            f.write(indent(".. cpp:function:: %s\n" % signature, 4))
            f.write("\n")
            f.write(indent(comment, 8))
            f.write("\n")
            f.write("\n")

    f.close()

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

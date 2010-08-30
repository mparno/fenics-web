# This script generates the reST input for the C++ Programmer's
# Reference from the DOLFIN C++ header files found in DOLFIN_DIR.
#
# The script looks for comments starting with /// and appearing
# before class and function declarations. The comments are copied
# verbatim to the generate reST files.
#
# Ths script also replaces _Foo_ with links to class pages.

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2010-08-26"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-08-30

import os, re

# Set output directory
output_dir = os.path.join("source", "programmers-reference", "cpp")

# Set directory for DOLFIN source code
if not "DOLFIN_DIR" in os.environ:
    raise RuntimeError, "You need to set the DOLFIN_DIR environment variable."
dolfin_dir = os.environ["DOLFIN_DIR"]

def generate_documentation(header, module):
    "Extract documentation for given header in given module"

    print "Generating documentation for %s..." % header

    # List of classes with documentation
    classnames = []
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
    f = open(os.path.join(dolfin_dir, "dolfin", module, header))
    for line in f:

        # Check for comment
        if "///" in line:

            # We may have either "///" and "/// "
            if "/// " in line:
                c = line.split("/// ")[1].rstrip()
            else:
                c = line.split("///")[1].rstrip()

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
            classnames.append(classname)
            documentation.append((classname, parent, comment, []))
            classname = None
            parent = None
            comment = None

        # Check for function signature
        elif comment is not None:
            s = line.strip()

            # Found start of new signature
            if signature is None:
                signature = s
                #indent = (len(s.split("(")[0]) + 1)*" "

            # Continuing signature on next line
            else:
                #signature += "\n" + indent + s
                signature += " " + s

            # Signature ends when we find ";" or "{"
            if ";" in s or "{" in s:

                # Strip out last part
                signature = signature.split(";")[0]
                signature = signature.split("{")[0]
                signature = signature.strip()

                # Remove stuff Spinx can't handle
                signature = signature.replace("virtual ", "")
                signature = signature.replace("inline ", "")

                # Remove ": stuff" for constructors
                if " : " in signature:
                    signature = signature.split(" : ")[0]

                # Skip destructors (not handled by Sphinx)
                destructor = "~" in signature

                # Get function name
                #function = signature.split("(")[0].split(" ")[-1]

                # Store documentation
                if len(documentation) > 0 and not destructor:
                    documentation[-1][-1].append((signature, comment))
                elif not destructor:
                    documentation = [(None, None, None, [(signature, comment)])]

                # Reset comment and signature
                comment = None
                signature = None

    # Close file
    f.close()

    # Sort documentation alphabetically within each class
    for (classname, parent, comment, function_documentation) in documentation:
        function_documentation.sort()

    return documentation, classnames

def write_documentation(documentation, header, module, classnames):
    "Write documentation for given header in given module"

    # For quick testing
    #if not header == "FunctionSpace.h":
    #    return

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
                comment = add_links(comment, classnames)
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
            comment = add_links(comment, classnames)

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

def add_links(text, classnames):
    "Add crosslinks for classes"
    for classname in classnames:
        p = re.compile(r"\b_%s_\b" % classname)
        text = p.sub(":cpp:class:`%s`" % classname, text)

    return text

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

def indent(string, num_spaces):
    "Indent given text block given number of spaces"
    return "\n".join(num_spaces*" " + l for l in string.split("\n"))

# Extract modules from dolfin.h
modules = []
f = open(os.path.join(dolfin_dir, "dolfin", "dolfin.h"))
for line in f:
    if line.startswith("#include <dolfin/"):
        module = line.split("/")[1]
        modules += [module]
f.close()

# Iterate over modules
documentation = {}
classnames = []
for module in modules:

    # Extract header files from dolfin_foo.h
    f = open(os.path.join(dolfin_dir, "dolfin", module, "dolfin_%s.h" % module))
    documentation[module] = []
    headers = []
    for line in f:

        # Generate documentation for header file
        if line.startswith("#include <dolfin/"):
            header = line.split("/")[2].split(">")[0]
            headers.append(header)
            doc, cls = generate_documentation(header, module)
            documentation[module].append((header, doc))
            classnames += cls

    # Generate index
    generate_index(module, headers)

# Write documentation
for module in documentation:
    for (header, doc) in documentation[module]:
        write_documentation(doc, header, module, classnames)

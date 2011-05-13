__authors__ = "Anders Logg <logg@simula.no>, Kristian Oelgaard <k.b.oelgaard@gmail.com> and Marie E. Rognes <meg@simula.no>"
__copyright__ = "Copyright (C) 2010 " + __authors__
__license__  = "GNU LGPL version 3 or any later version"

import os

# Import the documentation extractor tools from dolfin_utils
from dolfin_utils.documentation import extract_doc_representation
from dolfin_utils.documentation import indent, add_links

import indextemplates as templates

__all__ = ["generate_cpp_api_documentation"]

def write_class_doc(classname, parents, comment, function_documentation,
                    classnames):

    output = []
    if classname is not None:

        output += [".. cpp:class:: %s\n" % classname]

        if parents is not None:
            output += [indent("*Parent class(es)*\n", 4)]
            for parent in parents:
                output += [indent("* :cpp:class:`%s`\n" % parent, 8)]

        if comment is not None:
            comment = add_links(comment, classnames, ":cpp:class:") #NB
            output += [indent(comment, 4)] + ["\n"]

    else:
        output += ["\n"]

    for (signature, comment) in function_documentation:

        # Handle indentation of signature # Hm?
        if "\n" in signature:
            lines = signature.split("\n")
            signature = ("\n".join([lines[0]]
                                   + [(len(".. cpp.function::") + 1)*" " + l
                                      for l in lines[1:]]))
        # Add crosslinks in comment
        comment = add_links(comment, classnames, ":cpp:class:")

        output += [indent(".. cpp:function:: %s\n" % signature, 4)]
        output += [indent(comment, 8), "\n"]

    text = "\n".join(output)
    return text

def generate_header_rst_file(module, header, documentation, classnames):

    output = []
    for (classname, parents, comment, function_documentation) in documentation:
        output += [write_class_doc(classname, parents, comment,
                                   function_documentation, classnames)]

    text = "\n".join(output)
    note = indent(templates.note_template, 4)

    dictionary = {"module": module, "header": header,
                  "prefix": header.split(".")[0].lower(),
                  "title": header, "formatting": len(header)*"=",
                  "note": note, "text": text}
    template = templates.header_cpp_template % dictionary

    return template

def open_write_close(filename, contents):
    f = open(filename, "w")
    f.write(contents)
    f.close()

def generate_cpp_module(module, documentation, classnames, output_dir):

    print indent("Generating documentation for %s" % module, 0)

    # Generate index file for this module
    module_index = templates.generate_module_index(module)
    index_filename = os.path.join(output_dir, "index.rst")
    open_write_close(index_filename, module_index)

    # Generate .rst for each header file in this module
    for (header, doc) in documentation:

        print indent("Generating documentation for %s" % header, 2)
        rst_filename = os.path.join(output_dir, header.split(".")[0] + ".rst")
        header_rst = generate_header_rst_file(module, header, doc,
                                              classnames)
        open_write_close(rst_filename, header_rst)

    return

def generate_cpp_api_documentation(version, input_dir, output_dir):
    """
    Generate .rst files for C++ DOLFIN Programmer's Reference

    Arguments:
      * version
        version should be the version number, example 0.9.11 or dev
    """

    print "Generating .rst for C++ API documentation of DOLFIN-%s" % version
    print "Reading from %s" % input_dir
    print "Outputting to %s" % output_dir
    try:
        os.makedirs(output_dir)
    except:
        print "Dirs already existing %s" % output_dir

    # Get representation
    documentation, classnames = extract_doc_representation(input_dir)

    # Generate index file for library
    library_index = templates.generate_library_index(version)
    library_index_filename = os.path.join(output_dir, "index.rst")
    open_write_close(library_index_filename, library_index)

    # Generate documentation for each module
    for module in documentation:

        # Create directory for this module
        module_dir = os.path.join(output_dir, module)
        try:
            os.makedirs(module_dir)
        except:
            pass

        # Create documentation for this module
        generate_cpp_module(module, documentation[module], classnames,
                            module_dir)



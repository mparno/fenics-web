"""This script generates the reST input for the Python Programmer's
Reference from the DOLFIN Python module."""

__author__ = "Kristian B. Oelgaard <k.b.oelgaard@gmail.com>"
__date__ = "2010-09-15"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-10-08

import os, sys, types

# Set output directory
output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),\
                          os.pardir, "source", "doc/dolfin/programmers-reference", "python")
#output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "python-source")


# Import the dolfin and dolfindocstrings modules.
try:
    import dolfin
    from dolfin_utils.documentation import indent
except Exception as what:
    raise ImportError("Could not import the dolfin module \n  (error: %s),\n\
  update your PYTHONPATH variable?" % what)

index_string = \
"""
%s:

.. toctree::
    :maxdepth: 1

%s
"""

def get_modules(parent, loc, modules):
#    print "parent: ", parent
#    print "dir: ", loc
    for mod in os.listdir(loc):
        f = os.path.join(loc, mod)
        # Add modules (files) to global dict and to parent as submodules.
        if os.path.isfile(f):
            m, e = os.path.splitext(mod)
            if e == ".py" and m != "__init__":
#                print "mod: ", m
                new_mod = ".".join([parent, m])
                if new_mod in modules:
                    print new_mod, modules.keys()
                    raise RuntimeError("module already present???")
                modules[new_mod] = []
                modules[parent].append(m)
        # Add submodules (directories with '__init__.py' files) to global dict
        # and to parent as submodules.
        if os.path.isdir(f):
            if not "__init__.py" in os.listdir(f):
                continue
#            print "f: ", f
            new_mod = ".".join([parent, mod])
            if new_mod in modules:
                print new_mod, modules.keys()
                raise RuntimeError("module already present???")
            modules[new_mod] = []
            modules[parent].append(mod)

            # Recursively extract submodules.
            get_modules(new_mod, f, modules)

#def get_modules(mod, modules, top_module=""):
#    """Extract all modules defined in a module.

#    This function will not return external modules which are imported. To get
#    all modules, the function is called recursively."""
#    # This is the first call to the function, store name of module.
#    if top_module == "":
#        top_module = mod.__name__
#    print "\nmod: ", mod.__name__
#    # Get all sub modules.
#    for k, v in mod.__dict__.items():
#        if not isinstance(v, types.ModuleType):
#            continue
#        print "k: ", k
#        n = v.__name__
#        print "n: ", n
#        if n.split(".")[0] != top_module:
#            continue
#        # To avoid infinite recursion.
#        if n in modules:
#            continue
#        modules[n] = v
#        modules.update(get_modules(v, modules, top_module))
#    return modules

def get_objects(module, submodules):
    """Extract classes and functions defined in a module.
    The function will not return imported classes and functions."""

    modules = []
    # Special handling of 'cpp' module.
    for sub in submodules:
        if sub == "cpp":
            sub = "cpp (Swig autogenerated module) <cpp/index>"
            modules.append(sub)
        else:
            modules.append(sub + "/index")

    classes = []
    functions = []
    d = {}
    exec("from %s import *" % module, d)
    for key, val in d.items():
#        print "key: ", key
        if isinstance(val, (types.ClassType, types.TypeType)):
#            print "cls, mod: ", val.__module__
            if module == val.__module__:
                classes.append(key)
        elif isinstance(val, types.FunctionType):
#            print "fun, mod: ", val.__module__
            if module == val.__module__:
                functions.append(key)
        # Anything else we need to catch?
        else:
            pass

    return modules, classes, functions

def write_object(directory, module_name, name, obj_type):

    output = ".. Documentation for the %s %s\n\n" % (obj_type, module_name + "." + name)
    output += ".. _programmers_reference_python_%s:\n\n" % "_".join(module_name.split(".")[1:] + [name.lower()])
    output += name + "\n"
    output += "="*len(name) + "\n"
    output += "\n.. currentmodule:: %s\n\n" % module_name
    output += ".. auto%s:: %s\n" % (obj_type, name)
    outfile = os.path.join(directory, name + ".rst")
    f = open(outfile, "w")
    f.write(output)
    f.close()

def write_documentation(module, submodules):
    dirs = [output_dir]
    dirs += module.split(".")[1:]
    directory = os.path.sep.join(dirs)

    try:
        os.makedirs(directory)
    except:
        pass

    modules, classes, functions = get_objects(module, submodules)

    output = ".. Index file for the %s module.\n\n" % module
    output += ".. _programmers_reference_python_%s:\n\n" % "_".join(module.split(".")[1:] + ["index"])
    if module == "dolfin":
        output += """#############################
Python programmer's reference
#############################\n"""
    else:
        header = "%s module" % module
        stars = "*"*len(header)
        output += stars + "\n"
        output += header + "\n"
        output += stars + "\n"

    outfile = os.path.join(directory, "index.rst")
    f = open(outfile, "w")
    f.write(output)
    if modules:
        f.write(index_string % ("Modules", indent("\n".join(sorted(modules)), 4)))
    if classes:
        f.write(index_string % ("Classes", indent("\n".join(sorted(classes)), 4)))
    if functions:
        f.write(index_string % ("Functions", indent("\n".join(sorted(functions)), 4)))
    f.close()

    for o in classes:
#        if not o == "Mesh":
#            continue
        write_object(directory, module, o, "class")

    for o in functions:
        write_object(directory, module, o, "function")

#modules = get_modules(dolfin, {})
modules = {"dolfin":[]}
get_modules("dolfin", os.path.dirname(dolfin.__file__), modules)
for mod in sorted(modules.keys()):
    print "Writing files for module: ", mod
#    if not mod == "dolfin.mesh":
#    if not mod == "dolfin":
#        continue
    write_documentation(mod, modules[mod])


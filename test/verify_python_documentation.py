#!/usr/bin/env python
"""This script performs three different sanity tests on the docstrings found in
the Python module 'dolfin' from FEniCS.

Test 1.

  Verify that classes and functions in the dolfin module and dolfindocstrings
  module have the same docstrings. If not, it means that the doctrings in
  DOLFIN have to be regenerated.

Test 2.

  Verify that all classes and functions in the dolfin module are documented in
  the dolfindocstrings module. This test will ensure that we at least have some
  documentation for every feature available in the dolfin module.

Test 3.

  Verify that all classes and functions in the dolfindocstrings module are
  present in the dolfin module. This test will ensure that features which are
  deleted in dolfin will also be deleted from the documentation.

At the end of the test a summary will be printed with a list of errors (if any)
and what actions should be taken.

Note that the test performed by this script cannot determine if the
documentation is up to date or valid.
"""

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2010-06-23"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified on 2010-07-24

# Import modules necessary to perform the test
import types
from sys import path as sys_path
from os import chdir, pardir, getcwd, path
import inspect
import types

# Get location of this script and start there.
test_dir = sys_path[0]
chdir(test_dir)

# cd to parent directory and add documentation directory to path.
chdir(pardir)
doc_dir = path.join(getcwd(), "source", "programmers-reference", "python")
sys_path.append(doc_dir)

# Get location for log file
log_file = path.join(test_dir, "python_log")

# Import the dolfin and dolfindocstrings modules.
try:
    import dolfin
except:
    raise ImportError("Could not import the dolfin module, update your PYTHONPATH variable.")
import dolfindocstrings

# -----------------------------------------------------------------------------
# Functions to extract module information
# -----------------------------------------------------------------------------
def get_modules(mod, modules, top_module=""):
    """Extract all modules defined in a module.
  
    This function will not return external modules which are imported. To get
    all modules, the function is called recursively."""
    # This is the first call to the function, store name of module.
    if top_module == "":
        top_module = mod.__name__
    # Get all sub modules.
    for k, v in mod.__dict__.items():
        if not isinstance(v, types.ModuleType):
            continue
        n = v.__name__
        if n.split(".")[0] != top_module:
            continue
        # To avoid infinite recursion.
        if n in modules:
            continue
        modules[n] = v
        modules.update(get_modules(v, modules, top_module))
    return modules

def get_functions(mod, module_name):
    "Extract all functions from a module or a class."
    # Loop all function members and add to dictionary if the name starts with
    # module name. (is this check needed, we get a class which should be part
    # of the module already!)
    functions = {}
    if isinstance(mod, types.ModuleType):
        name = mod.__name__
    else:
        name = mod.__module__ + "." + mod.__name__

    # Get all function types.
    defs = [(k, v) for k, v in mod.__dict__.items() if inspect.isfunction(v)]
    defs.extend([(k, v) for k, v in mod.__dict__.items() if inspect.ismethod(v)])
    defs.extend([(k, v) for k, v in mod.__dict__.items() if inspect.ismemberdescriptor(v)])
    defs.extend([(k, v) for k, v in mod.__dict__.items() if inspect.ismethoddescriptor(v)])

    for k, v in defs:
        # Only add those functions which are actually defined in this module
        #  and not those which are simply imported.
        # NOTE: If we are checking class members don't try this since
        # 'staticmethod' does not know about __module__.
        #
        if module_name is not None and v.__module__ != module_name:
            continue
        n = name + "." + k
        if not n in functions:
            # Skip private member functions
            if k[0] == "_" and k[:2] != "__":
                continue
            functions[n] = v

    return functions

def parse_module(mod):
    "Extract functions and classes from module (recursively to get all classes)"

    # First get all modules.
    module_name = mod.__name__
    defines = get_modules(mod, {}, module_name)

    # Loop all modules
    for key, val in defines.items():
        # Handle all classes.
        for k, v in val.__dict__.items():
            # FIXME: Use new style classes in FEniCS so we can avoid ClassType
            if not isinstance(v, (types.TypeType, types.ClassType)):
                continue
            # If the class is not defined in this module (but imported) skip.
            if v.__module__ != key:
                continue
            n = v.__module__ + "." + v.__name__
            if n.split(".")[0] != module_name:
                continue
            defines[n] = v
            # Get all class member functions.
            defines.update(get_functions(v, None))

        # Get all functions defined in module.
        defines.update(get_functions(val, key))

    new_defines = {}
    for n in defines:
        nn = ".".join(n.split(".")[1:])
        new_defines[nn] = defines[n]

    return new_defines
#    return defines

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Test no. 1.
# -----------------------------------------------------------------------------
# Test if all docstrings from the doc module are equal to docstrings in the
# dolfin module on import (don't compare white space).
def equal_docs(mod, doc):
    """Test if two modules, classes or functions have equal docstrings"""
    non_equal = []
    for key, val in doc.items():
        if key in mod:
            try:
                vs = val.__doc__.split()
                ms = mod[key].__doc__.split()
                if vs != ms:
                    non_equal.append((key, "'%s' != '%s'" % (" ".join(vs), " ".join(ms))))
            except:
                non_equal.append((key, ""))

    return sorted(non_equal)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Test no. 2.
# -----------------------------------------------------------------------------
def missing_docs(mod, doc, excludes):
    "Test if documentation of modules, classes or functions in module is missing in documentation."
    missing = []
    for key, val in mod.items():
        if not key in doc and not key in excludes:
            missing.append(key)
        elif key in doc and not doc[key].__doc__:
            missing.append(key)
    return sorted(missing)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Test no. 3.
# -----------------------------------------------------------------------------
def obsolete_docs(mod, doc):
    "Test if modules, classes or functions in documentation are present in dolfin."
    obsolete = []
    for key, val in doc.items():
        if not key in mod:
            obsolete.append(key)
    return sorted(obsolete)

# -----------------------------------------------------------------------------

if __name__ == "__main__":

    # List of functions to exclude from tests
    excludes = []

    # Get function and class definitions in modules.
    dolfin_defines = parse_module(dolfin)
    docstrings_defines = parse_module(dolfindocstrings)

    # Run tests.
    non_equal = equal_docs(dolfin_defines, docstrings_defines)
    missing   = missing_docs(dolfin_defines, docstrings_defines, excludes)
    obsolete  = obsolete_docs(dolfin_defines, docstrings_defines)

    print
    print "*"*40 + " TEST SUMMARY " + "*"*40
    log_string = ""
    if non_equal:
        num_docs = len(docstrings_defines)
        s =  "\n    Test no. 1 failed."
        s += "\n    %d out of %d docstrings differ." % (len(non_equal), num_docs)
        print s + "\n    See python_log for details."
        s += "\n    Docstrings differ for the following modules, classes and functions:"
        s += "\n"
        s += "\n".join(["\n".join([" "*6 + c, " "*6 + v]) for c, v in non_equal])
        s += "\n"
        log_string += s
    if missing:
        num_docs = len(dolfin_defines) - len(excludes)
        s =  "\n    Test no. 2 failed."
        s += "\n    %d missing docstrings out of %d." % (len(missing), num_docs)
        print s + "\n    See python_log for details."
        s += "\n    The following modules, classes and functions lacks documentation:"
        s += "\n"
        s += "\n".join(["      " + c for c in missing])
        s += "\n"
        log_string += s
    if obsolete:
        num_docs = len(docstrings_defines)
        s =  "\n    Test no. 3 failed."
        s += "\n    %d out of %d docstrings are obsolete." % (len(obsolete), num_docs)
        print s + "\n    See python_log for details."
        s += "\n    Documentation for the following modules, classes and functions should be deleted:"
        s += "\n"
        s += "\n".join(["      " + c for c in obsolete])
        s += "\n"
        log_string += s

    if non_equal == missing == obsolete == []:
        print
        print "    All tests are OK!"
        print
    else:
        print
        # Write log
        f = open(log_file, "w")
        f.write(log_string)
        f.close()

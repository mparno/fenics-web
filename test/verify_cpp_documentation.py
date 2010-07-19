#!/usr/bin/env python
"""This script performs two different sanity tests on the documentation for
the C++ interface to the FEniCS project DOLFIN.

Test 1.

  Verify that documentation exists for classes and functions in DOLFIN.
  The script will only check if DOLFIN header files have matching reST files
  and that the reST files contain documentation for classes and functions found
  in the header file.
  This test will ensure that we at least have some documentation for all
  classes and functions in DOLFIN.

Test 2.

  Verify that all files, classes and functions in the reST files are present in
  DOLFIN. This test will ensure that features which are deleted in DOLFIN will
  also be deleted from the documentation.

At the end of the test a summary will be printed with a list of errors (if any)
and what actions should be taken.

Note that the test performed by this script cannot determine if the
documentation is up to date or valid.

NOTE: Does not find typedefs, enum and ignores namespaces.
"""

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2010-07-01"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified on 2010-07-19

# Import modules necessary to perform the test.
from os import path, chdir, system, listdir, curdir, pardir, getcwd, mkdir
from shutil import copytree, rmtree
from sys import path as sys_path
from commands import getstatusoutput
import xml.dom.minidom

# Get dolfin directory and add 'dolfin' (not interested in dolfin.h only contents in sub dirs).
error, dolfin_dir = getstatusoutput("pkg-config --variable=includedir dolfin")
dolfin_dir = path.join(dolfin_dir, "dolfin")
if error:
    raise ImportError("Could not locate DOLFIN, update your PKG_CONFIG_PATH variable.")

# Get location of this script and start there.
test_dir = sys_path[0]
chdir(test_dir)

# cd to parent directory and get documentation directory.
chdir(pardir)
doc_dir = path.join(getcwd(), "source", "programmers-reference", "cpp")

extensions = {"C++":".h", "reST":".rst", "XML":".xml"}
log_file = path.join(test_dir, "cpp_log")
# Create log file
f = open(log_file, "w")
f.close()

# -----------------------------------------------------------------------------
# Auxiliary functions for getting file names.
# -----------------------------------------------------------------------------
def get_files(file_type, directory, excludes):
    "Get source files recursively from a given path."
    files = []
    chdir(directory)
    for d in listdir(curdir):
        if path.isdir(d) and not path.islink(d) and not d in excludes:
            sub_dir = get_files(file_type, d, excludes)
            files.extend([path.join(directory, f) for f in sub_dir])
        name, ext = path.splitext(d)
        if ext == extensions[file_type] and not name in excludes:
            files.append(path.join(directory, name))
    chdir(pardir)
    return files

def print_files(files, header, top_dir=""):
    "Pretty print a list of files (strings)."
    if not files:
        return
    print "\nPrinting %d %s files:" % (len(files), header)
    for f in sorted(files):
        print "  " + f.replace(path.join(top_dir, ""), "")

# -----------------------------------------------------------------------------
# Auxiliary classes for C++ documentation.
# -----------------------------------------------------------------------------
class File(object):
    "Representation of a file object, containing classes and/or functions."
    def __init__(self, name, classes, members, finalized=False):
        self.name = name
        self.members = dict([(str(m), m) for m in members])
        self.classes = dict([(c.name, c) for c in classes])
        self.is_finalized = finalized

    def update_class(self, class_):
        "Update the classes defined in a file."
        if class_.name in self.classes:
            raise RuntimeError("Not expected.")
            self.classes[class_.name].update(class_)
        else:
            self.classes[class_.name] = class_

    def update_members(self, members):
        "Update the function definitions in a file."
        for member in members:
            if str(member) in self.members:
                raise RuntimeError("Should not happen.")
            else:
                self.members[str(member)] = member

    def empty(self):
        "Simple check if the file contained any definitions of interest."
        if not self.members and not self.classes:
            return True
        return False

    def finalize(self):
        """Nest all child classes and delete from class dict.
        Needed when parsing the Doxygen xml files."""
        if self.is_finalized:
            return
        for k,v in self.classes.items():
            new_inner = {}
            for i in v.classes:
                # This has to work, otherwise a class has disappeared!
                if i in self.classes:
                    new_inner[i] = self.classes[i]
                    del self.classes[i]
                else:
                    print "\nUnable to extract innerclass: %s from file: %s" % (i, self.name)
            v.classes = new_inner
        self.is_finalized = True

    def num_classes_func(self):
        "Return number of classes and functions defined/documented in file."
        self.finalize()
        num_clas = len(self.classes)
        num_func = len(self.members)
        for k, v in self.classes.items():
            nc, nf = v.num_classes_func()
            num_clas += nc
            num_func += nf
        return (num_clas, num_func)


    def __str__(self):
        "Pretty print for debugging."
        self.finalize()
        s = "\nfile: %s" % self.name
        nc, nf = self.num_classes_func()
        s += ("\nClasses: %d" % nc) + (", functions: %d" % nf)
        if self.members:
            s += "\n\n  members:\n"
            m = "\n".join(sorted(self.members.keys()))
            s += indent(m, 4)
        if self.classes:
            s += "\n\n  classes:\n"
            c = "\n\n".join(sorted([str(v) for k, v in self.classes.items()]))
            s += indent(c, 4)
            s += "\n"
        return s

class Class(object):
    "Representation of a C++ class containing innerclasses and/or member functions."
    def __init__(self, name, members, innerclasses):
        self.name = name
        self.members = dict([(str(m), m) for m in members])
        self.classes = innerclasses

    def num_classes_func(self):
        "Return number of classes and functions nested in class."
        num_clas = len(self.classes)
        num_func = len(self.members)
        for k, v in self.classes.items():
            nc, nf = v.num_classes_func()
            num_clas += nc
            num_func += nf
        return (num_clas, num_func)

    def __str__(self):
        "Pretty print for debugging."
        s = "class: %s" % self.name
        if self.members:
            s += "\n  members:\n"
            m = "\n".join(sorted(self.members.keys()))
            if self.classes:
                m += "\n"
            s += indent(m, 4)
        if self.classes:
            s += "\n  innerclasses:\n"
            c = "\n\n".join(sorted([str(v) for k, v in self.classes.items()]))
            s += indent(c, 4)
        return s

class Member(object):
    "Representation of a C++ (member-) function."
    def __init__(self, return_type, name, args_string):
        self.return_type = return_type
        self.name = name
        self.args_string = args_string
    def __str__(self):
        "Pretty print for debugging."
        if self.return_type:
            return "%s %s%s" % (self.return_type, self.name, self.args_string)
        return "%s%s" % (self.name, self.args_string)

def indent(s, ind):
    return "\n".join([" "*ind + l for l in s.split("\n")])

def format_string(s):
    "Some simple formatting because programmers might NOT write code like I do."
    # Put space between commas.
    s = ", ".join(s.split(","))
    # Simple hack to deal with reference and pointer formatting.
    if " &" in s:
        s = s.replace(" &", "& ")
    if " *" in s:
        s = s.replace(" *", "* ")
    # Simple hack to deal stuff like std::vector<double> --> std::vector< double >
    if "<" in s:
        s = s.replace("<", "< ")
    if ">" in s:
        s = s.replace(">", " >")

    # Remove any whitespace that was created.
    s = " ".join(s.split())
    return s

# -----------------------------------------------------------------------------
# Auxiliary functions for parsing Doxygen generated XML files.
# -----------------------------------------------------------------------------
def handleDocument(doc, files, top_dir):
    "Doxygen XML compounds: file, namespace, class and struct."
    for cdef in doc.getElementsByTagName("compounddef"):
        if cdef.getAttribute("kind") == "file":
            handleFile(cdef, files, top_dir)
        if cdef.getAttribute("kind") == "namespace":
            handleNamespace(cdef, files, top_dir)
        if cdef.getAttribute("kind") in ("class", "struct"):
            handleClass(cdef, files, top_dir)

def handleFile(element, files, top_dir):
    name = get_text(element.getElementsByTagName("compoundname")[0])
    file_name = get_location_file_name(element, top_dir)
    if name != path.split(file_name)[1]:
        raise RuntimeError("File names are different?!")
    members = []
    for mdef in element.getElementsByTagName("memberdef"):
        if mdef.getAttribute("kind") == "function":
            members.append(handleMember(mdef))
    file_name = file_name.replace(extensions["C++"], "")
    if file_name in files:
        files[file_name].update_members(members)
    else:
        files[file_name] = File(file_name, [], members)

def handleNamespace(element, files, top_dir):
    # Ignoring namespace for now, treating like file with no namespace.
    name = get_text(element.getElementsByTagName("compoundname")[0])
    name_space_members = {}
    for mdef in element.getElementsByTagName("memberdef"):
        if mdef.getAttribute("kind") == "function":
            file_name = get_location_file_name(mdef, top_dir)
            if file_name in name_space_members:
                name_space_members[file_name].append(handleMember(mdef))
            else:
                name_space_members[file_name] = [handleMember(mdef)]
    for file_name, members in name_space_members.items():
        file_name = file_name.replace(extensions["C++"], "")
        if file_name in files:
            files[file_name].update_members(members)
        else:
            files[file_name] = File(file_name, [], members)

def handleClass(element, files, top_dir):
    name = get_text(element.getElementsByTagName("compoundname")[0])
    # For now, strip namespaces until we know how we want it in Sphinx.
    # It will require quite some work to extract this as well and I don't think
    # we should include it in the documentation anyway.
    name = name.split("::")[-1]
    innerclasses = [get_text(icls) for icls in element.getElementsByTagName("innerclass")]
    innerclasses = [c.split("::")[-1] for c in innerclasses]
    members = []
    for mdef in element.getElementsByTagName("memberdef"):
        if mdef.getAttribute("kind") == "function":
            members.append(handleMember(mdef))
    class_ = Class(name, members, innerclasses)

    file_name = get_location_file_name(element, top_dir)
    file_name = file_name.replace(extensions["C++"], "")
    if file_name in files:
        files[file_name].update_class(class_)
    else:
        files[file_name] = File(file_name, [class_], [])

def handleMember(element):
    # We can also extract attributes like, 'inline', 'explicit' should we
    # decide that we need it in the documentation.
    return_type = get_text(element.getElementsByTagName("type")[0])
    return_type = format_string(return_type)
    name = get_text(element.getElementsByTagName("name")[0])
    args_string = get_text(element.getElementsByTagName("argsstring")[0])
    args_string = format_string(args_string)
    return Member(return_type, name, args_string)

def get_text(node):
    "Auxiliary function to get text recursively from a node."
    t = node.TEXT_NODE
    if node.nodeType == t:
        return node.nodeValue
    if node.hasChildNodes():
        s = ""
        for cn in node.childNodes:
            if cn.nodeType == t:
                s += cn.nodeValue
            else:
                s += get_text(cn)
        return s
    return ""

def get_location_file_name(element, top_dir=""):
    # Get all file names and make sure that they are the same, then pick one.
    file_names = [n.getAttribute("file") for n in element.getElementsByTagName("location")]
    if not all([i == j for i, j in zip(file_names[:-1], file_names[1:])]):
        raise RuntimeError("All file names should be the same")
    return file_names[0].replace(path.join(top_dir, ""), "")

def get_dolfin_classes(dolfin_files):
    """Extract classes and functions from DOLFIN header files and report empty files.

    Since parsing C++ turned out to be rather difficult, we rely on Doxygen to
    do the hard work for us to generate XML files which are in turn parsed by
    the xml.dom.minidom module.

    The function return a dictionary with relative file names (w.r.t. the
    dolfin directory struture) as keys and the values are File objects."""

    # cd to test directory
    chdir(test_dir)
    # Make sure to delete any old stuff.
    tmp = "tmp"
    if path.exists(tmp):
        rmtree(tmp)

    # Copy dolfin header files to tmp directory and run Doxygen.
    copytree(dolfin_dir, tmp)
    error, output = getstatusoutput("doxygen")
    if error:
        print output
        raise RuntimeError("Error while running Doxygen")

    # Get XML files and process.
    files = {}
    xml_files = get_files("XML", tmp, [])
    for f in xml_files:
        handleDocument(xml.dom.minidom.parse(f + extensions["XML"]), files, path.join(test_dir, tmp))

    # Log info
    log_body = test_output(dolfin_dir, dolfin_files, files)

    # Remove files which are not handled by some definition.
    remove_files(files)

    # Get log string and write
    log_string = format_files(files, "Output from DOLFIN files", log_body)
    f = open(log_file, "a")
    f.write(log_string)
    f.close()

    return files

# -----------------------------------------------------------------------------
# Auxiliary functions for parsing documetation reST files.
# -----------------------------------------------------------------------------
def handle_rst_Class(line_num, lines):
    """Get class name and extract members, handle child classes recursively.
    Relies heavily on strict reST formatting rules."""

    indent, line = lines[0]
    # Remove '.. cpp:class:: from line, the class name MUST be the last word!'
    name = line.split()[-1].strip()

    # Loop lines, add members and handle classes recursively.
    members = []
    classes = []
    local_num = 0
    for e, (i, l) in enumerate(lines[1:]):
        # Continue if the current line number already have been handled by
        # sub-processes.
        if e < local_num:
            continue
        # Break if indentation becomes identical to the level that the class
        # itself was defined on.
        if i <= indent:
            break
        # Add functions to file members.
        if "cpp:function::" in l:
            members.append(handle_rst_Member(l))
        if "cpp:class::" in l:
            d_num, class_ = handle_rst_Class(e, lines[e+1:])
            # Increase local_num so we skip the lines that the innerclass has handled.
            local_num += d_num
            classes.append(class_)
            continue
        # Increase the local line number since we've already handled this line
        local_num += 1
    # Increase the line number since we've handles the class line itself.
    local_num += 1
    return local_num, Class(name, members, dict([(c.name, c) for c in classes]))

def handle_rst_Member(line):
    """Get return type, member name and argument list from line.
    Relies heavily on strict reST formatting rules."""
    # Remove '.. cpp:function::' from line, '..' and 'cpp:function::' MUST be the
    # two first words.
    definition = " ".join(line.split()[2:])
    type_name, args_string = definition.split("(")
    args_string = "(" + args_string
    args_string = format_string(args_string)
    type_name = type_name.split()
    # The name is ALWAYS the last word before ()
    name = type_name[-1]
    return_type = " ".join(type_name[:-1])
    return_type = format_string(return_type)
    # TODO: Could crash program here if we don't have a return type
    # (does not apply to class constructors)
    return Member(return_type, name, args_string)


def get_doc_classes(doc_files):
    """Extract classes and functions from documentation reST files.

    This function will only look for the reST
    .. cpp:class:: and .. cpp:function:: directives and use the indentation
    level to determine which functions and classes belong where.

    The function return a dictionary with relative file names (w.r.t. the
    documentation directory struture) as keys and the values are File objects."""

    files = {}
    for f in doc_files:
        # Read contents of reST file
        rst_file = open(f + extensions["reST"], "r")
        lines = rst_file.read().split("\n")
        rst_file.close()
        # Remove the doc directory from file name to be able to compare names
        # with DOLFIN files directly.
        f = f.replace(path.join(doc_dir, ""), "")

        # Only get interesting lines and determine their indentation level.
        lines = [l for l in lines if "cpp:class::" in l or "cpp:function::" in l]
        indent = map(lambda l: len(l) - len(l.lstrip()), lines)
        indent_lines = zip(indent, lines)

        # Loop interesting lines, if func add to file members, if class call extract class
        line_num = 0
        members = []
        classes = []
        for e, (i, l) in enumerate(indent_lines):
            # Continue if the current line number already have been handled by
            # sub-processes.
            if e < line_num:
                continue
            # Add functions to file members.
            if "cpp:function::" in l:
                members.append(handle_rst_Member(l))
            if "cpp:class::" in l:
                d_num, class_ = handle_rst_Class(e, indent_lines[e:])
                # Increase local_num so we skip the lines that the innerclass has handled.
                line_num += d_num
                classes.append(class_)

        # It should be safe to add a File for every doc file name.
        files[f] = File(f, classes, members, True)

    # Log info
    log_body = test_output(doc_dir, doc_files, files)

    # Remove files which are not handled by some definition.
    remove_files(files)

    # Get log string and write
    log_string = format_files(files, "Output from Documentation files", log_body)
    f = open(log_file, "a")
    f.write(log_string)
    f.close()

    return files

# -----------------------------------------------------------------------------
# Auxiliary functions for output.
# -----------------------------------------------------------------------------
def remove_files(files):
    """Remove files that contain no classes or functions."""
    # TODO: Might want to simply delete the empty constructor and destructor since
    # they are not (yet??!) possible to implement in Sphinx.
    # Also need to handle inline functions, explicit etc.
    for k, v in files.items():
        if v.empty():
            del files[k]

def file_stats(files):
    "Extract statistics from a dictionary of files."

    num_files = len(files)
    num_classes = 0
    num_functions = 0
    for k in sorted(files.keys()):
        v = files[k]
        nc, nf = v.num_classes_func()
        num_classes += nc
        num_functions += nf
    return (num_files, num_classes, num_functions)

def format_files(files, header=None, body=None):
    "Create a string representation of a dictionary of files."

    s = "\n"
    # Get stats and add.
    num_files, num_classes, num_functions = file_stats(files)
    s += "\nTotal number of files: %d\
          \nTotal number of classes: %d\
          \nTotal number of functions: %d" % (num_files, num_classes, num_functions)

    if not body is None:
        s += body

    # Get string representation of file
    for k in sorted(files.keys()):
        s += str(files[k])

    # Prepend header and add footer.
    if not header is None:
        i = len(header)
        s = "\n"+"-"*(i+20)\
          + "\n" + " "*(10) + "%s" % header\
          + "\n" + "-"*(i+20)\
          + s\
          + "\n" + "-"*(i+20)
    s += "\n"
    return s

def test_output(top_dir, file_names, files):
    "Common test output function for reST and C++ headers."

    log_string = ""

    # Log any files that did not contain any classes and functions
    empty_files = [v.name for k, v in files.items() if v.empty()]
    if empty_files:
        log_string += "\nThe following %d C++ header files did not contain any classes or functions:\n" % len(empty_files)
        log_string += "\n".join(empty_files)
        log_string += "\n"

    # Check that we found all files.
    s0 = set(files.keys())
    s1 = set([f.replace(path.join(top_dir, ""), "") for f in file_names])
    if s0 != s1:
        log_string += "\nFiles disappeared:\n  file_from_classes - file_names:"
        log_string += "\n".join(list(s0 - s1))
        log_string += "\n  file_names - file_from_classes:"
        log_string += "\n".join(list(s1 - s0))
        log_string += "\n"

    if log_string == "":
        return None

# -----------------------------------------------------------------------------
# Test no. 1 and 2.
# -----------------------------------------------------------------------------
def missing_classes(object0, object1):
    """Test if classes and functions defined in object0 are present in object1
    where object# can be File or Class objects."""
    classes = []
    members = []
    for k0, v0 in object0.classes.items():
        if k0 not in object1.classes:
            classes.append(v0)
        # Test classes recursively.
        else:
            mc = missing_classes(v0, object1.classes[k0])
            if mc is not None:
                classes.append(mc)

    # Member functions are simple.
    for k0, v0 in object0.members.items():
        if k0 not in object1.members:
            members.append(v0)
    if classes or members:
        if isinstance(object0, File):
            return File(object0.name, classes, members, True)
        elif isinstance(object0, Class):
            return Class(object0.name, members, dict([(str(c), c) for c in classes]))
        else:
            print repr(object0)
            raise RuntimeError("This should be a Class.")
    return None

def test_missing(files0, files1, header):
    """Test if classes and functions defined in files0 are present in files1."""

    missing_files = {}
    for k0, v0 in files0.items():
        if not k0 in files1:
            missing_files[k0] = v0
        else:
            mc = missing_classes(v0, files1[k0])
            if mc is not None:
                missing_files[k0] = mc
    
    # Write log string to file
    f = open(log_file, "a")
    f.write(format_files(missing_files, header))
    f.close()

    return missing_files

# -----------------------------------------------------------------------------

if __name__ == "__main__":

    # Files and directories which should be excluded from tests
    dolfin_excludes = []
    doc_excludes = ["index"]

    # Get DOLFIN header files and extract classes and functions.
    dolfin_files = get_files("C++", dolfin_dir, dolfin_excludes)
    #print_files(dolfin_files, "DOLFIN", dolfin_dir)
    dolfin_classes = get_dolfin_classes(dolfin_files)


    # Get documentation files.
    doc_files = get_files("reST", doc_dir, doc_excludes)
    #print_files(doc_files, "reST", doc_dir)
    doc_classes = get_doc_classes(doc_files)

    # Test 1.
    missing_docs = test_missing(dolfin_classes, doc_classes,\
               "Documentation is missing for the following files and classes.")

    # Test 2.
    obsolete_docs  = test_missing(doc_classes, dolfin_classes,\
              "Documentation is obsolete for the following files and classes.")

    print
    print "*"*40 + " TEST SUMMARY " + "*"*40
    print

    if missing_docs:
        i, j, k = file_stats(missing_docs)
        print "    Test no. 1 failed."
        print "    Found %d files with %d classes and %d functions in DOLFIN that lacks documentation." % (i,j,k)
        print "    See %s for details\n" % log_file

    if obsolete_docs:
        i, j, k = file_stats(obsolete_docs)
        print "    Test no. 2 failed."
        print "    Found %d files with %d classes and %d functions in the documentation which are obsolete." % (i,j,k)
        print "    See %s for details\n" % log_file

    if missing_docs == obsolete_docs == {}:
        print "    All tests are OK! (see cpp_log for details)"
        print



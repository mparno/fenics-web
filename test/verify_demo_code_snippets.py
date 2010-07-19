#!/usr/bin/env python

"""This utility script will find all *.rst files in the source/demo
directory and checks that any code snippets highlighted by the .. code-block::
directive is legal in the sense that it is present in at least one of the
source files (.ufl, .py, .cpp) that is associated with the demo."""

from os import chdir, path, getcwd, curdir, pardir, listdir
from sys import stderr, path as sys_path

# Make sure we start where this test script is located.
chdir(sys_path[0])

# We currently only verify demo code.
chdir(path.join(pardir, "source", "demos"))

# We have C++ and Python versions of the demos.
directories = ["cpp", "python"]

# Dictionary of code blocks that has to be checked for each subdirectory
# including information about file types of the source.
block_source =  {"cpp":     {"c++": ".cpp", "python": ".ufl"},
                 "python":  {"python": ".py"}
                }

def verify_blocks(rst_file, source_files, source_dict):
    """Check that any code blocks in the rst file is present in at least one of
    the source files."""

    for block_type, source_type in source_dict.items():
        # Extract code blocks from rst file.
        blocks = get_blocks(rst_file, block_type)
        for block in blocks:
            # Check if block is in the list of files of correct type.
            block_in_source(block, [sf for sf in source_files\
                                    if path.splitext(sf)[-1] == source_type])

def get_blocks(rst_file, block_type):
    "Extract any code blocks of given type from the rst file."

    blocks = []

    # Open file and read lines.
    f = open(rst_file, "r")
    lines = f.read().split("\n")

    code = False
    block = []
    for l in lines:
        # Get start of code block.
        if "code-block::" in l and block_type in l:
            code = True
            block = []
            continue
        # The first line which is not an indented line terminates the code
        # block.
        if code and l and l[0] != " ":
            code = False
            # Join the block that we have and add to list of blocks.
            # Remove any whitespace.
            blocks.append(remove_whitespace("\n".join(block)))
        # If code is still True, then the line is part of the code block.
        if code:
            block.append(l)

    # Add block of code if found at the end of the rst file.
    if code:
        blocks.append(remove_whitespace("\n".join(block)))

    # Close file and return blocks.
    f.close()
    return blocks

def remove_whitespace(code):
    "Remove blank lines and whitespace in front of lines."
    return "\n".join([" ".join(l.split())\
                      for l in code.split("\n") if l != ""])

def block_in_source(block, source_files):
    """Check that the code block is present in at least one of
    the source files."""

    present = False
    code = ""
    for sf in source_files:
        f = open(sf, "r")
        # Read code and remove whitespace before comparing block and code.
        code = remove_whitespace(f.read())

        if block in code:
            present = True
        f.close()

        # If code is present, look no further.
        if present:
            return

    # Just crash the test, no need to proceed.
    if not present:
        if not source_files:
            print "\ncode block:\n", block
            raise RuntimeError("No source file!")

        print "\nError:"
        print "\ncode block:\n", block
        print "\nsource_files:\n", source_files
        print "\nin directory: ", getcwd()
        print
        raise RuntimeError("Illegal code block.")

if __name__ == "__main__":

    print "\nTesting that all code snippets are valid.\n"
    # Loop directories/categories/demos
    for directory in directories:
        chdir(directory)
        # Get all demo categories
        categories = [d for d in listdir(curdir) if path.isdir(d)]
        for category in categories:
            chdir(category)
            demos = [d for d in listdir(curdir) if path.isdir(d)]
            for demo in demos:
                chdir(demo)
                stderr.write(" "*2 + path.join(directory, category, demo))
                # Get files in demo directory and sort in rst and source files.
                files = listdir(curdir)
                rst_files = [f for f in files if path.splitext(f)[-1] == ".rst"]
                source_files = [f for f in files if path.splitext(f)[-1] in\
                                  (".py", ".ufl", ".cpp")]
                # Loop files, check if code blocks are present in source files.
                for rst_file in rst_files:
                    verify_blocks(rst_file, source_files, block_source[directory])
                stderr.write(", " + "OK.\n")
                chdir(pardir)
            chdir(pardir)
        chdir(pardir)

    print "\nOK.\n"


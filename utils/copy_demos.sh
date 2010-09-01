#!/bin/sh
#
# This script copies all demos from the DOLFIN demo directory as
# specified by DOLFIN_DIR. The script will try to only copy the
# relevant files, ignoring for example data output and files generated
# by CMake which may be present in the demo directories.
#
# Copyright (C) 2010 Anders Logg
# Licensed under the GNU GPL version 3 or any later version
#
# First added:  2010-09-01
# Last changed: 2010-09-01

# Set destination directory
DEST_DIR="source/demos"

# Set directories that should be copied
DEMO_DIRS="la pde undocumented"

# Check DOLFIN_DIR
if [ x"$DOLFIN_DIR" = "x" ]; then
    echo "You need to set the DOLFIN_DIR environment variable."
    exit 1
fi

# Copy files over using rsync
echo "Copying demos from $DOLFIN_DIR/demo"
for DEMO_DIR in $DEMO_DIRS; do

    echo "------------------------------------------------------------"
    echo "Copying demos for directory '$DEMO_DIR'"
    echo

    rsync -av \
        --exclude .vtu \
        --exclude .pvd \
        --exclude CMakeFiles \
        --exclude cmake_install.cmake \
        $DOLFIN_DIR/demo/$DEMO_DIR $DEST_DIR

done

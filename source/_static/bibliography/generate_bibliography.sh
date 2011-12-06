#!/bin/sh
#
# This script generates the bibliography in various formats
# from the common input file 'fenics.pub'.

# Check that we have publish
PUBLISH=`which publish`
if [ "x$PUBLISH" = "x" ]; then
    echo "The Publish program is required to run this script."
    echo "This may be installed by running the following commands:"
    echo ""
    echo "  hg clone ssh://hg@bitbucket.org/logg/publish"
    echo "  cd publish"
    echo "  sudo python setup.py install"
    echo ""
    exit 1
fi

# Export the two files
publish export fenics.bib database_filename=fenics.pub debug=1
publish export fenics.rst database_filename=fenics.pub debug=1

# Rename .rst file to get rid of Sphinx warning
mv fenics.rst fenics.rst_

# Copy file to citing directory
cp fenics.rst_ ../../citing

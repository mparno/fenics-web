#!/bin/sh

chmod +x build/html/index.html
rsync -avz --delete build/html/ fenics@fenicsproject.org:/home/fenics/www.fenicsproject.org/new/

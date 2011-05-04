#!/bin/sh

rsync -avz build/html/ fenics@fenicsproject.org:/home/fenics/www.fenicsproject.org/new/

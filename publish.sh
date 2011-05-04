#!/bin/sh

rsync -avz --delete build/html/ fenics@fenicsproject.org:/home/fenics/www.fenicsproject.org/new/

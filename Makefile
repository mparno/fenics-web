# Parent Makefile for Sphinx documentation, for convenience.
#

# Directories where Sphinx documentation has to be generated.
DIRS	= demos fenics-doc tutorial user-manual

.PHONY: help clean html dirhtml pickle json htmlhelp qthelp latex changes linkcheck doctest

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  dirhtml   to make HTML files named index.html in directories"
	@echo "  pickle    to make pickle files"
	@echo "  json      to make JSON files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  qthelp    to make HTML files and a qthelp project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  changes   to make an overview of all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  doctest   to run all doctests embedded in the documentation (if enabled)"

clean:
	rm -f *~
	-for d in $(DIRS); do (cd $$d; make clean; ); done

html:
	-for d in $(DIRS); do (cd $$d; make html; ); done

latex:
	-for d in $(DIRS); do (cd $$d; make latex; ); done

changes:
	-for d in $(DIRS); do (cd $$d; make changes; ); done

linkcheck:
	-for d in $(DIRS); do (cd $$d; make linkcheck; ); done

doctest:
	-for d in $(DIRS); do (cd $$d; make doctest; ); done


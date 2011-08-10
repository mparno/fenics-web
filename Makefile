# Makefile for FEniCS documentation

# You can set these variables from the command line
SPHINXOPTS  =
SPHINXBUILD = sphinx-build
PAPER       =

# Internal variables
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d build/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) source

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo ""
	@echo "  all       to build all targets: html, doc, publish"
	@echo "  web       to build web pages"
	@echo "  doc       to build (import) documentation for projects"
	@echo "  publish   to publish everything on www.fenicsproject.org"
	@echo "  clean     to clean out everything (build directory)"
	@echo ""
	@echo "In addition, the targets 'latex' and 'pdf' exist but are not used."

all:	web doc publish

web:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) build/html
	@echo
	@echo "Build finished. HTML generated in build/html."

doc:
	scripts/import_docs

publish:
	scripts/publish

clean:
	-rm -rf build

latex:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) build/latex
	@echo
	@echo "Build finished. LaTeX generated in build/latex."

pdf:
	make -C build/latex all-pdf

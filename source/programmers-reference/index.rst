.. The programmer's reference for FEniCS goes here.

.. _programmers_reference_index:

#############################
FEniCS Programmer's reference
#############################

From the Launchpad `blueprint
<https://blueprints.launchpad.net/fenics-doc/+spec/user-manual>`_:

* The FEniCS User Manual should be hand-written in reST and converted to HTML
  and PDF using Sphinx.

* Relevant parts should also be available as docstrings in the Python
  interface. This should be possible to fix with some clever scripting
  (assigning to __doc__).

* It should come in two different flavors, C++ and Python. As much as possible,
  material should be reused between the two.

* It should replace the current user manuals for DOLFIN, FFC, UFL and UFC.

* Emphasis should be on documenting the DOLFIN user interface. Details for FFC,
  UFL and UFC should be placed in Appendix.

* It should include numerous small examples (code snippets) that illustrate
  each function.

* It should contain a HOWTO chapter that answers questions like how do I set
  Neumann boundary conditions, how do I compute norms, how do I compute errors
  etc.

Contents:

.. toctree::
    :maxdepth: 1

    cpp/index
    python/index
    appendices/index

The user manuals are also available in PDF format, see
:download:`C++ version<../../build/latex/cpp_programmers_reference.pdf>` or
:download:`Python version<../../build/latex/python_programmers_reference.pdf>`.

Build PDFs for FFC, UFC and UFL manuals and link to those too.



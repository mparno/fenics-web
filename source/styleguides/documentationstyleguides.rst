.. Style guides for C++ and Python source code

.. _styleguides_sphinx_coding_style:

Sphinx coding style
===================

This style guide contains information on how to create documentation
for FEniCS using the Sphinx documentation tool. The first sections are
related to how the reST code should look. Then, a section on some
frequently used reST and Sphinx markup follows as a quick
reference. The last two sections are guides explaining in steps how to
document the programmer's reference and the demos respectively.

Code layout
-----------

Use spaces instead of tabs for indentation.

Use 4 spaces per indentation level. This does not apply to C++ code
examples (DOLFIN) where the 2 space indentation rule apply
(cf. :ref:`C++ indentation
<styleguides_cpp_coding_style_indentation>`).

The text width of normal text should not exceed 79 characters. This
rule can be disregarded if readability demands it, for instance for
code example tables.


Sections
--------

Section markers follow the convention from the `Python Documentation
<http://docs.python.org/documenting/rest.html>`_:

* ``#`` with overline, for parts
* ``*`` with overline, for chapters
* ``=``, for sections
* ``-``, for subsections
* ``^``, for subsubsections
* ``"``, for paragraphs

.. _styleguides_sphinx_cross_referencing:

Cross referencing
-----------------

Cross-references are created by placing labels at the location you
want to refer to and then using ``:ref:\`label-name\``` to create the
link. Example:

.. code-block:: rest

    .. _my-reference-label:

    Section to cross-reference
    --------------------------

    This is the text of the section.

    It refers to the section itself, see :ref:`my-reference-label`.

For this to work properly, the label names **must** be unique in the
entire documentation source.  To ensure this, the label names should
begin with the path to the file where the label is located relative to
the source directory. For example, the label for the C++ version of
the Poisson demo, which is located at the top of the
``demos/pde/poisson/cpp/documentation.rst`` file, should be given the name
``demos_pde_poisson_cpp_documentation``, while the label to this sub section is:

.. code-block:: rest

    .. _styleguides_sphinx_cross_referencing:

    Cross referencing
    -----------------

Frequently used markup (roles and directives)
---------------------------------------------

.. _styleguides_sphinx_code_snippets:

Code snippets
^^^^^^^^^^^^^

The FEniCS documentation makes heavy use of code snippets to
illustrate how the interfaces work. Code snippets are created using
the ``code-block`` directive (see `showing code examples
<http://sphinx.pocoo.org/markup/code.html>`_ for more details). This
makes it possible to show C++ and Python code snippets in the the
following way:

.. code-block:: rest

    .. code-block:: c++

        for (int i = 0; i < 10; i++)
          std::cout << i << std::endl;

and

.. code-block:: rest

    .. code-block:: python

        for i in range(10):
            print i

which results in the output:

.. code-block:: c++

    for (int i = 0; i < 10; i++)
      std::cout << i << std::endl;

and

.. code-block:: python

    for i in range(10):
        print i

respectively.

Math
^^^^

Writing FEniCS documentation often involves presenting mathematics,
especially when documenting demos. We use the ``math`` role and
directive to display inline math and equations respectively (see `math
support in Sphinx <http://sphinx.pocoo.org/ext/math.html>`_ for more
details).  The input markup for math is LaTeX, so the inline equation,
:math:`f(x) = x^2`, is typeset by

.. code-block:: rest

    :math:`f(x) = x^2`

and the equation

.. math::

    a(v,u) = \int \nabla v \cdot \nabla u \; \rm{d}\Omega

is typeset as

.. code-block:: rest

    .. math::

        a(v,u) = \int \nabla v \cdot \nabla u \; \rm{d}\Omega

.. note::

    You will need the package ``dvipng`` to display the math properly in HTML.

.. _styleguides_sphinx_download_files:

Download files
^^^^^^^^^^^^^^

To make a file available for download, use the ``download`` role (see
`inline markup <http://sphinx.pocoo.org/markup/inline.html>`_ for more
details) in the following way:

.. code-block:: rest

    See the :download:`main.cpp <../../source/main.cpp>` file.

Author comments
^^^^^^^^^^^^^^^

Please refrain from using the keywords *note*, *todo* and *fixme* in
comments such as

.. code-block:: rest

    .. note: Figure out how to present this in a better way
    .. todo: Add more text and equations
    .. fixme: The results look wrong, probably the boundary conditions

If you think a comment is required, use the ``note`` directive:

.. code-block:: rest

    .. note::

        Figure out how to present this in a better way

and ask on the fenics@lists.launchpad.net mailing list in order for
the issue to be resolved as quickly as possible. This helps keeping
the documentation in good shape.

.. _styleguides_sphinx_documenting_interface:

Documenting the FEniCS interface (programmer's reference)
---------------------------------------------------------

This short guide explains how to write documentation for the C++ and
Python interfaces to FEniCS.
Since the Python interface is (partially) generated automatically using
`Swig <http://www.swig.org/>`_ from the C++ implementation of DOLFIN the
directory/file structure of the documentation follows that of the C++ version
of DOLFIN.
In addition, we want the documentation for the Python version to be
available when using FEniCS with the Python interpreter.
To achieve this we write all documentation for the Python version in a
pseudo module which is an exact replication of the 'real' DOLFIN module and
then let the `Sphinx autodoc <http://sphinx.pocoo.org/ext/autodoc.html>`_
extension handle the rest.

To make matters more concrete let's consider the case of writing documentation
for the DOLFIN ``Mesh`` class and the ``closest_cell`` member function of this
class.

The ``Mesh`` class is defined in the file ``dolfin_dir/dolfin/mesh/Mesh.h``.
We therefore start by adding the files:

* ``programmers-reference/cpp/mesh/Mesh.rst``
* ``programmers-reference/python/mesh/Mesh.rst``

and updating the index files

* ``programmers-reference/cpp/index.rst``
* ``programmers-reference/cpp/mesh/index.rst``
* ``programmers-reference/python/index.rst``
* ``programmers-reference/python/mesh/index.rst``

appropriately.

We then proceed to add contents for the two different interfaces as described
in the following sections.

General remarks
^^^^^^^^^^^^^^^

To handle the documentation of two different languages in Sphinx, we
use `Sphinx Domains <http://sphinx.pocoo.org/domains.html>`_ to
distinguish between classes and functions belonging to the C++ and
Python interfaces.

As Sphinx does not allow sections in the markup for class/function
documentation, we use *italics* (``*italics*``) and definition lists
to group information.  The idea is to keep the markup as simple as
possible since the reST source for the Python documentation of classes
and functions will be used 'as is' in the docstrings of the DOLFIN
module.

Most information can be put in the three sections:

* *Arguments*, which are formatted using definition lists following this
  structure::

    *Arguments*
        <name>
            <type, description>
        <name2>
            <type, description>

  example::

      *Arguments*
          dim
              An integer, some dimension.
          d
              A double, some value.

* *Returns*, which is formatted in a similar fashion::

    *Returns*
        <return type>
            <description>

  example::

    *Returns*
        integer
            Some random integer.

* *Example*, a very small code snippet that shows how the
  class/function works. It does not necessarily have to be a
  stand-alone program.

Links to demos that use the feature being documented should be put in
a ``seealso`` directive.

The member functions of a class should be sorted alphabetically for
the C++ version. When using autodoc, Sphinx will sort the member
functions automatically for the Python module.

C++ interface
^^^^^^^^^^^^^^^^^

The code snippets presented in the following can be seen in their
complete form and context by clicking on the ``Show Source`` link on
the page containing the C++ documentation for the :cpp:class:`Mesh`
class.

The C++ documentation for the ``Mesh`` class is added to the
``programmers-reference/cpp/mesh/Mesh.rst`` file.

Defining the class
""""""""""""""""""

The beginning of the ``programmers-reference/cpp/mesh/Mesh.rst`` file
looks as follows:

.. code-block:: rest

    Mesh.h
    ======

    .. cpp:class:: Mesh

        *Parent class*

            * :cpp:class:`Variable`

        A Mesh consists of a set of connected and numbered mesh entities.

where only the first part of the ``Mesh`` class description has been included
for brevity.

We start with a section title ``Mesh.h`` since the ``Mesh.rst`` should
contain documentation for all classes and functions defined in
``Mesh.h`` and there might be multiple classes defined.  The ``Mesh``
class is defined by the Sphinx directive ``cpp:class::`` followed by
the name of the class.  Since the ``Mesh`` class derives from the
``Variable`` class, we list all parent classes explicitly where the
line ``:cpp:class:`Variable``` will create a link to the C++
documentation of the class ``Variable``.

.. note::

    In the future Sphinx might be clever enough to handle parent classes
    automatically, but until then this is how we do it.

Then follows a description of the purpose of the ``Mesh`` class before
the documentation of the member functions.

Constructors
""""""""""""

The constructors are documented as any other member function.
For the ``Mesh`` class we have two additional constructors besides the empty
constructor:

.. code-block:: rest

    .. cpp:class:: Mesh

        [snip]

        .. cpp:function:: Mesh(const Mesh& mesh)

            Copy constructor.

            *Arguments*
                mesh
                    A :cpp:class:`Mesh` object.

        .. cpp:function:: Mesh(std::string filename)

            Create mesh from data file.

            *Arguments*
                filename
                    A string, name of file to load.

The functions are defined in the class body such that they
automatically have the ``Mesh`` namespace.  The signature of the
functions (in this case the constructors ``Mesh(const Mesh& mesh)``
and ``Mesh(std::string filename)``) **must** be identical to that
found in the ``dolfin/mesh/Mesh.h`` file, otherwise subsequent testing
will report that the function is not documented (or obsolete).

.. note::

    It looks like the destructor ``~`` is not recognised, but we can skip
    documenting that until it is included in Sphinx.

    The empty constructor, in this case Mesh(), is implicitly created when
    defining the class (``.. cpp:class:: Mesh``).
    Explicitly defining it as one of the constructors will cause Sphinx to
    complain about multiple definitions.

closest_cell function
"""""""""""""""""""""

The documentation for the ``closest_cell`` function is added in the
same manner as the documentation for the constructors, but with
additional information about the return value and an example.

.. code-block:: rest

    .. cpp:function:: dolfin::uint closest_cell(const Point & point) const

        Computes the index of the cell in the mesh which is closest to the
        point query.

        *Arguments*
            point
                A :cpp:class:`Point` object.

        *Returns*
            integer
                The index of the cell in the mesh which is closest to point.

        *Example*
            .. code-block:: c++

                UnitSquare mesh(1, 1);
                Point point(0.0, 2.0);
                info("%d", mesh.closest_cell(point));

            output::

                1

Again, the function is defined in the class body, and the signature of
the function is identical to that found in the ``dolfin/mesh/Mesh.h``
file.

.. note::

    Since Sphinx does not yet handle overloaded functions that well, links to
    :cpp:func:`Mesh::closest_cell` (``:cpp:func:`Mesh::closest_cell```) from the
    index page will point to the class where it is defined instead of the
    actual function.
    This behaviour will hopefully change in the future.

Python interface
^^^^^^^^^^^^^^^^

The code snippets presented in the following can be seen in their
complete form and context by clicking on the ``Show Source`` link on
the page containing the Python documentation for the
:py:class:`dolfin.cpp.Mesh` class and in the
:download:`../programmers-reference/python/docstrings/dolfin/cpp.py` file
which contains the actual documentation for the Python ``Mesh`` class.

The Python ``Mesh`` class is defined in the :py:mod:`dolfin.cpp`
module. This module is automatically generated by Swig from the DOLFIN
C++ implementation.  In order to have the documentation for the Python
interface available from within the Python interpreter, we will put
the docstrings for the ``Mesh`` class in the appropriate file of the
pseudo module containing the dolfin docstrings namely
``programmers-reference/python/docstrings/dolfin/cpp.py``.

Defining the class
""""""""""""""""""

We simply define the class as it is defined in the 'real' DOLFIN module with
the exception that it only contains docstrings and no actual code.

.. code-block:: rest

    class Mesh(Variable):
        """
        A Mesh consists of a set of connected and numbered mesh entities.

        [snip]
        """

where only the first part of the ``Mesh`` class description has been
included for brevity.  Since the docstrings module requires correct
Python syntax the parent class ``Variable`` must of course be defined
also.

Construtors
"""""""""""

The constructors are documented in the docstring of the ``Mesh.__init__``
function like this:

.. code-block:: rest

    class Mesh(Variable):
        [snip]

        def __init__(self, *args):
            """
            **Overloaded versions**

            * Mesh\ **()**

              Create empty mesh.

            * Mesh\ **(mesh)**

              Copy constructor.

              *Arguments*
                  mesh
                      A :py:class:`Mesh` instance.

            * Mesh\ **(filename)**

              Create mesh from data file.

              *Arguments*
                  filename
                      A string, name of file to load.
            """

Since the constructor is overloaded, we use the argument list ``(self,
*args)`` in the function definition and the ``**Overloaded
versions**`` section to document the overloaded versions in a standard
list using ``*``.  For each constructor we define the argument list
using **bold face**. The ``\`` is needed to avoid adding space between
the function name (``Mesh``) and the argument list.  Each constructor
is then documented as any other function.

.. note::

    The above approach applies to any overloaded function, not just
    constructors and is necessary because Sphinx/Python does not handle
    overloaded functions (which we need to since the Python interface is
    generated from a C++ implementation).

closest_cell function
"""""""""""""""""""""

The documentation for the ``closest_cell`` function is added in the
same manner as the documentation for the constructors, but with
additional information about the return value and an example.

.. code-block:: rest

    class Mesh(Variable):
        [snip]

        def closest_cell(self, point):
            """
            Computes the index of the cell in the mesh which is closest to the
            point query.

            *Arguments*
                point
                    A :py:class:`Point` instance.

            *Returns*
                integer
                    The index of the cell in the mesh which is closest to point.

            *Example*
                >>> mesh = dolfin.UnitSquare(1,1)
                >>> point = dolfin.Point(0.0, 2.0)
                >>> mesh.closest_cell(point)
                1
            """

Using Sphinx autodoc
""""""""""""""""""""

To complete the Python documentation for the ``Mesh`` class, we simply add
the following to the ``programmers-reference/python/mesh/Mesh.rst`` file:

.. code-block:: rest

    Mesh
    ====

    .. currentmodule:: dolfin.cpp

    .. autoclass:: Mesh
        :members:
        :show-inheritance:
        :undoc-members:

We use the file ``programmers-reference/python/mesh/Mesh.rst`` to mirror the
structure of the DOLFIN source tree (see
:ref:`styleguides_sphinx_documenting_interface`).
The ``currentmodule`` directive tells Sphinx in which module to find the class
that should be documented.
The line ``.. autoclass:: Mesh`` automatically generates documentation for the
``Mesh`` class and the arguments specifies what information to be included.

Appendices
^^^^^^^^^^

Documentation for the FFC, UFC and UFL components of FEniCS is located
in the :ref:`appendix <programmers_reference_appendices_index>`.  The
structure of the documentation of a given module depends on the
file/class layout of the module and the content should be extracted
from the docstrings as is done for the Python interface to DOLFIN.
The layout of the docstrings should follow the same rules as outlined
in the above sections.

Testing the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^

When you are done writing documentation for the programmer's reference you
should run the two tests:

* ``test/verify_cpp_documentation.py``
* ``test/verify_python_documentation.py``

and build the documentation to ensure everything is in working order by running
the command:

.. code-block:: sh

    make all

in the top directory.


.. _styleguides_sphinx_documenting_demos:

Documenting demos
-----------------

This short guide explains the procedure for documenting a FEniCS demo.
As an example, we will demonstrate the steps involved to create
documentation for the :ref:`Poisson (C++) <demos_pde_poisson_cpp>` and
:ref:`Poisson (Python) <demos_pde_poisson_python>` demos.

Files
^^^^^

The demo documentation is located in the ``source/demos``
directory. This directory contains sub-directories for the various categories:

1. adaptivity
2. fem
3. function
4. la
5. mesh
6. ode
7. parameters
8. pde
9. plot
10. quadrature

.. warning::

    This might change in case we decide to reorganise the demos!

The Poisson demo mainly demonstrates how to solve a certain partial
differential equation (PDE), so we should add the following files:

``demos/pde/poisson/common.txt``
    Common information should be placed in this file, and the file
    should then be included in the C++ and Python versions (see
    :ref:`styleguides_sphinx_common_information`).

``demos/pde/poisson/cpp/documentation.rst``
    This file contains the reST source file with the documentation that is
    specific to the C++ version of the Poisson demo.

``demos/pde/poisson/cpp/main.cpp``
    This file contains the entire C++ source code for the solver and must be made
    available for :ref:`download <styleguides_sphinx_download_files>`.

``demos/pde/poisson/cpp/Poisson.ufl``
    This file contains the form file and must be made available for
    :ref:`download <styleguides_sphinx_download_files>`.
    If your demo contains multiple form files, all of these must be added.

``demos/pde/poisson/cpp/CMakeLists.txt``
    This file is necessary to compile the demo against DOLFIN and must be
    made available for :ref:`download <styleguides_sphinx_download_files>`.

``demos/pde/poisson/python/documentation.rst``
    This file contains the reST source file with the documentation
    that is specific to the Python version of the Poisson demo.

``demos/pde/poisson/python/demo.py``
    This file contains the entire Python source code for the solver and must
    be made available for :ref:`download
    <styleguides_sphinx_download_files>`.

Finally, add the demo to the index files to complete the setup of files.
This is done by adding the line ``poisson/cpp/documentation`` to the
``toctree`` of the ``demos/pde/index-cpp.rst`` file and the line
``poisson/python/documentation`` to the ``toctree`` of the
``demos/pde/index-python.rst`` file

The source code files should of course compile and run with the
versions of FEniCS software covered by the current documentation.

.. _styleguides_sphinx_common_information:

Common information
^^^^^^^^^^^^^^^^^^

Each demo should be available in both a C++ and a Python version.
However, the summary (describing what features are demonstrated) along
with the problem and method description are typically identical for
both versions.  It is therefore desirable to put this information in a
common source file to avoid code duplication.  This common code is
placed in the file ``demos/pde/poisson/common.txt``, which is
then included in the two files ``demos/pde/poisson/cpp/documentation.rst``
and ``demos/pde/poisson/python/documentation.rst`` using the ``include``
directive with the relative path to the file:

.. code-block:: rest

  .. include:: ../common.txt

C++ and Python specific contents
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each step of the solution procedure of a demo should be
explained. This can often be achieved by including
:ref:`styleguides_sphinx_code_snippets`.

.. note::

    It is important that the code snippets are exact copies of what can be
    found in the source files. This can be checked by the running the script
    ``test/verify_code_snippets.py`` in the top level directory.

As an example, the definition of the Dirichlet boundary is:

.. code-block:: c++

    class DirichletBoundary : public SubDomain
    {
      bool inside(const Array<double>& x, bool on_boundary) const
      {
        return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS;
      }
    };

for the C++ Poisson demo and

.. code-block:: python

    def boundary(x):
        return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

for the Python demo.

Additional information
^^^^^^^^^^^^^^^^^^^^^^

Use the ``note`` and ``warning`` directives to highlight important
information.  The ``seealso`` directive should be used when pointing
to alternative solutions or functions in the
:ref:`programmers_reference_index`.

Keywords should be added to the index, using the ``index`` directive to make
the documentation easier to navigate through.

See `the Sphinx documentation
<http://sphinx.pocoo.org/markup/para.html#index-generating-markup>`_
for how to use the above directives.

Testing the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^

When you are done writing documentation for the demos you should run the test:

* ``test/verify_demo_code_snippets.py``

and build the documentation to ensure everything is in working order by running
the command:

.. code-block:: sh

    make all

in the top directory.

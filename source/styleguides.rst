.. Style guides for C++, Python, and documentation



DOLFIN (C++) and FEniCS documentation (Sphinx).

.. _styleguides:

************
Style guides
************

.. _styleguides_cpp_coding_style:

C++ coding style
================

Naming conventions
------------------

Class names
^^^^^^^^^^^
Use camel caps for class names:

.. code-block:: c++

    class FooBar
    {
      ...
    };

Function names
^^^^^^^^^^^^^^

Use lower-case for function names and underscore to separate words:

.. code-block:: c++

    foo();
    bar();
    foo_bar(...);

Functions returning a value should be given the name of that value,
for example:

.. code-block:: c++

    class Array:
    {
    public:

      /// Return size of array (number of entries)
      uint size() const;

    };

In the above example, the function should be named ``size`` rather
than ``get_size``. On the other hand, a function not returning a
value but rather taking a variable (by reference) and assigning a value
to it, should use the ``get_foo`` naming scheme, for example:

.. code-block:: c++

    class Parameters:
    {
    public:

      /// Retrieve all parameter keys
      void get_parameter_keys(std::vector<std::string>& parameter_keys) const;

    };


Variable names
^^^^^^^^^^^^^^

Use lower-case for variable names and underscore to separate words:

.. code-block:: c++

    Foo foo;
    Bar bar;
    FooBar foo_bar;

Enum variables and constants
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Enum variables should be lower-case with underscore to separate words:

.. code-block:: c++

    enum Type {foo, bar, foo_bar};

We try to avoid using ``#define`` to define constants, but when
necessary constants should be capitalized:

.. code-block:: c++

    #define FOO 3.14159265358979

File names
^^^^^^^^^^

Use camel caps for file names if they contain the
declaration/definition of a class. Header files should have the
suffix ``.h`` and implementation files should have the
suffix ``.cpp``:

.. code-block:: c++

    FooBar.h
    FooBar.cpp

Use lower-case for file names that contain utilities/functions (not
classes).

Miscellaneous
-------------

.. _styleguides_cpp_coding_style_indentation:

Indentation
^^^^^^^^^^^

Indentation should be two spaces and it should be spaces. Do **not**
use tab(s).

Comments
^^^^^^^^

Comment your code, and do it often. Capitalize the first letter and
don't use punctuation (unless the comment runs over several
sentences). Here's a good example from ``TopologyComputation.cpp``:

.. code-block:: c++

    // Check if connectivity has already been computed
    if (connectivity.size() > 0)
      return;

    // Invalidate ordering
    mesh._ordered = false;

    // Compute entities if they don't exist
    if (topology.size(d0) == 0)
      compute_entities(mesh, d0);
    if (topology.size(d1) == 0)
      compute_entities(mesh, d1);

    // Check if connectivity still needs to be computed
    if (connectivity.size() > 0)
      return;

    ...

Integers and reals
^^^^^^^^^^^^^^^^^^

Use ``dolfin::uint`` instead of ``int`` (unless you really want to use
negative integers). Use ``dolfin::real`` instead of ``double`` only if
you are sure that you want to exploit arbitray precision:

.. code-block:: c++

    uint i = 0;
    double x = 0.0;

Placement of brackets
^^^^^^^^^^^^^^^^^^^^^

Curly brackets following multi-line control statements should appear
on the next line and should not be indented:

.. code-block:: c++

    for (uint i = 0; i < 10; i++)
    {
      ...
    }

For one line statements, ommit the brackets:

.. code-block:: c++

    for (uint i = 0; i < 10; i++)
      foo(i);

Header file layout
^^^^^^^^^^^^^^^^^^

Header files should follow the below template:

.. code-block:: c++

    // Copyright (C) 2008 Foo Bar.
    // Licensed under the GNU LGPL Version 2.1.
    //
    // Modified by Bar Foo, 2008.
    //
    // First added:  2008-01-01
    // Last changed: 2008-02-01

    #ifndef __FOO_H
    #define __FOO_H

    namespace dolfin
    {

      class Bar; // Forward declarations here

      /// Documentation of class

      class Foo
      {
      public:

        ...

      private:

        ...

      };

    }

    #endif

Implementation file layout
^^^^^^^^^^^^^^^^^^^^^^^^^^

Implementation files should follow the below template:

.. code-block:: c++

    // Copyright (C) 2008 Foo Bar.
    // Licensed under the GNU LGPL Version 2.1.
    //
    // Modified by Bar Foo, 2008.
    //
    // First added:  2008-01-01
    // Last changed: 2008-02-01

    #include <dolfin/Foo.h>

    using namespace dolfin;

    //-----------------------------------------------------------------------------
    Foo::Foo() : // variable initialization here
    {
      ...
    }
    //-----------------------------------------------------------------------------
    Foo::~Foo()
    {
      // Do nothing
    }
    //-----------------------------------------------------------------------------

The horizontal lines above (including the slashes) should be exactly 79
characters wide.

Including header files and using forward declarations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Do not use ``#include <dolfin.h>`` or ``#include``
``<dolfin/dolfin_foo.h>`` inside the DOLFIN source tree. Only include
the portions of DOLFIN you are actually using.

Include as few header files a possible and use forward declarations
whenever possible (in header files). Put the ``#include`` in the
implementation file.  This reduces compilation times and minimizes the
risk for cyclic dependencies.

Explicit constructors
^^^^^^^^^^^^^^^^^^^^^

Make all one argument constructors (except copy constructors) explicit:

.. code-block:: c++

    class Foo
    {
      explicit Foo(uint i);
    };

Virtual functions
^^^^^^^^^^^^^^^^^

Always declare inherited virtual functions as virtual in the subclasses.
This makes it easier to spot which functions are virtual.

.. code-block:: c++

    class Foo
    {
      virtual void foo();
      virtual void bar() = 0;
    };

    class Bar : public Foo
    {
      virtual void foo();
      virtual void bar();
    };

Use of libraries
----------------

Prefer ``C++`` strings and streams over old ``C``-style ``char*``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use ``std::string`` instead of ``const char*`` and use ``std::istream`` and
``std::ostream`` instead of ``FILE``. Avoid ``printf``,
``sprintf`` and other C functions.

There are some exceptions to this rule where we need to use old ``C``-style
function calls. One such exception is handling of command-line arguments
(``char* argv[]``).

Python coding style
===================

The FEniCS coding style for Python code adheres to the `PEP-8 style
guide <http://www.python.org/dev/peps/pep-0008/>`_ although it is not
strictly enforced.

.. _styleguides_sphinx_coding_style:

Sphinx coding style for FEniCS documentation
============================================

This style guide contains information on how to create documentation for
FEniCS using the Sphinx documentation tool. The first sections are related
to how the reST code should look like, then follows a section on frequently
used reST and Sphinx markup as a quick reference. The last two sections
are guides explaining in steps how to document the programmer's reference
and the demos respectively.

Code layout
-----------

Use spaces instead of tabs for indentation.

Use 4 spaces per indentation level. This does not apply to ``C++`` code
examples (DOLFIN) where the 2 space indentation rule apply.
See :ref:`C++ indentation <styleguides_cpp_coding_style_indentation>`.

The text width of normal text should not exceed 79 characters, but code example
tables and other cases where readability demands it this rule can be
disregarded.


Sections
--------

Section markers follow the convention from
`Python <http://docs.python.org/documenting/rest.html>`_:

* ``#`` with overline, for parts
* ``*`` with overline, for chapters
* ``=``, for sections
* ``-``, for subsections
* ``^``, for subsubsections
* ``"``, for paragraphs

.. _styleguides_sphinx_cross_referencing:

Cross referencing
-----------------

Cross-references are created by placing labels at the location which you want
to refer to and then use ``:ref:\`label-name\``` to create the link. Example:

.. code-block:: rest

    .. _my-reference-label:

    Section to cross-reference
    --------------------------

    This is the text of the section.

    It refers to the section itself, see :ref:`my-reference-label`.

For this to work properly, the label names **must** be unique in the entire
documentation source.  To ensure this, the label names should begin with
the path to the file where the label is located relative to the source
directory. As an example the label for the ``C++`` version of the Poisson
demo which is located at the top of the ``demos/cpp/pde/poisson/poisson.rst``
file should be given the name ``demos_cpp_pde_poisson`` while the label to
this sub section is:

.. code-block:: rest

    .. _styleguides_sphinx_cross_referencing:

    Cross referencing
    -----------------

Frequently used markup (roles and directives)
---------------------------------------------

.. _styleguides_sphinx_code_snippets:

Code snippets
^^^^^^^^^^^^^

The FEniCS documentation makes heavy use of code snippets to illustrate how the
interfaces work. Code snippets are created using the ``code-block`` directive
(see `showing code examples <http://sphinx.pocoo.org/markup/code.html>`_
for more details) which make it possible to show ``C++`` and ``Python`` code
snippets in the the following way:

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

Writing FEniCS documentation often involves presenting mathematics especially
when documenting demos. We use the ``math`` role and directive to display
inline math and equations respectively (see
`math support in Sphinx <http://sphinx.pocoo.org/ext/math.html>`_ for more
details).
The input markup for math is LaTeX so the inline equation,
:math:`f(x) = x^2`, is typeset by

.. code-block:: rest

    :math:`f(x) = x^2`

and the equation

.. math::

    a(v,u) = \int \nabla v \cdot \nabla u \; \rm{d}\Omega

is implemented as:

.. code-block:: rest

    .. math::

        a(v,u) = \int \nabla v \cdot \nabla u \; \rm{d}\Omega

.. note::

    You will need the package ``dvipng`` to display the math properly in HTML.

.. _styleguides_sphinx_download_files:

Download files
^^^^^^^^^^^^^^

To make a file available for download use the ``download`` role (see
`inline markup <http://sphinx.pocoo.org/markup/inline.html>`_ for more
details) in the following way:

.. code-block:: rest

    See the :download:`main.cpp <../../source/main.cpp>` file.

Author comments
^^^^^^^^^^^^^^^

Please refrain from using the keywords *note*, *todo* and *fixme* in comments
like:

.. code-block:: rest

    .. note: Figure out how to present this in a better way
    .. todo: Add more text and equations
    .. fixme: The results look wrong, probably the boundary conditions

If you feel a comment is in its place use the ``note`` directive:

.. code-block:: rest

    .. note::

        Figure out how to present this in a better way

and ask on the fenics@lists.launchpad.net mailing list, so we can resolve the
issue as quickly as possible to keep the documentation in good shape.

.. _styleguides_sphinx_documenting_interface:

Documenting the FEniCS interface (programmer's reference)
---------------------------------------------------------

This short guide explains how to write documentation for the ``C++`` and
``Python`` interfaces to FEniCS.
Since the ``Python`` interface is (partially) generated automatically using
Swig from the ``C++`` implementation of DOLFIN the directory/file structure of
the documentation follows that of the ``C++`` version of DOLFIN.
In addition, we want the documentation for the ``Python`` version to be
available when using FEniCS with the ``Python`` interpreter.
To achieve this we write all documentation for the ``Python`` version in a
pseudo module which is an exact replication of the 'real' DOLFIN module and
then let the `Sphinx autodoc <http://sphinx.pocoo.org/ext/autodoc.html>`_
extension handle the rest.

To make matters more concrete let's consider the case of writing documentation
for the DOLFIN ``Mesh`` class and the ``closest_cell`` member function of this
class.

The ``Mesh`` class is defined in the file ``dofin_dir/dolfin/mesh/Mesh.h``.
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

To handle the documentation of two different languages in Sphinx we use
`Sphinx Domains <http://sphinx.pocoo.org/domains.html>`_ to distinguish between
classes and functions belonging to the ``C++`` and ``Python`` interfaces.

Since Spinx does not allow sections in the markup for class/function
documentation we use *italics* (``*italics*``) and definition lists to group
information.
The idea is to keep the markup as simple as possible since the reST source for
the ``Python`` documentation of classes and functions will be used 'as is' in
the docstrings of the DOLFIN module.

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

* *Example*, a very small code snippet which shows how the class/function
  works. It does not necessarily have to be a stand-alone program.

Links to demos which use the feature being documented should be put in a
``seealso`` directive.

The member functions of a class should be sorted alphabetically for the
``C++`` version.
When using autodoc, Sphinx will sort the member functions automatically for the
``Python`` module.

``C++`` interface
^^^^^^^^^^^^^^^^^

The code snippets presented in the following can be seen in their complete
form and context by clicking on ``Show Source`` link on the page containing
the ``C++`` documentation for the :cpp:class:`Mesh` class.

The ``C++`` documentation for the ``Mesh`` class is added to the
``programmers-reference/cpp/mesh/Mesh.rst`` file.

Defining the class
""""""""""""""""""

The begining of the ``programmers-reference/cpp/mesh/Mesh.rst`` file looks
like this:

.. code-block:: rest

    Mesh.h
    ======

    .. cpp:class:: Mesh

        *Parent class*

            * :cpp:class:`Variable`

        A Mesh consists of a set of connected and numbered mesh entities.

where only the first part of the ``Mesh`` class description has been included
for brevity.

We start with a section title ``Mesh.h`` since the ``Mesh.rst`` should contain
documentation for all classes and functions defined in ``Mesh.h`` and there
might be multiple classes defined.
The ``Mesh`` class is defined by the Sphinx directive ``cpp:class::`` followed
by the name of the class.
Since the ``Mesh`` class derives from the ``Variable`` class we list all parent
classes explicitly where the line ``:cpp:class:`Variable``` will create a link
to the ``C++`` documentation of the class ``Variable``.

.. note::

    In the future Sphinx might be clever enough to handle parent classes
    automatically, but until then this is how we do it.

Then follows a description of the purpose of the ``Mesh`` class before th
documentation of the member functions.

Construtors
"""""""""""

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

The funtions are defined in the class body such that they automatically have the
``Mesh`` namespace.
The signature of the functions (in this case the constructors
``Mesh(const Mesh& mesh)`` and ``Mesh(std::string filename)``) **must** be
identical to that found in the ``dolfin/mesh/Mesh.h`` file, otherwise
subsequent testing will report that the function is not documented
(or obsolete).

.. note::

    It also looks like the destructor ``~`` is not recognised, but we can skip
    documenting that until it is included in Sphinx.

    The empty constructor, in this case Mesh(), is implicitly created when
    defining the class (``.. cpp:class:: Mesh``).
    Explicitly defining it as one of the constructors will cause Sphinx to
    complain about multiple definitions.

closest_cell function
"""""""""""""""""""""

The documentation for the ``closest_cell`` function is added like documentation
for the constructors with additional information about the return value and an
example.

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

Again, the funtion is defined in the class body, and the signature of the is
identical to that found in the ``dolfin/mesh/Mesh.h`` file.

.. note::

    Since Sphinx does not yet handle overloaded functions that well, links to
    :cpp:func:`Mesh::closest_cell` (``:cpp:func:`Mesh::closest_cell```) from the
    index page will point to the class where it is defined instead of the
    actual function.
    This behaviour will hopefully change in the future.

``Python`` interface
^^^^^^^^^^^^^^^^^^^^

The code snippets presented in the following can be seen in their complete
form and context by clicking on ``Show Source`` link on the page containing
the ``Python`` documentation for the :py:class:`dolfin.cpp.Mesh` class and in the
:download:`programmers-reference/python/docstrings/dolfin/cpp.py` file which
contains the actual documentation for the ``Python`` ``Mesh`` class.

Using Sphinx autodoc
""""""""""""""""""""

To complete the ``Python`` documentation for the ``Mesh`` class, we simply add
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
``Mesh`` class and the arguments

Appendices
^^^^^^^^^^

Documentation for the FFC, UFC and UFL components of FEniCS are located in
the :ref:`appendix <programmers_reference_appendices_index>`.
The structure of the documentation of a given module depends on the file/class
layout of the module and the content should be extracted from the docstrings
as is done for the ``Python`` interface to DOLFIN.
The layout of the docstrings should follow the same rules as outlined in the
above sections.

.. _styleguides_sphinx_documenting_demos:

Documenting demos
-----------------

This short guide explains the procedure for documenting a FEniCS demo.
As an example we will demonstrate the steps involved to create documentation
for the :ref:`Poisson (C++) <demos_cpp_pde_poisson>` and
:ref:`Poisson (Python) <demos_python_pde_poisson>` demos.

Files
^^^^^

The documentation is located in the ``source/demos`` directory which contains
the directories ``common``, ``cpp`` and ``python``.
First you must figure out which category your demo belongs to:

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

In our case the Poisson demo is a partial differential equation (PDE), so
we should add the following files:

``demos/common/pde/poisson/poisson.txt``
    put common information is this file and include in the ``C++`` and
    ``Python`` versions (see :ref:`styleguides_sphinx_common_information`).

``demos/cpp/pde/poisson/poisson.rst``
    this file contains the reST source file with the documentation which is
    specific to the ``C++`` version of the Poisson demo.

``demos/cpp/pde/poisson/main.cpp``
    this file contains the entire source code for the solver and must be made
    available for :ref:`download <styleguides_sphinx_download_files>`.

``demos/cpp/pde/poisson/Poisson.ufl``
    this file contains the form file and must be made available for
    :ref:`download <styleguides_sphinx_download_files>`.
    If your demo contains multiple form files they too must be added.

``demos/cpp/pde/poisson/SConstruct``
    this file is necessary to compile the demo against DOLFIN, and must be
    made available for :ref:`download <styleguides_sphinx_download_files>`.

``demos/python/pde/poisson/poisson.rst``
    this file contains the reST source file with the documentation which is
    specific to the ``Python`` version of the Poisson demo.

``demos/python/pde/poisson/demo.py``
    this file contains the entire solver writte in PyDOLFIN, and must be made
    available for :ref:`download <styleguides_sphinx_download_files>`.

Finally, add the demo to the index files ``demos/python/pde/index`` and
``demos/cpp/pde/index`` by adding ``poisson/poisson`` to the ``toctree`` to
complete the setup of files.

The source code files should of course be able to run with the versions of
FEniCS software e.g., DOLFIN, FFC and UFL, which is covered by the current
documentation.

.. _styleguides_sphinx_common_information:

Common information
^^^^^^^^^^^^^^^^^^

The demo has to be available in a ``C++`` and a ``Python`` version.
However, the summary (describing what features are demonstrated) along with the
problem and method description are typically identical for both versions.
It is therefore desirable to put this information in a common source file to
avoid code duplication.
In this case this code is put in the file
``demos/common/pde/poisson/poisson.txt`` which is then included in the two files
``demos/cpp/pde/poisson/poisson.rst`` and
``demos/python/pde/poisson/poisson.rst`` using the ``include`` directive with
the relative path to the file:

.. code-block:: rest

  .. include:: ../../../common/pde/poisson/poisson.txt

``C++`` and ``Python`` specific contents
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each step of the solution procedure of a demo should be explained. Often this
is achieved by including :ref:`styleguides_sphinx_code_snippets` which of
course must be given in the correct syntax depending on the demo version.

.. warning::

    It is important that the code snippets are exact copies of what can be
    found in the source files. The reason being that the source files will be
    compiled and tested against DOLFIN and if anything is broken the demos
    needs to be updated.

    Running the script ``test/verify_code_snippets.py`` in the top level directory
    will then locate all code snippets that needs to be updated to the new
    syntax.

As an example, the definition of the Dirichlet boundary is:

.. code-block:: c++

    class DirichletBoundary : public SubDomain
    {
      bool inside(const Array<double>& x, bool on_boundary) const
      {
        return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS;
      }
    };

for the ``C++`` Poisson demo and

.. code-block:: python

    def boundary(x):
        return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

for the ``Python`` demo.

Additional information
^^^^^^^^^^^^^^^^^^^^^^

Use the ``note`` and ``warning`` directives to highligt important information.
The ``seealso`` directive should be used when pointing to alternative
solutions or functions in the :ref:`programmers_reference_index`.

Keywords should be added to the index, using the ``index`` directive to make
the documentation easier to navigate through.

See `the Sphinx documentation
<http://sphinx.pocoo.org/markup/para.html#index-generating-markup>`_ on how to
use the above directives.

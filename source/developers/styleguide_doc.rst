.. _developers_styleguide_documentation:

FEniCS documentation coding style guide
=======================================

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

Only the first word of section headings should be capitalised unless, of
course, a word is a name which is normally capitalised. E.g., use

.. code-block:: rest

    FEniCS Python demos
    -------------------

instead of

.. code-block:: rest

    FEniCS Python Demos
    -------------------

.. _styleguides_sphinx_cross_referencing:

Cross referencing
-----------------

Cross-references are created by placing labels at the location you
want to refer to and then using ``:ref:`label-name``` to create the
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

    a(u, v) = \int \nabla u \cdot \nabla v \; \rm{d}\Omega

is typeset as

.. code-block:: rest

    .. math::

        a(u, v) = \int \nabla u \cdot \nabla v \; \rm{d}\Omega

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

.. _styleguides_sphinx_images:

Images
^^^^^^

To include an image, simply use the ``image`` directive (see
`Image <http://docutils.sourceforge.net/docs/ref/rst/directives.html#image>`_
for more details and options) in the following way:

.. code-block:: rest

    .. image:: image.png
        :scale: 75
        :align: center

.. _styleguides_sphinx_documenting_interface:

Documenting the FEniCS interface (programmer's reference)
---------------------------------------------------------

The reST files from which the :ref:`programmers_reference_index` is generated,
are automatically generated for the DOLFIN C++ library and the Python module
``dolfin``.
This short guide explains how to accomplish this and how to format the
documentation.

Design considerations
^^^^^^^^^^^^^^^^^^^^^

The procedure for writing documentation might seem cumbersome so let's have a
look at the design considerations which have led to this ostensible case of
overengineering.

The Python interface is (partially) generated automatically using
`Swig <http://www.swig.org/>`_ from the C++ implementation of DOLFIN.
Some classes are extended when building (see the ``dolfin/dolfin/swig/*post.i``
files) while others are added or extended manually in the Python layer defined
in ``dolfin/site-packages/dolfin``.
While this approach saves a lot of work when implementing the Python interface
it puts some constraints on the way the documentation can be handled.
In addition we have the following design goals for writing and maintaining the
documentation:

Avoid duplicate text
    In many cases the documentation of a feature will be virtually identical
    for the C++ and Python interfaces, and since the Python interface is
    generated from the C++ code, the documentation should be in the C++ source
    code.
    To avoid that the documentation on these pages and the comments in the
    source code (and the implementation itself) diverge, the documentation
    should be automatically generated from the C++ source code.
    Therefore the comments should be written using Sphinx markup.

Help in the Python interpreter
    The documentation of a class/function when running ``help(dolfin.foo)``
    in the Python interpreter should be identical to what can be found online.
    In practice this means that we have to generate the
    ``dolfin/dolfin/swig/docstrings.i`` file using the comments extracted from
    the C++ source before building the Python interface with Swig.

Simple markup
    Since the documentation is written directly in the C++ source code, we want
    markup to be simple such that we have 'code with comments' rather than
    'comments with code'.
    Another reason for preferring simple markup is that it is the raw docstring
    which will be available from the Python interpreter.

General remarks
^^^^^^^^^^^^^^^

As Sphinx does not allow sections in the markup for class/function
documentation, we use *italics* (``*italics*``) and definition lists to group
information.
This is to keep the markup as simple as possible since the reST source for the
Python documentation of classes and functions will be used 'as is' in the
docstrings of the DOLFIN module.

Most information can be put in the three sections:

* *Arguments*, which are formatted using definition lists following this
  structure::

    *Arguments*
        <name> (<type>)
            <description>
        <name2> (<type>)
            <description>

  example::

      *Arguments*
          dim (int)
              some dimension.
          d (double)
              some value.

* *Returns*, which is formatted in a similar fashion::

    *Returns*
        <return type>
            <description>

  example::

      *Returns*
          int
              Some random integer.

* *Example*, a very small code snippet that shows how the
  class/function works. It does not necessarily have to be a
  stand-alone program.

Links to demos that use the feature being documented should be put in
a ``seealso`` directive.

Documenting a feature
^^^^^^^^^^^^^^^^^^^^^

To make matters more concrete let's consider the case of writing documentation
for the member function ``closest_cell`` of the DOLFIN ``Mesh`` class.
The Python interface to this class is generated by Swig and it is not extended
in the Python layer.
Writing documentation for other classes and functions in DOLFIN which are not
extended or added in the Python layer follow a similar procedure.

Adding docstrings to source files
"""""""""""""""""""""""""""""""""

The ``Mesh::closest_cell`` function is defined in the file
``dolfin_dir/dolfin/mesh/Mesh.h``, and the comment lines and function
definition look as follows:

.. code-block:: c++

    /// Computes the index of the cell in the mesh which is closest to the
    /// point query.
    ///
    /// *Arguments*
    ///     point (_Point_)
    ///         A _Point_ object.
    ///
    /// *Returns*
    ///     uint
    ///         The index of the cell in the mesh which is closest to point.
    ///
    /// *Example*
    ///     .. code-block:: c++
    ///
    ///         UnitSquare mesh(1, 1);
    ///         Point point(0.0, 2.0);
    ///         info("%d", mesh.closest_cell(point));
    ///
    ///     output::
    ///
    ///         1
    dolfin::uint closest_cell(const Point& point) const;

Note that the documentation of a function or class is placed above the
definition in the source code.
The structure and content follow the guidelines in the previous section.

The Point object is a class like Mesh and it is defined in the FEniCS interface.
To insert a link to the documentation of this class use leading and trailing
underscore i.e., ``_Point_``.
When parsing the comment lines this string will be substituted with either
``:cpp:class:`Point``` or ``:py:class:`Point``` depending on whether
documentation for the C++ or Python interface is being generated.
The return type, in this case ``dolfin::uint``, will automatically be mapped to
the correct Python type when generating the documentation for the Python
interface.

.. note::

    If you are writing documentation for one of the functions/classes which are
    added to the Python layer manually you have to add manually the correct
    links and types. In the above case ``:py:class:`Point``` and ``int``
    respectively.

The example code uses C++ syntax because it is located in the C++ header file.
Translating this code to a correct Python equivalent is rather difficult.
It is therefore necessary to add example code using the Python syntax manually.
This code should be put in the ``dolfin/dolfin/swig/codeexamples.py`` which
contains a simple dictionary of example code.
The dictionary containing only the example code for the example above should
look as follows:

.. code-block:: python

    codesnippets = {
    "Mesh":{
    "dolfin::uint closest_cell(const Point& point) const":
    """
    .. code-block:: python

        >>> mesh = dolfin.UnitSquare(1, 1)
        >>> point = dolfin.Point(0.0, 2.0)
        >>> mesh.closest_cell(point)
        1
    """}
    }

The first dictionary contains dictionaries for all classes with code examples
for each function.
Note that the full C++ function signature has been used to identify the
function to which the code example belongs.

After adding the documentation to the ``Mesh.h`` file and Python code example
to the ``codeexamples.py`` file, you have to run the script
``dolfin/dolfin/swig/generate.py`` to generate the
``dolfin/dolfin/swig/docstrings.i`` file and then build DOLFIN to update the
docstrings in the ``dolfin`` Python module.

Generating the documentation
""""""""""""""""""""""""""""

To generate the documentation pages for the C++ interface, you need to run the
script ``fenics-doc/utils/generate_cpp_doc.py``.
This will create reST files containing the documentation from all header files
found in DOLFIN.
The contents from ``Mesh.h`` will be saved in the
``programmers-reference/cpp/mesh/Mesh.rst`` file which can be seen in its
complete form and context by clicking on the ``Show Source`` link on the
:cpp:class:`Mesh` class page.
You have to set the ``DOLFIN_DIR`` variable first which should point to the
directory where the DOLFIN version which you want to document is located.

To generate the documentation pages for the Python interface, you need to run
the script ``fenics-doc/utils/generate_python_doc.py``.
The Python ``dolfin`` module has to be in your ``PYTHONPATH`` for this to work.
Since the Python ``Mesh`` class is defined in the :py:mod:`dolfin.cpp`
module generated by Swig the output reST file for the Mesh class is
``programmers-reference/python/cpp/Mesh.rst`` which can be seen by clicking on
the ``Show Source`` link on the :py:class:`dolfin.cpp.Mesh` class page.
This file contain very little markup since we rely on the
`Sphinx autodoc <http://sphinx.pocoo.org/ext/autodoc.html>`_ extension to
extract the Python docstrings automatically.

Finally, build the documentation by running::

    make all

in the ``fenics-doc`` directory.

Summary
^^^^^^^

In summary, to update/generate the documentation follow the below procedure:

* Make appropriate changes to the DOLFIN source code.
* If you made changes to C++ header files or docstrings in
  ``dolfin/dolfin/swig/*.i`` you should update the
  ``dolfin/dolfin/swig/codeexamples.py`` file with an example snippet if
  applicable and run the script ``dolfin/dolfin/swig/generate.py``
  to update the ``dolfin/dolfin/swig/docstrings.i`` file.

* Build DOLFIN to update the ``dolfin`` Python module.
* Update your ``PYTHONPATH`` variable to point to the ``dolfin`` module
  and set the ``DOLFIN_DIR`` variable to point to the DOLFIN directory.
* Run the scripts ``fenics-doc/utils/generate_cpp_doc.py`` and
  ``fenics-doc/utils/generate_python_doc.py``.
* Build the documentation by running::

    make all

  in the top directory.

Python modules
^^^^^^^^^^^^^^

Describe how to write documentation (docstrings and autodoc) for Python modules
UFL, FFC etc.

.. note::

    This section is incomplete because we have not yet started to migrate the
    old manuals yet.

Appendices
^^^^^^^^^^

Documentation for the FFC, UFC and UFL components of FEniCS is located
in the :ref:`appendix <programmers_reference_appendices_index>`.  The
structure of the documentation of a given module depends on the
file/class layout of the module and the content should be extracted
from the docstrings as is done for the Python interface to DOLFIN.
The layout of the docstrings should follow the same rules as outlined
in the above sections.

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

.. note::

    Only the files ``poisson/common.txt``, ``poisson/cpp/documentation.rst``
    and ``poisson/python/documentation.rst`` should be added manually.
    A template with these files is :download:`available <foo.tar>`.

    The source files should be added automatically but running the script
    ``utils/copy_demos.sh`` from the top directory, note that you need to set
    the environment variable ``DOLFIN_DIR`` to the given version of DOLFIN
    which you are documenting.

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

For some demos a picture showing the solution obtained from the input is
appropriate. The picture can be saved by pressing the ``i`` key in the Viper
plot window. This file e.g., ``simulation0000.png``, should then be copied to
the same location as the ``common.txt`` file (but use a different name
for instance ``poisson_u.png`` to avoid name conflicts) and included like this:

.. code-block:: rest

    .. image:: ../poisson_u.png
        :scale: 75
        :align: center

(see :ref:`styleguides_sphinx_images` for more details).
Notice that the parent directory is included in the path although the two files
``common.txt`` and ``poisson_u.png`` are located in the same directory.
This is necessary because the ``common.txt`` file will be included in files in
the two sub directories ``cpp`` and ``python``. The picture should be included
immediately after the input information as seen for instance in the
:ref:`Poisson (C++) <demos_pde_poisson_cpp>` demo.

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

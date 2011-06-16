.. _styleguide_documentation:

FEniCS documentation coding style guide
=======================================

This style guide contains general information on how to write
documentation for FEniCS using the Sphinx documentation tool. This
page mainly contains a section on frequently used reST and Sphinx
markup for reference.

See :ref:`styleguides_sphinx_documenting_interface`
and :ref:`styleguides_sphinx_documenting_demos` for special instructions
regarding the DOLFIN documentation.

General
-------

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
entire documentation source.

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

In code-blocks as above, use spaces instead of tabs for
indentation. For Python code examples, use 4 spaces per indentation
level. For C++ code examples (DOLFIN), the 2 space indentation rule
apply (cf. :ref:`C++ indentation
<styleguides_cpp_coding_style_indentation>`).


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

.. .. note::

..     You will need the package ``dvipng`` to display the math properly in HTML.

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

.. code-block:: rest

    .. todo::

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


.. toctree::
   :hidden:

   documenting_dolfin_api
   documenting_dolfin_demos

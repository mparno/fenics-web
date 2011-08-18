.. _styleguides_sphinx_documenting_demos:

========================
Documenting DOLFIN demos
========================

The documentation for the DOLFIN demos is written by hand and located
together with the demos in the DOLFIN source tree. To document a (new)
DOLFIN demo located in the directory ``foo`` (for instance
``pde/poisson``), follow the below 2 steps. In general, the simplest
way is probably to look at one of the documented demos for instance
(``demo/pde/poisson/``) and follow the same setup.

#. Add these 3 files

    * ``foo/common.txt``

         This file should contain common information such as the main
         features the demo illustrates and, if applicable, a
         mathematical description of the differential equation that is
         solved. This file should then be included in the C++ and
         Python versions (see
         :ref:`styleguides_sphinx_demo_common_information`).

    * ``foo/cpp/documentation.rst``

         This file should contain the reST source file with the
         documentation that is specific to the C++ version of the demo
         (see :ref:`styleguides_sphinx_specific_demo_information`).

    * ``foo/python/documentation.rst``

         This file should contain the reST source file with the
         documentation that is specific to the Python version of the
         demo (see
         :ref:`styleguides_sphinx_specific_demo_information`).  .

   If either the C++ or the Python version of the demo does not exist,
   add the version and continue.

#. Move the directory ``foo`` from the directory ``undocumented`` to
   the suitable directory (for instance ``pde`` or ``la``).


.. note::

   The demo documentation is automatically included in the complete
   DOLFIN documentation when running ``make doc`` after building
   DOLFIN. While documenting a demo, it may be handy to only run
   ``make doc_demo`` and then ``make doc_html_[python|cpp]``.

.. note::

   Tests for the validity of the code snippets used in the demo
   documentation are included in the standard DOLFIN tests.

.. _styleguides_sphinx_demo_common_information:

Common information
^^^^^^^^^^^^^^^^^^

Each demo should be available in both a C++ and a Python version.
However, the summary (describing what features are demonstrated) along
with the problem and method description are typically identical for
both versions.  It is therefore desirable to put this information in a
common source file to avoid code duplication.  This common code is
placed in the file ``foo/common.txt``, which is then included in the
two files ``foo/cpp/documentation.rst`` and
``foo/python/documentation.rst`` using the ``include`` directive with
the relative path to the file:

.. code-block:: rest

  .. include:: ../common.txt

Including images
^^^^^^^^^^^^^^^^

For some demos a picture showing the solution obtained from the input
is appropriate. The picture can be saved by pressing the ``i`` key in
the Viper plot window. This file e.g., ``simulation0000.png``, should
then be copied to the same location as the ``common.txt`` file (but
use a different name for instance ``foo_u.png`` to avoid name
conflicts) and included like this:

.. code-block:: rest

    .. image:: ../foo_u.png
        :scale: 75
        :align: center

(see :ref:`styleguides_sphinx_images` for more details).
Notice that the parent directory is included in the path although the two files
``common.txt`` and ``poisson_u.png`` are located in the same directory.
This is necessary because the ``common.txt`` file will be included in files in
the two sub directories ``cpp`` and ``python``. The picture should be included
immediately after the input information.

.. _styleguides_sphinx_specific_demo_information:

C++ and Python specific contents
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The C++ and Python documentation reST source files should

#. Explain each step of the solution procedure. Do this by including
   and explaining :ref:`styleguides_sphinx_code_snippets` from the
   demo source code.

#. Include links to the API documentation using the ``:cpp:class:``
   and ``:py:class:`` directives. Note that for the Python classes,
   the full module path is required (for instance
   ``py:class:dolfin.cpp.NewtonSolver``)

#. Include the complete set of files needed to run the demo using the
   ``include`` directive.

.. todo:      MER: Add better links above

.. Additional information
.. ^^^^^^^^^^^^^^^^^^^^^^

.. Use the ``note`` and ``warning`` directives to highlight important
.. information.  The ``seealso`` directive should be used when pointing
.. to alternative solutions or functions in the
.. :ref:`programmers_reference_index`.

.. Keywords should be added to the index, using the ``index`` directive to make
.. the documentation easier to navigate through.

.. See `the Sphinx documentation
.. <http://sphinx.pocoo.org/markup/para.html#index-generating-markup>`_
.. for how to use the above directives.

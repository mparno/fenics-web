..  Style guides for FEniCS documentation (Sphinx) and DOLFIN (C++).

.. _styleguides:

************
Style guides
************

.. _styleguides_cpp_coding_style:

C++ coding style for DOLFIN
===========================

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

    vofoo();
    vobar();
    vofoo_bar(...);

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

Indentation
^^^^^^^^^^^

Indentation should be two spaces and it should be spaces, **not**
tab(s).

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
      computeEntities(mesh, d0);
    if (topology.size(d1) == 0)
      computeEntities(mesh, d1);

    // Check if connectivity still needs to be computed
    if (connectivity.size() > 0)
      return;

    ...

Integers and reals
^^^^^^^^^^^^^^^^^^

Use ``dolfin::uint`` instead of ``int`` (unless you really
want to use negative integers which is rare) and ``dolfin::real``
instead of ``double``:

.. code-block:: c++

    uint i = 0;
    double x = 0.0;

These are typedefs for the standard C++ types ``unsigned int``
and ``double`` (defined in ``dolfin/common/types.h``).

Placement of brackets
^^^^^^^^^^^^^^^^^^^^^

Curly brackets following a control statement should appear in the next
line and not be indented:

.. code-block:: c++

    for (uint i = 0; i < 10; i++)
    {
      ...
    }

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

Including header files
^^^^^^^^^^^^^^^^^^^^^^

Don't use ``#include <dolfin.h>`` or ``#include``
``<dolfin/dolfin_foo.h>`` inside the DOLFIN kernel. Only include the
portions of DOLFIN you are actually using.

Forward declarations
^^^^^^^^^^^^^^^^^^^^

Actually, try to include as little as possible and use forward
declarations whenever possible (in header files). Put the
``#include`` in the implementation file.

Explicit constructors
^^^^^^^^^^^^^^^^^^^^^

Make all constructors (except copy constructors) explicit if there is no
particular reason not to do so:

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

    class Bar
    {
      virtual void foo();
      virtual void bar();
    };

Use of libraries
----------------

Prefer C++ strings and streams to old C-style ``char*``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use ``std::string`` instead of ``const char*`` and use ``std::istream`` and
``std::ostream`` instead of ``FILE``. Avoid ``printf``,
``sprintf`` and the like.

There are exceptions to this rule where we need to use old C-style
function calls. One such exception is handling of command-line
arguments (``char* argv[]``).

.. _styleguides_sphinx_coding_style:

Sphinx coding style for FEniCS documentation
============================================

Use note for doc-authors in case of missing documentation, things to be
considered etc.

Put the style guide and some notes on how to document classes, functions,
write tutorials/examples etc.

..  Section markers from http://docs.python.org/documenting/rest.html
    # with overline, for parts
    * with overline, for chapters
    =, for sections
    -, for subsections
    ^, for subsubsections
    ", for paragraphs


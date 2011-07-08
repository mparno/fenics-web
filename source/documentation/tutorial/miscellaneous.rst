.. Automatically generated reST file from Doconce source
   (http://code.google.com/p/doconce/)


Miscellaneous Topics
====================

Glossary
--------


.. index:: self

.. index:: FEniCS

.. index:: DOLFIN

.. index:: Viper

.. index:: UFL

.. index:: class

.. index:: instance


.. index:: method (class)

.. index:: attribute (class)


Below we explain some key terms used in this tutorial.

  FEniCS: name of a software suite composed of many individual software
          components (see ``fenics.org``). Some components are DOLFIN and
	  Viper, explicitly referred to in this tutorial. Others are
          FFC and FIAT, heavily used by the programs appearing in this tutorial,
          but never explicitly used from the programs.

  DOLFIN: a FEniCS component, more precisely a C++ library, with
          a Python interface, for performing important actions in finite element
          programs. DOLFIN makes use of many other FEniCS components and
          many external software packages.

  Viper:  a FEniCS component for quick visualization of finite element
          meshes and solutions.

  UFL:    a FEniCS component implementing the *unified form language*
          for specifying finite element forms in FEniCS programs.
          The definition of the forms, typically called ``a`` and ``L`` in
          this tutorial, must have legal UFL syntax. The same applies to
          the definition of functionals (see the section :ref:`tut:poisson1:functionals`).

  Class (Python): a programming construction for creating objects
          containing a set of variables and functions. Most
          types of FEniCS objects are defined through the class concept.

  Instance (Python): an object of a particular type, where the type is
          implemented as a class. For instance,
          ``mesh = UnitInterval(10)`` creates
          an instance of class ``UnitInterval``, which is reached by the
          name ``mesh``. (Class ``UnitInterval`` is actually just
          an interface to a corresponding C++ class in the DOLFIN C++ library.)

  Class method (Python): a function in a class, reached by dot
          notation: ``instance_name.method_name``

  argument ``self`` (Python): required first parameter in class methods,
         representing a particular object of the class.
         Used in method definitions, but never in calls to a method.
         For example, if ``method(self, x)`` is the definition of
         ``method`` in a class ``Y``, ``method`` is called as
         ``y.method(x)``, where ``y`` is an instance of class ``X``.
         In a call like ``y.method(x)``, ``method`` is invoked with
         ``self=y``.

  Class attribute (Python): a variable in a class, reached by
         dot notation: ``instance_name.attribute_name``


Overview of Objects and Functions
---------------------------------

Most classes in FEniCS have an explanation of the purpose and usage
that can be seen by using the general documentation command
``pydoc`` for Python objects. You can type

.. index:: pydoc


.. code-block:: py


        pydoc dolfin.X

to look up documentation of a Python class ``X`` from the DOLFIN
library (``X`` can be ``UnitSquare``, ``Function``,
``Viper``, etc.). Below is an overview of the most important classes
and functions
in FEniCS programs, in the order they typically appear within programs.

``UnitSquare(nx, ny)``: generate mesh over the unit square
:math:`[0,1]\times [0,1]` using ``nx`` divisions in :math:`x` direction and
``ny`` divisions in :math:`y` direction. Each of the ``nx*ny`` squares
are divided into two cells of triangular shape.

``UnitInterval``, ``UnitCube``, ``UnitCircle``, ``UnitSphere``,
``Interval``, ``Rectangle``, and ``Box``: generate mesh over
domains of simple geometric shape, see the section :ref:`tut:prepro`.

``FunctionSpace(mesh, element_type, degree)``:
a function space defined over a mesh, with a given element type
(e.g., ``'CG'`` or ``'DG'``), with basis functions as polynomials of
a specified degree.

``Expression(formula)``: a scalar- or vector-valued function, given as a
mathematical expression ``formula`` (string) written in C++ syntax.

``Function(V)``: a scalar- or vector-valued finite element field in
the function space ``V``. If ``V`` is a ``FunctionSpace`` object,
``Function(V)`` becomes a scalar field, and with ``V`` as a
``VectorFunctionSpace`` object, ``Function(V)`` becomes a
vector field.

``SubDomain``: class for defining a subdomain, either a part of the
boundary, an internal boundary, or a part of the domain.
The programmer must subclass ``SubDomain`` and implement the
``inside(self, x, on_boundary)`` function
(see the section :ref:`tut:poisson1:impl`) for telling whether a
point ``x`` is inside the subdomain or not.

``Mesh``: class for representing a finite element mesh, consisting of
cells, vertices, and optionally faces, edges, and facets.

``MeshFunction``: tool for marking parts of the domain or the boundary.
Used for variable coefficients ("material properties", see
the section :ref:`tut:possion:2D:2mat:problem`) or for
boundary conditions (see the section :ref:`tut:poisson:mat:neumann`).

``DirichletBC(V, value, where)``: specification of Dirichlet (essential)
boundary conditions via a function space ``V``, a function
``value(x)`` for computing the value of the condition at a point ``x``,
and a specification ``where`` of the boundary, either as a
``SubDomain`` subclass instance, a plain function, or as a
``MeshFunction`` instance.
In the latter case, a 4th argument is provided to describe which subdomain
number that describes the relevant boundary.

``TestFunction(V)``: define a test function on a space ``V`` to be used
in a variational form.

``TrialFunction(V)``: define a trial function on a space ``V`` to be used
in a variational form to represent the unknown in a finite element problem.

``assemble(X)``: assemble a matrix, a right-hand side, or a functional,
given a from ``X`` written with UFL syntax.

``assemble_system(a, L, bc)``: assemble the matrix and the right-hand
side from a bilinear (``a``) and linear (``L``) form written with UFL
syntax. The ``bc`` parameter holds one or more ``DirichletBC`` objects.

``VariationalProblem(a, L, bc)``: define and solve a variational problem,
given a bilinear (``a``) and linear (``L``) form, written with UFL
syntax, and one or more ``DirichletBC`` objects stored in ``bc``.
A 4th argument, ``nonlinear=True``, can be given to define and solve
nonlinear variational problems (see the section :ref:`tut:nonlinear:Newton:auto`).

``solve(A, U, b)``: solve a linear system with ``A`` as coefficient
matrix (``Matrix`` object), ``U`` as unknown (``Vector`` object),
and ``b`` as right-hand side (``Vector`` object).
Usually, ``U`` is replaced by ``u.vector()``, where
``u`` is a ``Function`` object representing the unknown finite
element function of the problem, while
``A`` and ``b`` are computed by calls to ``assemble``
or ``assemble_system``.

``plot(q)``: quick visualization of a mesh, function, or mesh function
``q``, using the Viper component in FEniCS.

``interpolate(func, V)``: interpolate a formula or finite
element function ``func`` onto the function space ``V``.

``project(func, V)``: project a formula or finite element function ``func``
onto the function space ``V``.


.. _tut:app:install:

Installing FEniCS
-----------------

The FEniCS software components are available for Linux, Windows and Mac OS
X platforms. Detailed information on how to get FEniCS running on such
machines are available at the ``fenics.org`` website.
Here are just some quick descriptions and recommendations by the author.

To make the installation of FEniCS as painless and reliable as
possible, the reader is strongly recommended to use Ubuntu Linux.
(Even though Mac users now can get FEniCS by a one-click install, I
recommend using Ubuntu on Mac, unless you have high Unix competence
and much experience with compiling and linking C++ libraries on Mac OS
X.)  Any standard PC can easily be equipped with Ubuntu Linux, which
may live side by side with either Windows or Mac OS X or another Linux
installation.  Basically, you download Ubuntu from
``www.ubuntu.com/getubuntu/download``, burn the file on a CD, reboot the
machine with the CD, and answer some usually straightforward questions
(if necessary). The graphical user interface (GUI) of Ubuntu is quite
similar to both Windows 7 and Mac OS X, but to be efficient when doing
science with FEniCS this author recommends to run programs in a
terminal window and write them in a text editor like Emacs or Vim. You
can employ integrated development environment such as Eclipse, but
intensive FEniCS developers and users tend to find terminal windows
and plain text editors more user friendly.

Instead of making it possible to boot your machine with the Linux
Ubuntu operating system, you can run Ubuntu in a separate window in
your existing operation system. On Mac, you can use the VirtualBox
software available from ``http://www.virtualbox.org`` to run Ubuntu, or
you can buy a commercial tool like *VMWare Fusion* or *Parallels*.  On
Windows, Wubi makes a tool that automatically installs Ubuntu on the
machine. Just give a username and password for the Ubuntu
installation, and Wubi performs the rest. You can also use VirtualBox
on Windows machines.

Once the Ubuntu window
is up and running, FEniCS is painlessly installed by

.. code-block:: console

        sudo apt-get install fenics

Sometimes the FEniCS software in a standard Ubuntu installation lacks
some recent features and bug fixes. Visiting ``fenics.org`` and copying just
five Unix commands is all you have to do to install a newer version of
the software.


.. _tut:appendix:books:

Books on the Finite Element Method
----------------------------------

There are a large number of books on the finite element method.  The
books typically fall in either of two categories: the abstract
mathematical version of the method and the engineering "structural
analysis" formulation. FEniCS builds heavily on concepts in the
abstract mathematical exposition.  An easy-to-read book, which
provides a good general background for using FEniCS, is Gockenbach
[Gockenbach2006]. The book by Donea and Huerta
[DoneaHuerta2003] has a similar style, but aims at readers with
interest in fluid flow problems. Hughes [Hughes1987] is also
highly recommended, especially for those interested in solid mechanics
and heat transfer applications.

Readers with background in the engineering "structural analysis"
version of the finite element method may find Bickford
[Bickford1994] as an attractive bridge over to the abstract
mathematical formulation that FEniCS builds upon.  Those who have a
weak background in differential equations in general should consult a
more fundamental book, and Eriksson {\em et
al}. [ErikssonEstepHansboEtAl1996] is a very good choice.  On the
other hand, FEniCS users with a strong background in mathematics and
interest in the mathematical properties of the finite element method,
will appreciate the texts by Brenner and Scott [BrennerScott2008],
Braess [Braess2007], Ern and Guermond [ErnGuermond2004],
Quarteroni and Valli [QuarteroniValli1994], or Ciarlet [Ciarlet2002].


.. _tut:appendix:pybooks:

Books on Python
---------------

Two very popular introductory books on Python are "Learning Python" by
Lutz [Lutz2007] and "Practical Python" by Hetland [Hetland2002].  More
advanced and comprehensive books include "Programming Python" by Lutz
[Lutz2006], and "Python Cookbook" [MartelliAscher2005] and "Python in
a Nutshell" [Martelli2006] by Martelli.  The web page
``http://wiki.python.org/moin/PythonBooks`` lists numerous additional
books.  Very few texts teach Python in a mathematical and numerical
context, but the references [Langtangen2008], [Langtangen2009a]
[Kiusalaas2005] are exceptions.



.. _tut:app:cpp:functions:

User-Defined Functions
----------------------


When defining a function in terms of a mathematical expression inside
a string formula, e.g.,

.. code-block:: python

        myfunc = Expression('sin(x[0])*cos(x[1])')

the expression contained in the first argument
will be turned into a C++ function
and compiled to gain efficiency. Therefore,
the syntax used in the expression must be valid C++ syntax.
Most Python syntax for mathematical expressions are also valid C++ syntax,
but power expressions make an exception: ``p**a`` must be written as
``pow(p,a)`` in C++ (this is also an alternative Python syntax).
The following mathematical functions can be used directly
in C++
expressions when defining ``Expression`` objects:
``cos``, ``sin``, ``tan``, ``acos``, ``asin``,
``atan``, ``atan2``, ``cosh``, ``sinh``, ``tanh``, ``exp``,
``frexp``, ``ldexp``, ``log``, ``log10``, ``modf``,
``pow``, ``sqrt``, ``ceil``, ``fabs``, ``floor``, and ``fmod``.
Moreover, the number :math:`\pi` is available as the symbol ``pi``.
All the listed functions are taken from the ``cmath`` C++ header file, and
one may hence
consult documentation of ``cmath`` for more information on the
various functions.

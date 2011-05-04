.. UFL user documentation

.. _ufl_user_documentation:

######################
UFL user documentation
######################

*****************
About this manual
*****************

Intended audience
=================

This manual is written both for the beginning and the advanced user.
There is also some useful information for developers. More advanced
topics are treated at the end of the manual or in the appendix.

Typographic conventions
=======================

Code is written in monospace ``like this``. Commands that should be
entered in a Unix shell::

  # ./configure
  # make

Commands are written in the dialect of the ``bash`` shell. For other
shells, such as ``tcsh`, appropriate translations may be needed.



Enumeration and list indices
============================

Throughout this manual, elements :math:`x_i` of sets :math:`\{x_i\}`
of size :math:`n` are enumerated from :math:`i = 0` to :math:`i =
n-1`. Derivatives in :math:`\mathbb{R}^n` are enumerated similarly
:math:`\frac{\partial}{\partial x_0}, \frac{\partial}{\partial x_1}, \ldots, \frac{\partial}{\partial x_{n-1}}`.

Contact
=======

Comments, corrections and contributions to this manual are most welcome
and should be sent to ufl@lists.launchpad.net


************
Introduction
************

The Unified Form Language (UFL) is a domain specific language for
defining discrete variational forms and functionals in a notation close
to pen-and-paper formulation.

UFL is part of the FEniCS Project, and is usually used in combination
with other components from this project to compute solutions to partial
differential equations. The form compilers FFC and SFC use UFL as their
end-user interface, producing implementations of the UFC interface as
their output. See the DOLFIN manual for more details about using UFL in
an integrated problem solving environment.

This manual is intended for different audiences.  If you are an end user
and all you want to do is to solve your PDEs with the FEniCS framework,
Chapters XX and XX are for
you. These two chapters explain how to use all operators available in
the language and present a number of examples to illustrate the use of
the form language in applications. The rest of the chapters contain more
technical details intended for developers who need to understand what
is happening behind the scenes and modify or extend UFL in the future.

Chapter XX details the implementation of the language, in particular
how expressions are represented internally by UFL.  This can also be
useful knowledge to understand error messages and debug errors in your
form files.

Chapter XX explains many algorithms to work with UFL expressions,
mostly intended to aid developers of form compilers.  The algorithms
available includes helper functions for easy and efficient iteration
over expression trees, formatting tools to present expressions as text or
images of different kinds, utilities to analyse properties of expressions
or checking their validity, automatic differentiation algorithms, as
well as algorithms to work with the computational graphs of expressions.

*************
Form language
*************

UFL consists of a set of operators and atomic expressions that can be
used to express variational forms and functionals.  Below we will define
all these operators and atomic expressions in detail.

UFL is built on top of, or embedded in, the high level language Python.
Since the form language is built on top of Python, any Python code is
valid in the definition of a form (but not all Python code defines a
multilinear form).  In particular, comments (lines starting with ``#``
and functions (keyword \texttt{def}, see Section \ref{sec:user-defined}
below) are useful in the definition of a form.  However, it is usually a
good idea to avoid using advanced Python features in the form definition,
to stay close to the mathematical notation.

The entire form language can be imported in Python with the line

.. code-block:: python

   from ufl import *

which is assumed in all examples below and can be omitted in ``.ufl``
files.  This can be useful for experimenting with the language in an
interactive Python interpreter.


Forms and integrals
===================

UFL is designed to express forms in the following generalized format:

.. math::

   a(v_1, \ldots, v_r; w_1, \ldots,  w_n)
      =
          \sum_{k=1}^{n_c} \int_{\Omega_k}
                I^c_k(v_1, \ldots, v_r; w_1, \ldots w_n) dx
         +     \sum_{k=1}^{n_e} \int_{\partial\Omega_k}
                I^e_k(v_1, \ldots, v_r; w_1, \ldots,  w_n) ds
         +     \sum_{k=1}^{n_i} \int_{\Gamma_k}
                I^i_k(v_1, \ldots, v_r; w_1, \ldots, w_n) dS.

Here the form :math:`a` depends on the *form arguments* :math:`v_1,
\ldots, v_r` and the *form coefficients* :math:`w_1, \ldots, w_n`,
and its expression is a sum of integrals.  Each term of a valid form
expression must be a scalar-valued expression integrated exactly once. How
to define form arguments and integrand expressions is detailed in the
rest of this chapter.

Integrals are expressed through multiplication with a measure,
representing an integral over either of

* the interior of the domain :math:`\Omega` (:math:`dx`, cell integral);

* the boundary :math:`\partial\Omega` of :math:`\Omega` (:math:`ds`,
  exterior facet integral);

* the set of interior facets :math:`\Gamma` (:math:`dS`, interior facet
  integral).

UFL declares the measures ``dx`` :math:`\leftrightarrow dx`, ``ds``
:math:`\leftrightarrow ds`, and ``dS`` :math:`\leftrightarrow dS`.

As a basic example, assume ``v`` is a scalar-valued expression and
consider the integral of ``v`` over the interior of :math:`\Omega`. This
may be expressed as::

  a = v*dx

and the integral of ``v`` over :math:`\partial\Omega` is written as::

  a = v*ds

Alternatively, measures can be redefined to represent numbered subsets of
a domain, such that a form can take on different expressions on different
parts of the domain.  If ``c``, ``e0`` and ``e1`` are scalar-valued
expressions, then::

  a = c*dx + e0*ds(0) + e1*ds(1)

represents

.. math::

   a = \int_\Omega c\,dx + \int_{\partial\Omega_0} e_0 \, ds + \int_{\partial\Omega_1} e_1 \, ds.

where

.. math::

   \partial\Omega_0 \subset \partial\Omega, \qquad \partial\Omega_1 \subset \partial\Omega.

Generalizing this further we end up with the expression \eqref{eq:form_integrals}.
Note that the domain :math:`\Omega` and its subdomains and boundaries
are not known to UFL. These will not enter the stage until
you start using UFL in a problem solving environment like DOLFIN.

.. topic:: Advanced usage

  A feature for advanced users is attaching metadata to integrals.
  This can be used to define different quadrature degrees for different
  terms in a form, and to override other form compiler specific options
  separately for different terms::

    a = c0*dx(0, metadata0) + c1*dx(1, metadata1)

  The convention is that metadata should be a dict with any of the
  following keys:

  * ``"integration_order"``: Integer defining the polynomial order
    that should be integrated exactly. This is a compilation hint, and the
    form compiler is free to ignore this if for example exact integration
    is being used.

  * ``"ffc"``: A dict with further FFC specific options, see the
    FFC manual.

  * ``"sfc"``: A dict with further SFC specific options, see the
    SFC manual.

  * Other string: A dict with further options specific to some other
    external code.

  Other standardized options may be added in later versions. ::

    metadata0 = {"ffc": {"representation": "quadrature"}}
    metadata1 = {"integration_order": 7,
             "ffc": {"representation": "tensor"}}

    a = v*u*dx(0, metadata1) + f*v*dx(0, metadata2)


Finite element spaces
=====================

Before we can explain how form arguments are declared, we need to show how
to define function spaces.  UFL can represent very flexible general
hierarchies of mixed finite elements, and has predefined names for most
common element families.


Cells
-----

A polygonal cell is defined by a basic shape and a degree (note
that the other components of FEniCS do not yet handle cells of higher
degree, so this will only be useful in the future), written like::

  cell = Cell(shape, degree)

Valid shapes are "interval", "triangle", "tetrahedron", "quadrilateral",
and "hexahedron".  Some examples::

  # Cubic triangle cell
  cell = Cell("triangle", 3)

  # Quadratic tetrahedron cell
  cell = Cell("tetrahedron", 2)

Objects for linear cells of all basic shapes are predefined::

  # Predefined linear cells
  cell = interval
  cell = triangle
  cell = tetrahedron
  cell = quadrilateral
  cell = hexahedron

In the rest of this document, a variable name ``cell`` will be used where
any cell is a valid argument, to make the examples dimension independent
wherever possible.  Using a variable ``cell`` to hold the cell type used
in a form is highly recommended, since this makes most form definitions
dimension independent.


Element families
----------------

UFL predefines a set of names of known element families.  When defining
a finite element below, the argument ``family`` is a string and its
possible values include:

* ``"Lagrange"`` or ``"CG"``, representing standard scalar
  Lagrange finite elements (continuous piecewise polynomial functions);

* ``"Discontinuous Lagrange"`` or ``"DG"``, representing
  scalar discontinuous Lagrange finite elements (discontinuous piecewise
  polynomial functions);

* ``"Crouzeix-Raviart"`` or ``"CR"``, representing scalar
  Crouzeix--Raviart elements;

* ``"Brezzi-Douglas-Marini"`` or ``"BDM"``, representing
  vector-valued Brezzi--Douglas--Marini :math:`H(\mathrm{div})` elements;

* ``"Brezzi-Douglas-Fortin-Marini`` or ``"BDFM"``, representing
  vector-valued Brezzi--Douglas--Fortin--Marini :math:`H(\mathrm{div})`
  elements;

* ``"Raviart-Thomas"`` or ``"RT"``, representing
  vector-valued Raviart--Thomas :math:`H(\mathrm{div})` elements.

* ``"Nedelec 1st kind H(div)"`` or ``"N1div"``,
  representing vector-valued Nedelec :math:`H(\mathrm{div})` elements
  (of the first kind).

* ``"Nedelec 2st kind H(div)"`` or ``"N2div"``,
  representing vector-valued Nedelec :math:`H(\mathrm{div})` elements
  (of the second kind).

* ``"Nedelec 1st kind H(curl)"`` or ``"N1curl"``, representing
  vector-valued Nedelec :math:`H(\mathrm{curl})` elements
  (of the first kind).

* ``"Nedelec 2st kind H(curl)"`` or ``"N2curl"``,
  representing vector-valued Nedelec :math:`H(\mathrm{curl})` elements
  (of the second kind).

  %\item
  %  \texttt{"Bubble"} or \texttt{"B"}, representing FIXME;

* ``"Quadrature"`` or ``"Q"``, representing artificial ``finite elements``
  with degrees of freedom being function evaluation at quadrature points;

* ``"Boundary Quadrature"`` or ``"BQ"``, representing artificial
  ``finite elements'' with degrees of freedom being function evaluation
  at quadrature points on the boundary;


.. topic:: Advanced usage

  New elements can be added dynamically by the form compiler using the
  function ``register_element``. See the docstring for details.
  To see which elements are registered (including the standard built in
  ones listed above) call the function ``show_elements``.


Basic elements
--------------

A ``FiniteElement``, sometimes called a basic element, represents a
finite element in some family on a given cell with a certain polynomial
degree. Valid families and cells are explained above.
The notation is::

  element = FiniteElement(family, cell, degree)

Some examples::

  element = FiniteElement("Lagrange", interval, 3)
  element = FiniteElement("DG", tetrahedron, 0)
  element = FiniteElement("BDM", triangle, 1)


Vector elements
---------------

A ``VectorElement`` represents a combination of basic elements such that
each component of a vector is represented by the basic element. The size
is usually omitted, the default size equals the geometry dimension.
The notation is::

  element = VectorElement(family, cell, degree[, size])

Some examples::

  element = VectorElement("CG", triangle, 2)
  element = VectorElement("DG", tetrahedron, 0, size=6)


Tensor elements
---------------

A ``TensorElement`` represents a combination of basic elements such that
each component of a tensor is represented by the basic element. The
shape is usually omitted, the default shape is (d, d) where d is the
geometry dimension. The notation is::

  element = TensorElement(family, cell, degree[, shape, symmetry])

Any shape tuple consisting of positive integers is valid,
and the optional symmetry can either be set to ``True``
which means standard matrix symmetry (like :math:`A_{ij} = A_{ji}`),
or a ``dict`` like ``{ (0,1):(1,0), (0,2):(2,0) }``
where the ``dict`` keys are index tuples that are
represented by the corresponding ``dict`` value.

Examples::

  element = TensorElement("CG", cell, 2)
  element = TensorElement("DG", cell, 0, shape=(6,6))
  element = TensorElement("DG", cell, 0, symmetry=True)
  element = TensorElement("DG", cell, 0, symmetry={(0,0): (1,1)})


Mixed elements
--------------

A ``MixedElement` represents an arbitrary combination of other elements.
``VectorElement`` and ``TensorElement`` are special cases of a
``MixedElement`` where all sub-elements are equal.

General notation for an arbitrary number of subelements::

  element = MixedElement(element1, element2[, element3, ...])

Shorthand notation for two subelements::

  element = element1 * element2

Note: Multiplication is a binary operator, such that ::

  element = element1 * element2 * element3

represents ``(e1 * e2) * e3}, i.e. this is a mixed element with two
sub-elements ``(e1 * e2)`` and ``e3``.

See section~\ref{sec:formarguments} for details on how defining
functions on mixed spaces can differ from functions on other
finite element spaces.

Examples::

  # Taylor-Hood element
  V = VectorElement("Lagrange", cell, 2)
  P = FiniteElement("Lagrange", cell, 1)
  TH = V * P

  # A tensor-vector-scalar element
  T = TensorElement("Lagrange", cell, 2, symmetry=True)
  V = VectorElement("Lagrange", cell, 1)
  P = FiniteElement("DG", cell, 0)
  ME = MixedElement(T, V, P)

EnrichedElement
---------------

The data type ``EnrichedElement`` represents the vector sum of two
(or more) finite elements.

Example: The Mini element can be constructed as::

  P1 = VectorElement("Lagrange", "triangle", 1)
  B  = VectorElement("Bubble", "triangle", 3)
  Q  = FiniteElement("Lagrange", "triangle", 1)

  Mini = (P1 + B) * Q

Form arguments
==============

Form arguments are divided in two groups, basis functions and
functions (the term *function* in UFL maps to the term
*coefficient* in UFC).  A ``BasisFunction`` represents an
arbitrary basis function in a given discrete finite element space,
while a ``Function`` represents a function in a discrete finite
element space that will be provided by the user at a later stage. The
number of ``BasisFunction``\ s that occur in a ``Form`` equals
the arity of the form.

Basis functions
---------------

The data type ``BasisFunction`` represents a basis function on a
given finite element. A ``BasisFunction`` must be created for a
previously declared finite element (simple or mixed)::

  v = BasisFunction(element)

Note that more than one ``BasisFunction`` can be declared for the same
``FiniteElement``. Basis functions are associated with the arguments of
a multilinear form in the order of declaration.

For a ``MixedElement``, the function ``BasisFunctions`` can be used to
construct tuples of ``BasisFunction``\ s, as illustrated here for a mixed
Taylor--Hood element::

  v, q = BasisFunctions(TH)
  u, p = BasisFunctions(TH)

For a ``BasisFunction`` on a ``MixedElement`` (or ``VectorElement``
or ``TensorElement``), the function ``split`` can be used to extract
basis function values on subspaces, as illustrated here for a mixed
Taylor--Hood element::

  vq = BasisFunction(TH)
  v, q = split(up)

A shorthand for this is in place called ``BasisFunctions``::

  v, q = BasisFunctions(TH)

For convenience, ``TestFunction`` and ``TrialFunction`` are special
instances of ``BasisFunction`` with the property that a ``TestFunction``
will always be the first argument in a form and ``TrialFunction`` will
always be the second argument in a form (order of declaration does
not matter).  Their usage is otherwise the same as for ``BasisFunction``::

  v = TestFunction(element)
  u = TrialFunction(element)
  v, q = TestFunctions(TH)
  u, p = TrialFunctions(TH)


Coefficient functions
---------------------

The data type ``Function`` represents a function belonging to a given
finite element space, that is, a linear combination of basis functions
of the finite element space. A ``Function`` must be declared for a
previously declared ``FiniteElement``::

  f = Function(element)

Note that the order in which ``Function``\ s are declared is important,
directly reflected in the ordering they have among the arguments to each
``Form`` they are part of.

``Function`` is used to represent user-defined functions, including, e.g.,
source terms, body forces, variable coefficients and stabilization terms.
UFL treats each ``Function`` as a linear combination of unknown basis
functions with unknown coefficients, that is, UFL knows nothing about
the concrete basis functions of the element and nothing about the value
of the function.

Note that more than one function can be declared for the same
``FiniteElement``. The following example declares two ``BasisFunction``_s
and two ``Function``\ s for the same ``FiniteElement``::

  v = BasisFunction(element)
  u = BasisFunction(element)
  f = Function(element)
  g = Function(element)

For a ``Function` on a ``MixedElement`` (or ``VectorElement`` or
``TensorElement``), the function ``split`` can be used to extract function
values on subspaces, as illustrated here for a mixed Taylor--Hood element::

  up = Function(TH)
  u, p = split(up)

A shorthand for this is in place called ``Functions``::

  u, p = Function(TH)

Spatially constant (or discontinuous piecewise constant) functions can
conveniently be represented by ``Constant``, ``VectorConstant``, and
``TensorConstant``::

  c0 = Constant(cell)
  v0 = VectorConstant(cell)
  t0 = TensorConstant(cell)

These three lines are equivalent with first defining
DG0 elements and then defining a ``Function``
on each, illustrated here::

  DG0 = FiniteElement("Discontinuous Lagrange", cell, 0)
  DG0v = VectorElement("Discontinuous Lagrange", cell, 0)
  DG0t = TensorElement("Discontinuous Lagrange", cell, 0)

  c1 = Function(DG0)
  v1 = Function(DG0v)
  t1 = Function(DG0t)

Basic Datatypes
===============

UFL expressions can depend on some other quantities in addition to the
functions and basis functions described above.

Literals and geometric quantities
---------------------------------

Some atomic quantities are derived from the cell.  For example, the
(global) spatial coordinates are available as a vector valued expression
``cell.x``::

  # Linear form for a load vector with a sin(y) coefficient
  v = TestFunction(element)
  x = cell.x
  L = sin(x[1])*v*dx

Another quantity is the (outwards pointing) facet normal ``cell.n``.
The normal vector is only defined on the boundary, so it can't be used
in a cell integral.

Example functional ``M``, an integral of the normal component of a
function ``g`` over the boundary::

  n = cell.n
  g = Function(VectorElement("CG", cell, 1))
  M = dot(n, g)*ds

Python scalars (int, float) can be used anywhere a scalar expression
is allowed. Another literal constant type is ``Identity`` which
represents an :math:`n\times n` unit matrix of given size :math:`n`, as in this example::

  # Geometric dimension
  d = cell.d

  # d x d identiy matrix
  I = Identity(d)

  # Kronecker delta
  delta_ij = I[i,j]

.. note: Advanced usage

  Note that there are some differences from FFC.
  In particular, using ``FacetNormal`` or ``cell.n``
  does not implicitly add another coefficient Function to the form,
  the normal should be automatically computed in UFC code.
  Note also that ``MeshSize`` has been removed because the
  meaning is ambiguous (does it mean min, max, avg, cell radius?),
  so use a ``Constant`` instead.


Indexing and tensor components
==============================

UFL supports index notation, which is often a convenient way to
express forms. The basic principle of index notation is that summation
is implicit over indices repeated twice in each term of an expression.
The following examples illustrate the index notation, assuming that
each of the variables ``i`` and ``j`` have been declared as
a free ``Index``:

.. math::

   \mbox{``v[i]*w[i]``} &\leftrightarrow& \sum_{i=0}^{n-1} v_i w_i = \mathbf{v}\cdot\mathbf{w}, \\
   \mbox{``Dx(v, i)*Dx(w, i)``} &\leftrightarrow&
   \sum_{i=0}^{d-1}
   \frac{\partial v}{\partial x_i}
   \frac{\partial w}{\partial x_i} = \nabla v \cdot \nabla w, \\
   \mbox{``Dx(v[i], i)``} &\leftrightarrow& \sum_{i=0}^{d-1}
   \frac{\partial v_i}{\partial x_i} = \nabla \cdot v, \\
   \mbox{``Dx(v[i], j)*Dx(w[i], j)``} &\leftrightarrow&
   \sum_{i=0}^{n-1} \sum_{j=0}^{d-1}
   \frac{\partial v_i}{\partial x_j} \frac{\partial w_i}{\partial x_j} = \nabla \mathbf{v} : \nabla \mathbf{w}.

Here we'll try to very briefly summarize the basic concepts of tensor
algebra and index notation, just enough to express the operators in UFL.

Assuming an Euclidean space in :math:`d` dimensions with :math:`1 \le
d 3`, and a set of orthonormal basis vectors :math:`\mathbf{i}_i` for :math:`i
\in {0, \ldots, d-1 }`, we can define the dot product of any two basis
functions as

.. math::

   \mathbf{i}_{i} \cdot \mathbf{i}_{j} = \delta_{ij},

where :math:`\delta_{ij}` is the Kronecker delta

.. math::

   \delta_{ij}
   \equiv
   \begin{cases}
   1, \quad i = j, \\
   0, \quad \text{otherwise}.
   \end{cases}

A rank 1 tensor (vector) quantity :math:`\mathbf{v}` can be represented in
terms of unit vectors and its scalar components in that basis.  In tensor
algebra it is common to assume implicit summation over indices repeated
twice in a product::

.. math::

   \mathbf{v} = v_k \mathbf{i}_k \equiv \sum_k v_k \mathbf{i}_k.

Similarly, a rank two tensor (matrix) quantity :math:`\mathbf{A}` can
be represented in terms of unit matrices, that is outer products of
unit vectors:

.. math::

   \mathbf{A} = A_{ij} \mathbf{i}_i \mathbf{i}_j \equiv \sum_i \sum_j A_{ij} \mathbf{i}_i \mathbf{i}_j .

This generalizes to tensors of arbitrary rank:

.. math::

   \mathcal{C} &= C_\iota \mathbf{i}_{\iota_0} \otimes \cdots \otimes \mathbf{i}_{\iota_{r-1}} \\
    &\equiv \sum_{\iota_0} \cdots \sum_{\iota_{r-1}}
    C_\iota \mathbf{i}_{\iota_0}\otimes\cdots \otimes \mathbf{i}_{\iota_{r-1}},

where :math:`\mathcal{C}` is a rank :math:`r` tensor and :math:`\iota`
is a multi-index of length :math:`r`.

%TODO: More about tensor algebra concepts to better support
%      the following sections? I don't know how much we can
%      assume people knows about this?

When writing equations on paper, a mathematician can easily switch
between the :math:`\mathbf{v}` and :math:`v_i` representations without
stating it explicitly. This is possible because of flexible notation
and conventions. In a programming language, we can't use the boldface
notation which associates :math:`\mathbf{v}` and :math:`v` by convention,
and we can't always interpret such conventions unambiguously.  Therefore,
UFL requires that an expression is explicitly mapped from its tensor
representation (:math:`\mathbf{v}`, :math:`\mathbf{A}`) to its component
representation (:math:`v_i`, :math:`A_{ij}`) and back.  This is done using
``Index`` objects, the indexing operator (``v[i]``), and the function
``as_tensor``.  More details on these follow.

In the following descriptions of UFL operator syntax, i-l and p-s are
assumed to be predefined indices, and unless otherwise specified the name
v refers to some vector valued expression, and the name A refers to some
matrix valued expression.  The name C refers to a tensor expression of
arbitrary rank.

Defining indices
----------------

A set of indices ``i``, ``j``, ``k``, ``l`` and ``p``, ``q``, ``r``,
``s`` are predefined, and these should be enough for many applications.
Examples will usually use these objects instead of creating new ones to
conserve space.

The data type ``Index`` represents an index used for subscripting
derivatives or taking components of non-scalar expressions.
To create indices, you can either make a single using ``Index()``
or make several at once conveniently using ``indices(n)``::

  i = Index()
  j, k, l = indices(3)

Each of these represents an ``index range`` determined by the context;
if used to subscript a tensor-valued expression, the range is given
by the shape of the expression, and if used to subscript a derivative,
the range is given by the dimension :math:`d` of the underlying shape
of the finite element space.  As we shall see below, indices can be a
powerful tool when used to define forms in tensor notation.


.. note: Advanced usage

  If using UFL inside PyDOLFIN or another larger programming environment,
  it is a good idea to define your indices explicitly just before your
  form uses them, to avoid name collisions.  The definition of the
  predefined indices is simply::

    i, j, k, l = indices(4)
    p, q, r, s = indices(4)

.. note: Advanced usage

  Note that in the old FFC notation, the definition ::

    i = Index(0)

  meant that the value of the index remained constant.  This does not mean
  the same in UFL, and this notation is only meant for internal usage.
  Fixed indices are simply integers instead::

    i = 0


Taking components of tensors
----------------------------
% TODO: Explain in more words

Basic fixed indexing of a vector valued expression v or matrix valued
expression A:

* ``v[0]``: component access, representing the scalar value of the first
  component of v

* ``A[0,1]``: component access, representing the scalar value of the
  first row, second column of A


Basic indexing:
* ``v[i]``: component access, representing the scalar value of some
  component of v
* ``A[i,j]``: component access, representing the scalar value of some
  component i,j of A

More advanced indexing:

* ``A[i,0]``: component access, representing the scalar value of some
  component i of the first column of A

* ``A[i,:]``: row access, representing some row i of A, i.e. rank(A[i,:]) == 1

* ``A[:,j]``: column access, representing some column j of A,
  i.e. rank(A[:,j]) == 1

* ``C[...,0]``: subtensor access, representing the subtensor of A with
  the last axis fixed, e.g., A[...,0] == A[:,0]

* ``C[j,...]``: subtensor access, representing the subtensor of A with
  the last axis fixed, e.g., A[j,...] == A[j,:]


Making tensors from components
------------------------------

If you have expressions for scalar components of a tensor and wish to
convert them to a tensor, there are two ways to do it. If you have a
single expression with free indices that should map to tensor axes,
like mapping :math:`v_k` to :math:`\mathbf{v}` or :math:`A_{ij}` to
:math:`\mathbf{A}`, the following examples show how this is done::

  vk = Identity(cell.d)[0,k]
  v = as_tensor(vk, (k,))

  Aij = v[i]*u[j]
  A = as_tensor(Aij, (i,j))

Here ``v`` will represent unit vector :math:`\mathbf{i}_0`, and ``A``
will represent the outer product of ``v`` and ``u``.

If you have multiple expressions without indices, you can build tensors
from them just as easily, as illustrated here::

  v = as_vector([1.0, 2.0, 3.0])
  A = as_matrix([[u[0], 0], [0, u[1]]])
  B = as_matrix([[a+b for b in range(2)] for a in range(2)])

Here ``v``, ``A`` and ``B`` will represent the expressions

.. math::

   \mathbf{v} &= \mathbf{i}_0 + 2 \mathbf{i}_1 + 3 \mathbf{i}_2, \\
   \mathbf{A} &= \begin{bmatrix} u_0 & 0 \\ 0 & u_1 \end{bmatrix}, \\
   \mathbf{B} &= \begin{bmatrix} 0 & 1 \\ 1 & 2 \end{bmatrix}.

Note that the function ``as_tensor`` generalizes from vectors to tensors
of arbitrary rank, while the alternative functions ``as_vector`` and
``as_matrix`` work the same way but are only for constructing vectors
and matrices.  They are included for readability and convenience.


Implicit summation
------------------

Implicit summation can occur in only a few situations.  A product of
two terms that shares the same free index is implicitly treated as a
sum over that free index:

* ``v[i]*v[i]``: :math:`\sum_i v_i v_i`
* ``A[i,j]*v[i]*v[j]``: :math:`\sum_j (\sum_i A_{ij} v_i) v_j`

A tensor valued expression indexed twice with the same free index is
treated as a sum over that free index:

* ``A[i,i]``: :math:`\sum_i A_{ii}`
* ``C[i,j,j,i]``: :math:`\sum_i \sum_j C_{ijji}`

The spatial derivative, in the direction of a free index, of an expression
with the same free index, is treated as a sum over that free index:

* ``v[i].dx(i)``: :math:`\sum_i v_i`
* ``A[i,j].dx(i)``: :math:`\sum_i \frac{d(A_{ij})}{dx_i}`

Note that these examples are some times written :math:`v_{i,i}` and
:math:`A_{ij,i}` in pen-and-paper index notation.


Basic algebraic operators
=========================

The basic algebraic operators ``+``, ``-``, ``*``, ``/`` can be used
freely on UFL expressions. They do have some requirements on their
operands, summarized here:

Addition or subtraction, ``a + b`` or ``a - b``:

* The operands a and b must have the same shape.
* The operands a and b must have the same set of free indices.

Division, ``a / b``:

* The operand b must be a scalar expression.

* The operand b must have no free indices.

* The operand a can be non-scalar with free indices, in which division
  represents scalar division of all components with the scalar b.

Multiplication, ``a * b``:

* The only non-scalar operations allowed is scalar-tensor,
  matrix-vector and matrix-matrix multiplication.

* If either of the operands have any free indices, both must be scalar.

* If any free indices are repeated, summation is implied.


Basic nonlinear functions
=========================

Some basic nonlinear functions are also available, their meaning mostly
obvious.

* ``abs(f)``: the absolute value of f.

* ``sign(f)``: the sign of f (+1 or -1).

* ``pow(f, g)`` or ``f**g``

* ``sqrt(f)``

* ``exp(f)``

* ``ln(f)``

* ``cos(f)``

* ``sin(f)``

These functions do not accept non-scalar operands or operands with free
indices or ``BasisFunction`` dependencies.


Tensor algebra operators
========================

``transpose``
-------------

The transpose of a matrix A can be written as::

  AT = transpose(A)
  AT = A.T
  AT = as_matrix(A[i,j], (j,i))

The definition of the transpose is
\begin{align}
  \mbox{\texttt{AT[i,j]}} \leftrightarrow (\AA^{\top})_{ij} = \AA_{ji}.
\end{align}

For transposing higher order tensor expressions, index notation can
be used::

  AT = as_tensor(A[i,j,k,l], (l,k,j,i))

``tr``
------

The trace of a matrix A is the sum of the diagonal entries.  This can
be written as::

  t = tr(A)
  t = A[i,i]

The definition of the trace is

.. math::

  \mbox{\texttt{tr}(A)} \leftrightarrow \mathrm{tr} \mathbf{A} = A_{ii} = \sum_{i=0}^{n-1} A_{ii}.


``dot``
-------

The dot product of two tensors a and b can be written::

  # General tensors
  f = dot(a, b)

  # Vectors a and b
  f = a[i]*b[i]

  # Matrices a and b
  f = as_matrix(a[i,k]*b[k,j], (i,j))

The definition of the dot product of unit vectors is (assuming an
orthonormal basis for a Euclidean space):

.. math::

  \mathbf{i}_i \cdot \mathbf{i}_j = \delta_{ij}

where :math:`\delta_{ij}` is the Kronecker delta as explained earlier.
The dot product of higher order tensors follow from this, as illustrated
with the following examples.

An example with two vectors

.. math::

   \mathbf{v} \cdot \mathbf{u} = (v_i \mathbf{i}_i) \cdot (u_j \mathbf{i}_j)
        = v_i u_j (\mathbf{i}_i \cdot \mathbf{i}_j) = v_i u_j \delta_{ij} = v_i u_i

An example with a tensor of rank two

.. math::

  \mathbf{A} \cdot \mathbf{B}
  &= (A_{ij} \mathbf{i}_i \mathbf{i}_j) \cdot (B_{kl} \mathbf{i}_k \mathbf{i}_l) \\
  &= (A_{ij}B_{kl}) \mathbf{i}_i(\mathbf{i}_j \cdot \mathbf{i}_k) \mathbf{i}_l \\
  &= (A_{ij}B_{kl}\delta_{jk}) \mathbf{i}_i \mathbf{i}_l \\
  &= A_{ik}B_{kl} \mathbf{i}_i \mathbf{i}_l.

This is the same as to matrix-matrix multiplication.

An example with a vector and a tensor of rank two

.. math::

   \mathbf{v} \cdot \mathbf{A}
   &= (v_j \mathbf{i}_j) \cdot (A_{kl} \mathbf{i}_k \mathbf{i}_l) \\
   &= (v_j A_{kl}) (\mathbf{i}_j \cdot \mathbf{i}_k) \mathbf{i}_l \\
   &= (v_j A_{kl}\delta_{jk}) \mathbf{i}_l \\
   &= v_k A_{kl} \mathbf{i}_l

This is the same as to vector-matrix multiplication.

% TODO: Is 'contraction' or 'axis' obvious and exactly correctly used?
%       Get a better formulation from somewhere?
This generalizes to tensors of arbitrary rank:
%The dot product is a contraction over the
The dot product applies to the last axis of a and the first axis of b.
The tensor rank of the product is rank(a)+rank(b)-2.

%and generalized to tensors of arbitrary rank you get this crazy expression:
%\begin{align}
%(A_\iota^a \ii_{\iota^a_0}\otimes\cdots\otimes\ii_{\iota^a_{r-1}})
%\cdot
%(B_\iota^b \ii_{\iota^b_0}\otimes\cdots\otimes\ii_{\iota^b_{r-1}})
%\\
%=
%(A_\iota^a B_\iota^b)
%(\ii_{\iota^a_{r-1}} \cdot \ii_{\iota^b_0})
%\ii_{\iota^a_0}\otimes\cdots\otimes\ii_{\iota^a_{r-2}}
%\otimes
%\ii_{\iota^b_1}\otimes\cdots\otimes\ii_{\iota^b_{r-1}}
%\\
%=
%(A_\iota^a B_\iota^b \delta_{\iota^a_{r-1} \iota^b_0})
%\ii_{\iota^a_0}\otimes\cdots\otimes\ii_{\iota^a_{r-2}}
%\otimes
%\ii_{\iota^b_1}\otimes\cdots\otimes\ii_{\iota^b_{r-1}}
%\end{align}
%\begin{align}
%\ii_i\cdot\ii_j = \delta_{ij}
%\end{align}

``inner``
---------

The inner product is a contraction over all axes of a and b, that is
the sum of all componentwise products.  The operands must have the exact
same dimensions.  For two vectors it is equivalent to the dot product.

If $\mathbf{A}$ and $\mathbf{B}$ are rank two tensors and $\mathcal{C}$
and $\mathcal{D}$ are rank 3 tensors their inner products are

.. math::
   \mathbf{A} : \mathbf{B}   &= A_{ij} B_{ij}
   \\
   \mathcal{c} : \mathcal{D} &= C_{ijk} D_{ijk}

Using UFL notation, the following pairs of declarations are equivalent::

  # Vectors
  f = inner(a, b)
  f = v[i]*b[i]

  # Matrices
  f = inner(A, B)
  f = A[i,j]*B[i,j]

  # Rank 3 tensors
  f = inner(C, D)
  f = C[i,j,k]*D[i,j,k]


``outer``
---------

The outer product of two tensors a and b can be written::

  A = outer(a, b)

The general definition of the outer product of two tensors
:math:`\mathcal{C}` of rank :math:`r` and :math:`\mathcal{D}` of rank
:math:`s` is

.. math::

   \Cc \otimes \Dc
    =
    C_{\iota^a_0\ldots\iota^a_{r-1}} D_{\iota^b_0\ldots\iota^b_{s-1}}
    \ii_{\iota^a_0}\otimes\cdots\otimes\ii_{\iota^a_{r-2}}
    \otimes
    \ii_{\iota^b_1}\otimes\cdots\otimes\ii_{\iota^b_{s-1}}

Some examples with vectors and matrices are easier to understand

.. math::

   \vv \otimes \uu = v_i u_j \ii_i\ii_j, \\
   \vv \otimes \BB = v_i B_{kl} \ii_i\ii_k\ii_l, \\
   \AA \otimes \BB = A_{ij} B_{kl} \ii_i\ii_j\ii_k\ii_l .

The outer product of vectors is often written simply as

.. math::

   \vv \otimes \uu = \vv\uu,

which is what we've done with $\ii_i\ii_j$ above.

The rank of the outer product is the sum of the ranks of the operands.

``cross``
---------

The operator ``cross`` accepts as arguments two logically vector-valued
expressions and returns a vector which is the cross product (vector
product) of the two vectors:

.. math::
   \mbox{\texttt{cross(v, w)}} \leftrightarrow \vv \times \ww
  = (v_1 w_2 - v_2 w_1, v_2 w_0 - v_0 w_2, v_0 w_1 - v_1 w_0).

Note that this operator is only defined for vectors of length three.

``det``
-------

The determinant of a matrix A can be written::

  d = det(A)

``dev``
-------

The deviatoric part of matrix A can be written::

  B = dev(A)

``sym``
-------

The symmetric part of A can be written::

  B = sym(A)

The definition is

.. math::

  \mathop{sym} \AA = \frac 1 2 (\AA + \AA^T)

``skew``
--------

The skew symmetric part of A can be written::

  B = skew(A)

The definition is

.. math::

   \mathop{skew} \AA = \frac 1 2 (\AA - \AA^T)

``cofac``
---------

The cofactor of a matrix A can be written::

  B = cofac(A)

The definition is

.. math::

   \mathop{cofac} \AA = \mathop{det}(\AA) \AA^{-1}

The implementation of this is currently rather crude, with a hardcoded
symbolic expression for the cofactor.  Therefore, this is limited to 1x1,
2x2 and 3x3 matrices.

``inv``
-------

The inverse of matrix A can be written::

  Ainv = inv(A)

The implementation of this is currently rather crude, with a hardcoded
symbolic expression for the inverse.  Therefore, this is limited to 1x1,
2x2 and 3x3 matrices.


Differential Operators
======================

Three different kinds of derivatives are currently supported: spatial
derivatives, derivatives w.r.t. user defined variables, and derivatives
of a form or functional w.r.t. a function.


Basic spatial derivatives
-------------------------

Spatial derivatives hold a special place in partial differential equations
from physics and there are several ways to express those. The basic way is::

  # Derivative w.r.t. x_2
  f = Dx(v, 2)
  f = v.dx(2)
  # Derivative w.r.t. x_i
  g = Dx(v, i)
  g = v.dx(i)

% TODO: Document this below
%# Derivative w.r.t. x
%x = cell.x
%h = diff(v, x)

If ``v`` is a scalar expression, ``f`` here is the scalar derivative of
``v`` w.r.t. spatial direction z.  If ``v`` has no free-indices, ``g``
is the scalar derivative w.r.t. spatial direction :math:`x_i`, and ``g``
has the free-index ``i``.  Written as formulas, this can be expressed
compactly using the :math:`v_{,i}` notation:

.. math::

   f = \frac{\partial v}{\partial x_2} = v_{,2}, \\
   g = \frac{\partial v}{\partial x_i} = v_{,i}.

Note the resemblance of :math:`v_{,i}` and :math:`v.dx(i)`.

If the expression to be differentiated w.r.t. :math:`x_i` has ``i``
as a free-index, implicit summation is implied::

  # Sum of derivatives w.r.t. x_i for all i
  g = Dx(v[i], i)
  g = v[i].dx(i)

Here ``g`` will represent the sum of derivatives
w.r.t. :math:`x_i` for all ``i``, that is

.. math::

   g = \sum_i \frac{\partial v}{\partial x_i} = v_{i,i}.

Note the compact index notation :math:`v_{i,i}` with implicit summation.


Compound spatial derivatives
----------------------------

UFL implements several common differential operators.  The notation is
simple and their names should be self explaining::

  Df = grad(f)
  df = div(f)
  cf = curl(v)
  rf = rot(f)

The operand ``f`` can have no free indices.

%NB! Although their general meaning is well defined, pay
%attention to their exact definition here, because there
%are different traditional ways to interpret them.
%For example, \ufl{} transposes the gradient of a vector
%compared to the old \ffc{} notation. The definition used in
%\ufl{} is consistent with the traditions of fluid mechanics,
%allowing for example writing the convection term
%$\ww\cdot\nabla\uu\cdot\vv$ in a natural fashion like
%\begin{code}
%dot( dot(w, grad(u)), v )
%\end{code}

%The definition of these operators follow from
%the vector of partial derivatives
%\[
%\nabla \equiv \frac{\partial}{\partial x_k} \mathbf{i}_k
%    = \sum_{k=0}^{d-1} \frac{\partial}{\partial x_k} \mathbf{i}_k,
%\]
%and the definition of the dot product, outer product,
%and cross product from the previous section.

Gradient
--------

The gradient of a scalar $u$ is defined as
\begin{align}
\mathop{grad}(u) \equiv \nabla u =
    \sum_{k=0}^{d-1} \frac{\partial u}{\partial x_k} \ii_k,
\end{align}
which is a vector of all spatial partial derivatives of $u$.

The gradient of a vector $\vv$ is defined as
\begin{align}
\mathop{grad}(\vv) \equiv \nabla \vv
    = \frac{\partial v_i}{\partial x_j} \ii_i \ii_j,
\end{align}
which written componentwise is
\begin{align}
\AA = \nabla\vv, \qquad A_{ij} = v_{i,j}
\end{align}
In general for a tensor $\AA$ of rank $r$ the definition is
\begin{align}
\mathop{grad}(\AA) \equiv \nabla \AA
    = (\frac{\partial}{\partial x_i}) (A_\iota\ii_{\iota_0}\otimes\cdots\otimes\ii_{\iota_{r-1}}) \otimes \ii_i
    = \frac{\partial A_\iota}{\partial x_i} \ii_{\iota_0}\otimes\cdots\otimes\ii_{\iota_{r-1}}\otimes\ii_i,
\end{align}
where $\iota$ is a multiindex of length $r$.

In \ufl{}, the following pairs of declarations are equivalent:
\begin{code}
Dfi = grad(f)[i]
Dfi = f.dx(i)

Dvi = grad(v)[i, j]
Dvi = v[i].dx(j)

DAi = grad(A)[..., i]
DAi = A.dx(i)
\end{code}
for a scalar expression \texttt{f}, a vector expression \texttt{v},
and a tensor expression \texttt{A} of arbitrary rank.

\subsection{Divergence}
%div(f) = f[i,...].dx(i)

The divergence of any nonscalar (vector or tensor) expression $\AA$ is
defined as the contraction of the partial derivative over the last
axis of the expression.

TODO: Detailed examples like for gradient.

In \ufl{}, the following declarations are equivalent:
\begin{code}
dv = div(v)
dv = v[i].dx(i)

dA = div(A)
dA = A[..., i].dx(i)
\end{code}
for a vector expression v and a tensor expression A.

\subsection{Curl and rot}
\index{rotation}
\index{curl}
\index{\texttt{curl}}
\index{\texttt{rot}}

The operator \texttt{curl} accepts as argument a
vector-valued expression and returns its curl:
\begin{equation}
  \mbox{\texttt{curl(v)}} \leftrightarrow \mathrm{curl} \, \vv = \nabla \times \vv
  = (\frac{\partial v_2}{\partial x_1} - \frac{\partial v_1}{\partial x_2},
  \frac{\partial v_0}{\partial x_2} - \frac{\partial v_2}{\partial x_0},
  \frac{\partial v_1}{\partial x_0} - \frac{\partial v_0}{\partial x_1}).
\end{equation}
Note that this operator is only defined for vectors of length three.

%Alternatively, the name \texttt{rot} can be used for this operator.
%TODO: Define rot.

\subsection{Variable derivatives}

\ufl{} also supports differentiation with respect
to user defined variables. A user defined variable
can be any\footnote{TODO: There are probably some things that don't make sense.}
expression that is defined as a variable.

The notation is illustrated here:
\begin{code}
# Define some arbitrary expression
u = Function(element)
w = sin(u**2)

# Annotate expression w as a variable that can be used in diff
w = variable(w)

# This expression is a function of w
F = I + diff(u, x)

# The derivative of expression f w.r.t. the variable w
df = diff(f, w)
\end{code}
Note that the variable \texttt{w} still represents the same expression.

This can be useful for example to implement
material laws in hyperelasticity where the stress
tensor is derived from a Helmholtz strain energy function.

Currently, \ufl{} does not implement time in any particular way,
but differentiation w.r.t. time can be done without this support
through the use of a constant variable t:
\begin{code}
t = variable(Constant(cell))
f = sin(x[0])**2 * cos(t)
dfdt = diff(f, t)
\end{code}

\subsection{Functional derivatives}
The third and final kind of derivatives are derivatives
of functionals or forms w.r.t. to a \texttt{Function}.
This is described in more detail in section \ref{subsec:AD}
about form transformations.

%------------------------------------------------------------------------------
\section{DG operators}
\index{DG operators}
\index{discontinuous Galerkin}
\index{jump}
\index{avg}
\index{restriction}

\ufl{} provides operators for implementation of discontinuous Galerkin
methods. These include the evaluation of the jump and average
of a function (or in general an expression) over the interior facets
(edges or faces) of a mesh.

\subsection{Restriction: \texttt{v('+')} and \texttt{v('-')}}

When integrating over interior facets (\texttt{*dS}), one may restrict
expressions to the positive or negative side of the facet:
\begin{code}
element = FiniteElement("Discontinuous Lagrange",
                        "tetrahedron", 0)

v = TestFunction(element)
u = TrialFunction(element)

f = Function(element)

a = f('+')*dot(grad(v)('+'), grad(u)('-'))*dS
\end{code}

Restriction may be applied to functions of any finite element space
but will only have effect when applied to expressions that are
discontinuous across facets.

\subsection{Jump: \texttt{jump(v)}}

The operator \texttt{jump} may be used to express the jump of a
function across a common facet of two cells. Two versions of the
\texttt{jump} operator are provided.

If called with only one argument, then the \texttt{jump} operator
evaluates to the difference between the restrictions of the given
expression on the positive and negative sides of the facet:
\begin{equation}
  \mbox{\texttt{jump(v)}}
  \leftrightarrow
  \llbracket v \rrbracket = v^+ - v^-.
\end{equation}
If the expression \texttt{v} is scalar, then \texttt{jump(v)} will
also be scalar, and if \texttt{v} is vector-valued, then \texttt{jump(v)}
will also be vector-valued.

If called with two arguments, \texttt{jump(v, n)} evaluates to the
jump in \texttt{v} weighted by \texttt{n}. Typically, \texttt{n} will
be chosen to represent the unit outward normal of the facet (as seen
from each of the two neighboring cells). If \texttt{v} is scalar, then
\texttt{jump(v, n)} is given by
\begin{equation}
  \mbox{\texttt{jump(v, n)}}
  \leftrightarrow
  \llbracket v \rrbracket_n = v^+ n^+ + v^- n^-.
\end{equation}
If \texttt{v} is vector-valued, then \texttt{jump(v, n)} is given by
\begin{equation}
  \mbox{\texttt{jump(v, n)}}
  \leftrightarrow
  \llbracket v \rrbracket_n = v^+ \cdot n^+ + v^- \cdot n^-.
\end{equation}
Thus, if the expression \texttt{v} is scalar, then \texttt{jump(v, n)} will
be vector-valued, and if \texttt{v} is vector-valued, then
\texttt{jump(v, n)} will be scalar.

\subsection{Average: \texttt{avg(v)}}

The operator \texttt{avg} may be used to express the average
of an expression across a common facet of two cells:
\begin{equation}
  \mbox{avg(v)}
  \leftrightarrow
  \langle v \rangle = \frac{1}{2} (v^+ + v^-).
\end{equation}
The expression \texttt{avg(v)} has the same value shape as the expression \texttt{v}.

%------------------------------------------------------------------------------
\section{Conditional Operators}
\index{conditional operators}

\subsection{Conditional}
\ufl{} has limited support for branching, but for some PDEs it is needed.
The expression \texttt{c} in
\begin{code}
c = conditional(condition, true_value, false_value)
\end{code}
evaluates to \texttt{true\_value} at run-time if \texttt{condition}
evaluates to true, or to \texttt{false\_value} otherwise.

This corresponds to the C++ syntax \texttt{(condition ? true\_value: false\_value)},
or the Python syntax \texttt{(true\_value if condition else false\_value)},

\subsection{Conditions}
\begin{itemize}
\item \texttt{eq(a, b)} represents the condition that a == b
\item \texttt{ne(a, b)} represents the condition that a != b
\item \texttt{le(a, b)} represents the condition that a <= b
\item \texttt{ge(a, b)} represents the condition that a >= b
\item \texttt{lt(a, b)} represents the condition that a <  b
\item \texttt{gt(a, b)} represents the condition that a >  b
\end{itemize}

TODO: This is rather limited, probably need the operations
"and" and "or" as well, the syntax will be rather convoluted...
Can we improve? Low priority though.

\begin{advancedenv}
Because of details in the way Python behaves, we cannot overload
the builtin comparison operators for this purpose, hence these named operators.
\end{advancedenv}

%------------------------------------------------------------------------------
\section{User-defined operators}
\label{sec:user-defined}
\index{user-defined operators}

A user may define new operators, using standard Python syntax. As an
example, consider the strain-rate operator $\epsilon$ of linear elasticity,
defined by
\begin{equation}
  \epsilon(v) = \frac{1}{2} (\nabla v + (\nabla v)^{\top}).
\end{equation}
This operator can be implemented as a function using the Python \texttt{def}
keyword:
\begin{code}
def epsilon(v):
    return 0.5*(grad(v) + grad(v).T)
\end{code}
Alternatively, using the shorthand \texttt{lambda} notation, the
strain operator may be defined as follows:
\begin{code}
epsilon = lambda v: 0.5*(grad(v) + grad(v).T)
\end{code}
\index{def}\index{lambda}

%------------------------------------------------------------------------------
\section{Form Transformations}
\index{form transformations}

When you have defined a \texttt{Form}, you can derive new related
forms from it automatically. UFL defines a set of common
form transformations described in this section.

\subsection{Replacing arguments of a Form}
The function \texttt{replace} lets you replace terminal
objects with other values, using a mapping defined
by a Python dict. This can be used for example to
replace a \texttt{Function} with a fixed value for
optimized runtime evaluation.

\begin{code}
f = Function(element)
g = Function(element)
c = Constant(cell)
a = f*g*v*dx
b = replace(a, { f: 3.14, g: c })
\end{code}

The replacement values must have the same basic properties
as the original values, in particular value shape and
free indices.


\subsection{Action of a form on a function}
The action of a bilinear form $a$ is defined as
\[
b(v; w) = a(v, w),
\]
The action of a linear form $L$ is defined as
\[
f(;w) = L(w)
\]
This operation is implemented in UFL simply by replacing the rightmost
basis function (trial function for $a$, test function for $L$)
in a \texttt{Form}, and is used like this:
\begin{code}
L = action(a, w)
f = action(L, w)
\end{code}
To give a concrete example, these declarations are equivalent:
\begin{code}
a = inner(grad(u), grad(v))*dx
L = action(a, w)

a = inner(grad(u), grad(v))*dx
L = inner(grad(w), grad(v))*dx
\end{code}

If a is a rank 2 form used to assemble the matrix A,
L is a rank 1 form that can be used to assemble the vector $b = Ax$ directly.
This can be used to define both the form of a matrix and the form
of its action without code duplication,
and for the action of a Jacobi matrix computed using derivative.

If L is a rank 1 form used to assemble the vector b,
f is a functional that can be used to assemble the scalar
value $f = b\cdot w$ directly. This operation is sometimes
used in, e.g., error control with L being the residual equation
and w being the solution to the dual problem.
(However, the discrete vector for the assembled residual equation
will typically be available, so doing the dot product using
linear algebra would be faster than using this feature.)
FIXME: Is this right?


\subsection{Energy norm of a bilinear Form}
The functional representing the energy norm $|v|_A = v^T A v$ of
a matrix A assembled from a form $a$ can be computed like this
\begin{code}
f = energy_norm(a, w)
\end{code}
which is equivalent to
\begin{code}
f = action(action(a, w), w)
\end{code}


\subsection{Adjoint of a bilinear Form}
The adjoint $a'$ of a bilinear form $a$ is defined as
\[
a'(u,v) = a(v,u).
\]
This operation is implemented in UFL simply by swapping
test and trial functions in a \texttt{Form}, and is used like this:
\begin{code}
aprime = adjoint(a)
\end{code}


\subsection{Linear and bilinear parts of a Form}
Some times it is useful to write an equation on the format
\[
a(v,u) - L(v) = 0.
\]
Before we can assemble the linear equation
\[
A u = b,
\]
we need to extract the forms corresponding to the left hand
side and right hand side. This corresponds
to extracting the bilinear and linear terms of the form
respectively, or the terms that depend on both a test
and a trial function on one side and the terms that
depend on only a test function on the other.

This is easily done in UFL using \texttt{lhs} and \texttt{rhs}:
\begin{code}
b = u*v*dx - f*v*dx
a, L = lhs(b), rhs(b)
\end{code}
Note that \texttt{rhs} multiplies the extracted terms by $-1$,
corresponding to moving them from left to right, so this
is equivalent to
\begin{code}
a = u*v*dx
L = f*v*dx
\end{code}

As a slightly more complicated example, this formulation
\begin{code}
  F = v*(u - w)*dx + k*dot(grad(v), grad(0.5*(w + u)))*dx
  a, L = lhs(F), rhs(F)
\end{code}
is equivalent to
\begin{code}
  a = v*u*dx + k*dot(grad(v), 0.5*grad(u))*dx
  L = v*w*dx - k*dot(grad(v), 0.5*grad(w))*dx
\end{code}

%\subsection{Splitting integrals by polygon degree}
%TODO: Implement "integrals = split_by_degree(a)"


%\subsection{Blocks of a form on mixed element spaces}
%TODO: Implement A, B, C, D = blocks(a)


\subsection{Automatic Functional Differentiation}\label{subsec:AD}

\ufl{} can compute derivatives of functionals
or forms w.r.t. to a \texttt{Function}.
This functionality can be used for example to
linearize your nonlinear residual equation automatically,
or derive a linear system from a functional,
or compute sensitivity vectors w.r.t. some coefficient.

A functional can be differentiated to obtain a linear form,
\[
F(v; w) = \frac{d}{dw} f(;w)
\]
and a linear form
  \footnote{Note that by ``linear form'' we only mean a form that is linear
  in its test function, not in the function you differentiate with respect to.}
can be differentiated to obtain the bilinear form
corresponding to its Jacobi matrix:
\[
J(v, u; w) = \frac{d}{dw} F(v; w).
\]
The UFL code to express this is (for a simple functional $f(w)=\int_\Omega \frac 1 2 w^2\,dx$)
\begin{code}
f = (w**2)/2 * dx
F = derivative(f, w, v)
J = derivative(F, w, u)
\end{code}
which is equivalent to:
\begin{code}
f = (w**2)/2 * dx
F = w*v*dx
J = u*v*dx
\end{code}

Assume in the following examples that:
\begin{code}
v = TestFunction(element)
u = TrialFunction(element)
w = Function(element)
\end{code}

The stiffness matrix can be computed from the functional
$\int_\Omega \nabla w : \nabla w \, dx$,
by the lines
\begin{code}
f = inner(grad(w), grad(w))/2 * dx
F = derivative(f, w, v)
J = derivative(F, w, u)
\end{code}
which is equivalent to:
\begin{code}
f = inner(grad(w), grad(w))/2 * dx
F = inner(grad(w), grad(v)) * dx
J = inner(grad(u), grad(v)) * dx
\end{code}
Note that here the basis functions are provided explicitly,
which is some times necessary, e.g., if part of the form
is linearlized manually like in
(\emph{TODO: An example that makes sense would be nicer, this is just a random form.})
\begin{code}
g = Function(element)
f = inner(grad(w), grad(w))*dx
F = derivative(f, w, v) + dot(w-g,v)*dx
J = derivative(F, w, u)
\end{code}

Derivatives can also be computed w.r.t. functions in mixed spaces.
Consider this example, an implementation of the harmonic map
equations using automatic differentiation.
\begin{code}
X = VectorElement("Lagrange", cell, 1)
Y = FiniteElement("Lagrange", cell, 1)

x = Function(X)
y = Function(Y)

L = inner(grad(x), grad(x))*dx + dot(x,x)*y*dx

F = derivative(L, (x,y))
J = derivative(F, (x,y))
\end{code}
Here \texttt{L} is defined as a functional with two coefficient functions
\texttt{x} and \texttt{y} from separate finite element spaces.
However, \texttt{F} and \texttt{J} become
linear and bilinear forms respectively with basis functions
defined on the mixed finite element
\begin{code}
M = X + Y
\end{code}
There is a subtle difference between defining \texttt{x} and \texttt{y}
separately and this alternative implementation
(reusing the elements \texttt{X},\texttt{Y},\texttt{M}):
\begin{code}
u = Function(M)
x, y = split(u)

L = inner(grad(x), grad(x))*dx + dot(x,x)*y*dx

F = derivative(L, u)
J = derivative(F, u)
\end{code}
The difference is that the forms here have \emph{one}
coefficient function \texttt{u} in the mixed space, and the forms above
have \emph{two} coefficient functions \texttt{x} and \texttt{y}.


TODO: Move this to implementation part?
If you wonder how this is all done, a brief explanation follows.
Recall that a \texttt{Function} represents a
sum of unknown coefficients multiplied with unknown
basis functions in some finite element space.
\begin{align}
w(x) = \sum_k w_k \phi_k(x)
\end{align}
Also recall that a \texttt{BasisFunction} represents any
(unknown) basis function in some finite element space.
\begin{align}
v(x) = \phi_k(x), \qquad \phi_k \in V_h .
\end{align}
A form $L(v; w)$ implemented in \ufl{} is intended for discretization like
\begin{align}
b_i = L(\phi_i; \sum_k w_k \phi_k), \qquad \forall \phi_i \in V_h .
\end{align}
The Jacobi matrix $A_{ij}$ of this vector can be obtained by
differentiation of $b_i$ w.r.t. $w_j$, which can be written
\begin{align}
A_{ij} = \frac{d b_i}{d w_j} = a(\phi_i, \phi_j; \sum_k w_k \phi_k), \qquad \forall \phi_i \in V_h, \quad \forall \phi_j \in V_h ,
\end{align}
for some form $a$. In \ufl{}, the form $a$ can be obtained by differentiating $L$.
To manage this, we note that as long as the domain $\Omega$ is
independent of $w_j$, $\int_\Omega$ commutes with $\frac{d}{d w_j}$,
and we can differentiate the integrand expression instead, e.g.,
\begin{align}
L(v; w) = \int_\Omega I_c(v; w) \, dx + \int_{\partial\Omega} I_e(v; w) \, ds, \\
\frac{d}{d w_j} L(v; w) = \int_\Omega \frac{d I_c}{d w_j} \, dx + \int_{\partial\Omega} \frac{d I_e}{d w_j} \, ds.
\end{align}
In addition, we need that
\begin{align}
\frac{d w}{d w_j} = \phi_j, \qquad \forall \phi_j \in V_h ,
\end{align}
which in \ufl{} can be represented as
\begin{align}
w &= \text{\texttt{Function(element)}}, \\
v &= \text{\texttt{BasisFunction(element)}}, \\
\frac{dw}{d w_j} &= v,
\end{align}
since $w$ represents the sum and $v$ represents any and all basis functions in $V_h$.

Other operators have well defined derivatives, and by repeatedly applying
the chain rule we can differentiate the integrand automatically.

\emph{The notation here has potential for improvement, feel free
to ask if something is unclear, or suggest improvements.}



\subsection{Combining form transformations}
Form transformations can be combined freely.
Note that to do this, derivatives are usually
be evaluated before applying e.g. the action
of a form, because \texttt{derivative} changes
the arity of the form.

\begin{code}
element = FiniteElement("CG", cell, 1)
w = Function(element)
f = w**4/4*dx(0) + inner(grad(w), grad(w))*dx(1)
F = derivative(f, w)
J = derivative(F, w)
Ja = action(J, w)
Jp = adjoint(J)
Jpa = action(Jp, w)
g = Function(element)
Jnorm = energy_norm(J, g)
\end{code}

TODO: Find some more examples, e.g. from error control!

%------------------------------------------------------------------------------
\section{Tuple Notation}
\index{tuple notation}

In addition to the standard integrand notation described above, UFL
supports a simplified \emph{tuple notation} by which $L^2$ inner
products may be expressed as tuples. Consider for example the
following bilinear form as part of a variational problem for a
reaction--diffusion problem:
\begin{displaymath}
  \begin{split}
  a(v, u)
  &= \int_{\Omega} \nabla v \cdot \nabla u + v u \dx \\
  &= (\nabla v, \nabla u) + (v, u)
  \end{split}
\end{displaymath}
In standard UFL notation, this bilinear form may be expressed as
\begin{code}
a = inner(grad(v), grad(u))*dx + v*u*dx
\end{code}
In tuple notation, this may alternatively be expressed as
\begin{code}
a = (grad(v), grad(u)) + (v, u)
\end{code}

In general, a form may be expressed as a sum of tuples or triples of the form
\begin{code}
(v, w)
(v, w, dm)
\end{code}
where \texttt{v} and \texttt{w} are expressions of matching rank (so that
\texttt{inner(v, w)} makes sense), and \texttt{dm} is a measure. If the
measure is left out, it is assumed that it is \texttt{dx}.

The following example illustrates how to express a form containing
integrals over subdomains and facets:
\begin{code}
a = (grad(v), grad(u)) + (v, b*grad(u), dx(2))
  + (v, u, ds) + (jump(v), jump(u), dS)
\end{code}

The following caveats should be noted:
\begin{itemize}
\item
  The only operation allowed on a tuple is addition. In particular,
  tuples may not subtracted. Thus,
  \texttt{a = (grad(v), grad(u)) - (v, u)} must be expressed as
  \texttt{a = (grad(v), grad(u)) + (-v, u)}.
\item
  Tuple notation may not be mixed with standard UFL integrand
  notation. Thus, \texttt{a = (grad(v), grad(u)) + inner(v, u)*dx} is not
  valid.
\end{itemize}

\begin{advancedenv}
Tuple notation is strictly speaking not a part of the form
language, but tuples may be converted to UFL forms using the function
\texttt{tuple2form} available from the module \texttt{ufl.algorithms}. This
is normally handled automatically by form compilers, but the
\texttt{tuple2form} utility may useful when working with UFL from a
Python script. Automatic conversion is also carried out by UFL form
operators such as \texttt{lhs} and \texttt{rhs}.
\end{advancedenv}

%------------------------------------------------------------------------------
\section{Form Files}
\index{form files}
\index{ufl files}

\ufl{} forms and elements can be collected in a \emph{form file}
with the extension \texttt{.ufl}. Form compilers will typically
execute this file with the global \ufl{} namespace available,
and extract forms and elements that are defined after execution.
The compilers do not compile all forms and elements that are
defined in file, but only those that are \emph{exported}.
A finite element with the variable name \texttt{element} is
exported by default, as are forms with the names \texttt{M},
\texttt{L}, and \texttt{a}. The default form names are intended
for a functional, linear form, and bilinear form respectively.

To export multiple forms and elements or use other names,
an explicit list with the forms and elements to export
can be defined. Simply write
\begin{code}
elements = [V, P, TH]
forms = [a, L, F, J, L2, H1]
\end{code}
at the end of the file to export the elements and forms
held by these variables.



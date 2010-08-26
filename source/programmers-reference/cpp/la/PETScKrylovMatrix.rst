.. Documentation for the header file dolfin/la/PETScKrylovMatrix.h

.. _programmers_reference_cpp_la_Mesh:

PETScKrylovMatrix.h
===================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: PETScKrylovMatrix

    *Parent class*
    
        * :cpp:class:`PETScBaseMatrix`
        
        This class represents a matrix-free matrix of dimension m x m.
        It is a simple wrapper for a PETSc shell matrix. The interface
        is intentionally simple. For advanced usage, access the PETSc
        Mat pointer using the function mat() and use the standard PETSc
        interface.
        
        The class PETScKrylovMatrix enables the use of Krylov subspace
        methods for linear systems Ax = b, without having to explicitly
        store the matrix A. All that is needed is that the user-defined
        PETScKrylovMatrix implements multiplication with vectors. Note that
        the multiplication operator needs to be defined in terms of
        PETSc data structures (Vec), since it will be called from PETSc.

    .. cpp:function:: PETScKrylovMatrix()
    
        Constructor

    .. cpp:function:: PETScKrylovMatrix(const PETScVector& x, const PETScVector& y)
    
        Create a virtual matrix matching the given vectors

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: virtual void mult(const PETScVector& x, PETScVector& y) const = 0
    
        Compute product y = Ax

    .. cpp:function:: virtual ~PETScKrylovMatrix()
    
        Destructor

    .. cpp:function:: void init(const PETScVector& x, const PETScVector& y)
    
        Initialize virtual matrix matching the given vectors

    .. cpp:function:: void resize(int M, int N)
    
        Resize virtual matrix, all on one processor


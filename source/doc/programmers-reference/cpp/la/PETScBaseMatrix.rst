.. Documentation for the header file dolfin/la/PETScBaseMatrix.h

.. _programmers_reference_cpp_la_petscbasematrix:

PETScBaseMatrix.h
=================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: PETScMatrixDeleter

.. cpp:class:: PETScBaseMatrix

    *Parent class*
    
        * :cpp:class:`PETScObject,`
        
    This class is a base class for matrices that can be used in
    PETScKrylovSolver.

    .. cpp:function:: PETScBaseMatrix()
    
        Constructor

    .. cpp:function:: PETScBaseMatrix(boost::shared_ptr<Mat> A)
    
        Constructor

    .. cpp:function:: void resize(uint m, uint n) = 0
    
        Resize virtual matrin

    .. cpp:function:: uint size(uint dim) const
    
        Return number of rows (dim = 0) or columns (dim = 1) along dimension dim

    .. cpp:function:: boost::shared_ptr<Mat> mat() const
    
        Return PETSc Mat pointer

    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)

